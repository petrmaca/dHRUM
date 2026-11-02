# dHRUM — project context

R package implementing the **dHRUM distributed hydrological model** (C++ core via Rcpp).
Repo: https://github.com/petrmaca/dHRUM

## Architecture

- `src/`: C++ model core.
  - `dHRUM.cpp` — basin-level class (`dHRUM`), owns `std::vector<single_HMunit>`, basin aggregate `data_HB_1d basinDta`; OpenMP loops.
  - `single_HMunit.cpp` (~4,800 lines) — per-HRU hydrological balance; enum-dispatched component models (soil storages, groundwater reservoirs, interception, surface, fast response, snow, fast/slow division).
  - `data_HB_1d.cpp` — time-series container, calendar, PET.
  - `params.cpp` — parameter values and bounds (huge switch statements over parameter enums).
  - `dHRU_*_R_Rcpp.cpp` — thin Rcpp glue files, one per concern (inputs, params, outputs, pet, gwTypes, soilSTypes, …).
- `inst/include/`: C++ headers (`dHRUM.h`, `single_HMunit.h`, `params.h`, `parStructSels.h`, `numberSel.h`, …).
- `R/`: only auto-generated `RcppExports.R` plus `lai.r`. **There is no hand-written modeling layer** — everything user-facing is an Rcpp wrapper.
- `tests/`: loose runnable scripts, **not testthat** (run with `source()`; `tests/test.r` = main 1-HRU smoke test, `pet_regression.r` self-checks with `stop()`). Dev-notes `tests/HowToAdd*` document the repo's conventions for extending the C++ core (new HM structure, new output variable, new parameter).
- `Calibrations/`: standalone optimization scripts (not package code): Amalie (`SCEoptim` + DEoptim versions), Camel, POH/Tepla, pond, wetland. Templates for any calibration work.

## Key conventions

- NAMESPACE uses `exportPattern("^[[:alpha:]]+")` — new exported R functions need no NAMESPACE change.
- DESCRIPTION Imports only `Rcpp`; add optimizer deps (hydroGOF, RcppDE, SoilHyP) to `Suggests`.
- Build with `devtools::install()` or `R CMD INSTALL .` after any C++ change (committed `src/*.o` may be stale). After adding/moving Rcpp `[[Rcpp::export]]` functions, regenerate wrappers with `Rcpp::compileAttributes()`.
- C++17, `-fopenmp`. New C++ features follow the `tests/HowToAdd*` checklists (enum + name-vector in `parStructSels.h`, switches in `data_HB_1d`/`single_HMunit`/`params.cpp`, glue in a matching `dHRU_*_R_Rcpp.cpp` file).
- Model config is sequential: init pointer → set component types → set PT(L) inputs → set params → `calcPetToAllHrus` → `dHRUMrun`. Key output column: `TOTR` (total runoff); 41 output vars incl. `DIRS` (total direct runoff storage, from the master merge).
- Scripted runs fail in RStudio with "figure margins too large" if the plot pane is small (base-R plots everywhere) — verify numerically instead of sourcing plot-heavy scripts.

## Known pitfalls

- **Do not calibrate parameter `"SMAXpdm"`** — the string→enum map in `src/dHRU_params_R_Rcpp.cpp` (lines ~58, 377) silently maps it to `B_SOIL`.
- `allParNames` has a duplicated `"CAN_ST"`.
- **Use `nthreads = 1` for calibration/reproducibility** — `gatherTsFromHrus` has a known OpenMP race making basin sums nondeterministic (fix planned in REFACTORING_PLAN Stage 1.2).
- OS-temp note: use the assistant scratch dir, not the user's R `tempdir()`.

## Refactoring context (REFACTORING_PLAN.md)

Staged behavior-preserving C++ refactor (Stages 0–5: golden-output harness, copy-ctor bug fix `init_CanS`, OMP race fix, enum-indexed storage, params table, calendar/PET extraction). Constraints:
- **No change to the R-facing API**; `Calibrations/` scripts must keep working.
- Behavior-preserving apart from Stage-1 bug fixes; outputs bit-for-bit identical otherwise.
- Build golden-output harness BEFORE touching C++.

## Optimization (`optim` branch) work — current state and plan

There is currently **no C++ optimization/objective code** anywhere in the repo (verified across `src/`, `inst/include/`, and all branches) — all calibration runs happen in R via `Calibrations/` scripts.
Goal on this branch: add a package-level calibration layer.

### Stage A — R-only calibration layer (do first, no C++ changes)

- New file `R/calibration.r`: objective-closure factory + `calibrate_dhurm()` wrapper.
- Pulls initial values/bounds from `getCurdHRUpars()` ($Cur_par/$Up_bound/$Low_bound).
- Each evaluation: `setParamsToAlldHrus` → `calcHBInAlldHrus` → `gatherHBdata` → `getOutput` → extract `TOTR` → compute criterion.
- Switchable backend (`optim` L-BFGS-B / `RcppDE::DEoptim` / `SoilHyP::SCEoptim`); objectives NSE/KGE/SSE/FDC-composite with warmup and optional log transform.
- Precedent: `Calibrations/Amalie/dHRUM_SCEoptimBP.r` (best structure; note the `sse()` bug in `dHRUM_DeoptimBP.r` — ignores its `myPar` arg; Camel script does it correctly). Add `tests/test_calibration.r` in the repo's script style.
- Survives REFACTORING_PLAN because it only touches the R-facing API.

### Stage B — C++ optimization layer (preferred design)

Wait until the REFACTORING_PLAN stages land (avoid conflicts with `single_HMunit`/`dHRUM`/`params` rewrites). Then implement the optimization machinery as a **self-contained C++ library layer on top of the model core, independent of R** — usable from C++ main(), Python (pybind11 later), or the existing Rcpp glue. R bindings are thin wrappers only.

**New files** (follow the `tests/HowToAdd*` conventions):

- `inst/include/objective.h` + `src/objective.cpp` — objective functions and their interface.
- `inst/include/optimizer.h` + `src/optimizer.cpp` — pluggable optimizer implementations.
- `src/dHRU_objective_R_Rcpp.cpp` — thin Rcpp glue exposing both to R (one-file-per-concern convention).

**Objective-function interface — the key new abstraction**: allow each calibration to specify **what is observed vs what is simulated, per HRU**, so any time series can be calibrated against any model output:

```cpp
// pairing: which observed series is matched against which model output, for which HRU
struct ObsSimPair {
  unsigned hruIndex;        // HRU (0-based); use basin aggregate via a reserved index or a separate basin flag
  ts_type   simTsType;      // simulated series, e.g. ts_type::TOTR, ts_type::SOIS
  const hdata* obs;         // observed series (or store by owning object); length = sim length
  double    weight;         // weight in a composite objective
};

class ObjectiveBase {
 public:
  virtual ~ObjectiveBase() = default;
  virtual double evaluate(dHRUM& model, const std::vector<ObsSimPair>& pairs,
                          unsigned warmup) const = 0;
};
// Implementations: NseObjective, KgeObjective, SseObjective, FdcObjective, WeightedCompositeObjective.
// Each reads the simulated hdata via model.get_HbDta(pairs[i].simTsType) (or
// model.getSingleHruTsDta-equivalent exposed on dHRUM) and scores against pairs[i].obs
// after skipping `warmup` time steps; optional log/sqrt transform flag.
```

This generalizes from "extract `TOTR`" to **calibration against any selected time series in the hydrological data — and a different series per HRU** (e.g. `TOTR` from HRU_3 against gauge 1, `SOIS` from HRU_7 against soil-moisture observations), not just basin runoff.

**Optimizer engine interface** (enables adding future optimizers and new functions without touching the model):

```cpp
class OptimizerBase {
 public:
  virtual ~OptimizerBase() = default;
  virtual OptimizationResult optimize(ObjectiveBase& obj,
                                      const std::vector<double>& init,
                                      const std::vector<double>& lower,
                                      const std::vector<double>& upper,
                                      const OptimizerOptions& opts) = 0;
};
// Initial implementations: DeOptimizer (SCE-UA / DE as in Calibrations use),
// plus a simple NelderMead or CMA-ES. OptimizationResult holds best params, criterion,
// n_evals, convergence trace.
```

`ObjectiveBase` and `OptimizerBase` are pure interfaces with virtual destructors — **new objective functions and new optimizers are added by subclassing, requiring zero changes to `dHRUM`/`single_HMunit`/`params`.** The engine only calls existing public model methods (`setParamsToAlldHrus`, `calcHbToAllHrus`, `gatherTsFromHrus`, `get_HbDta`).

**Build / dependency rules**:
- `objective.cpp` / `optimizer.cpp` link against `dHRUM.h`/`inst/include` but **must not include Rcpp/R headers**; keep them Rcpp-free so the layer also compiles standalone (e.g. a future `main()` driver or pybind11 module).
- Add `#include "objective.h"` etc. only in the Rcpp glue; `src/Makevars` unchanged (already `-I../inst/include`, CXX17, OpenMP).
- After adding `[[Rcpp::export]]` functions in `dHRU_objective_R_Rcpp.cpp`, run `Rcpp::compileAttributes()`.

**Performance note**: per-iteration cost today is dominated by shipping an N×41 matrix through `getOutput()`. The C++ objective avoids this entirely — it works on the internal valarrays via `get_HbDta()`. If params-setting proves hot, add a `setParams` overload taking `(enum index, value)` pairs (planned anyway in the params-table refactor) to skip the string→enum lookup per iteration.

**Threading caveat**: the OMP race in `gatherTsFromHrus` (REFACTORING_PLAN Stage 1.2) must be fixed before running multi-HRU objectives with `nthreads > 1`; until then C++ objectives are deterministic only at `nthreads = 1`.

**Verification**: reuse the golden-output harness (Stage 0) — a C++ objective must return bit-for-bit identical simulated series vs the R path (`getOutput()` → column extract) for the same parameter vector, and the C++ NSE must match `hydroGOF::NSE` within tight tolerance.

### Lumped vs distributed calibration — single logical interface

Lumped (1 HRU) is the degenerate case of distributed (N HRUs): `lumped == ScalarObjective({{0, TOTR, obs, 1}}, crit)`. The optimizer never special-cases — it only calls `obj.evaluate(model, warmup)`.

#### Step 1 — extend `dHRUM` (additive only, no change to existing public API)

```cpp
// additions to inst/include/dHRUM.h — hydrology code untouched
class dHRUM {
 public:
  // per-HRU simulated series (thin wrapper over dHruVec[i].getSingleHruTsDta)
  hdata  getHruTsDta(const ts_type& t, unsigned hruIndex) const;

  // enum-indexed param vector get/set — no string lookup per iteration.
  // Complements (does not replace) setParamsToAlldHrus(pair-vector); the params-table
  // refactor (Stage 4) already wants this, so it aligns rather than conflicts.
  std::vector<numberSel> getParamVec(unsigned hru) const;                     // exists: get_param_vec
  void setParamVec   (unsigned hru, const std::vector<numberSel>& v);         // NEW
  void setParamVecAll(const std::vector<numberSel>& v);                       // NEW, lumped fast path

  unsigned getNumHrus() const { return dimHM; }
  unsigned getNumPars(unsigned hru) const;                                    // exists: get_singleHRUnumPars
};
```

The only genuinely new model-side capability is exposing per-HRU time series — `single_HMunit::getSingleHruTsDta` already exists (src/single_HMunit.cpp:3494) but isn't reachable from `dHRUM`.

#### Step 2 — objective layer (Rcpp-free C++, links only against `inst/include`)

```cpp
// observed target = pair (observed series, simulated series chosen per HRU)
struct CalTarget {
  unsigned  hruIndex;   // 0..dimHM-1; lumped always 0. Basin aggregate: hruIndex==dimHM or isBasin flag
  ts_type   simType;    // simulated series, e.g. ts_type::TOTR / SOIS
  const hdata* obs;     // caller-owned observed series, same length as sim
  double    weight = 1;
};

class Criterion { // single-pair criterion after `warmup` skipped: NSE, KGE, SSE, FDC, logNSE…
 public:
  virtual ~Criterion() = default;
  virtual double value(const hdata& sim, const hdata& obs, unsigned warm) const = 0;
};

class Objective {
 public:
  virtual ~Objective() = default;
  virtual double evaluate(dHRUM& model, unsigned warmup) const = 0;
};

// ONE concrete objective covers both lumped and distributed:
class ScalarObjective final : public Objective {
 public:
  ScalarObjective(std::vector<CalTarget> targets,       // 1 target = lumped; M = distributed
                  std::unique_ptr<Criterion> crit,
                  bool logTransform = false);
  double evaluate(dHRUM& m, unsigned warm) const override {
    // distribute targets over OpenMP threads when safe (see "OpenMP use" below):
    // #pragma omp parallel for reduction(+:acc) schedule(static) if(m.getNumHrus()>1)
    //   sim = m.getHruTsDta(t.simType, t.hruIndex); optional log;
    //   acc += t.weight * crit->value(sim, *t.obs, warm);
    // return acc / totalWeight;
  }
};
```

#### Step 3 — `ParamMapper` (where lumped/distributed really lives) + optimizer

```cpp
// maps flattened optimizer vector -> per-HRU parameter slots
class ParamMapper {
 public:
  virtual ~ParamMapper() = default;
  virtual void   apply(dHRUM& m, const std::vector<numberSel>& flat) = 0;
  virtual size_t dim() const = 0;                                      // optimizer's dimension
};

class IdentityMapper : public ParamMapper { /* distributed: per-HRU params, dim = sum_hru nPars(hru) */ };
class BroadcastMapper : public ParamMapper { /* lumped & spatially-uniform distributed: 1 param set -> all HRUs */ };

class Optimizer {
 public:
  virtual ~Optimizer() = default;
  virtual OptResult minimize(Objective& obj, dHRUM& model,
                    ParamMapper& mapper,
                    const std::vector<numberSel>& init,
                    const std::vector<numberSel>& lo,
                    const std::vector<numberSel>& hi,
                    unsigned warmup, const OptParams& ctrl) = 0;
};
// Implementations: DeOptimizer, SceUaOptimizer, NelderMeadOptimizer… (subclass to extend, zero model changes)
```

**Lumped** = `ScalarObjective({{0,TOTR,obs,1}}, crit)` + `BroadcastMapper` (dim = nPars). **Distributed** = `ScalarObjective(many CalTargets, crit)` + `IdentityMapper` (dim = Σ nPars). Fine-grained mixtures (some shared, some per-HRU) = a third mapper subclass; the optimizer never changes.

#### OpenMP use

- The model already parallelizes per-HRU work inside `setParamsToAlldHrus`, `loadPT*DatToAllHrus`, `calcHbToAllHrus`, `gatherTsFromHrus` with `#pragma omp parallel for num_threads(num_threads)`, `num_threads` set via `initdHruModel(..., nthreads)`.
- **Objective scoring** (`ScalarObjective::evaluate`, `Criterion::value`) is parallel-safe: targets share no state, use `#pragma omp parallel for reduction(+:acc) schedule(static)` over `CalTarget`s — benefits multi-target distributed objectives.
- **`setParamVecAll`** can relaunch the same pattern in the new setter (loop over HRUs is embarrassingly parallel).
- **Optimizer population evaluation** is the highest-value parallelism: run each trial's `setParamVec` + `calcHbToAllHrus` on a *cloned* `dHRUM` per thread (deep copy exists via copy-ctor; needed since model state is mutable). Requires thread-safe loggers; evaluation = private clones, objective = reduction.
- **Caveat — do not nest**: inner `#pragma omp parallel for` inside an outer one is invalid unless `omp_set_nested` + explicit nested teams; if the optimizer parallelizes over trials first, set the model's `num_threads` to 1 for those worker instances (an argument for a per-instance `set_num_treads(1)` before cloning).
- **`gatherTsFromHrus` race** (REFACTORING_PLAN Stage 1.2) still gates *any* multi-HRU objective relying on basin sums at `nthreads>1`. Workarounds: `nthreads==1`, or fix the race first (private per-thread accumulator + reduction like above). Lumped (`dimHM==1`) is unaffected but still passes through the same loop.

#### Build order

1. `dHRUM` getters/setters (additive, Rcpp-free) → golden-harness check (lumped bit-for-bit vs R path).
2. `Criterion` + `ScalarObjective` (+ `LogTransform` flag) → validate against `hydroGOF::NSE`, `Calibrations/Amalie/dHRUM_SCEoptimBP.r` outputs.
3. `BroadcastMapper` (lumped works end-to-end, zero OMP risk) → `tests/test_calibration.r`.
4. `IdentityMapper` + `DeOptimizer` → matches `Calibrations/Amalie/dHRUM_DeoptimBP.r` terminal NSE, `nthreads=1`.
5. Fix `gatherTsFromHrus` race (Stage 1.2) → re-run distributed at `nthreads>1`, enable trial-parallel optimizer clones.

Two open design choices (must settle before coding):
- Observations as raw `const hdata*` (simple, caller-owned) vs owning `ObservedSeries` (adds warmup mask / NA handling / calendar metadata, duplicates some state).
- Basin aggregate as `hruIndex == dimHM` (compact, matches `get_HbDta`) vs explicit `bool isBasin` (clearer; prefers `isBasin` since it avoids a "magic index" special case in `getHruTsDta`).
