# dHRUM Refactoring Plan

Scope: `single_HMunit`, `dHRUM`, `data_HB_1d`, `params` (plus small changes to `numberSel.h` / `parStructSels.h`).

Guiding constraints:

- **No change to the R-facing API.** All Rcpp glue in `src/dHRU_*_R_Rcpp.cpp` keeps compiling unchanged. Calibration scripts under `Calibrations/` keep working.
- **Behavior-preserving.** Except for the two bug fixes in Stage 1 (which fix wrong behavior), model outputs must be bit-for-bit identical before and after each stage. Build a golden-output harness before touching anything (Stage 0).
- **One stage per commit.** Rebuild and run the harness + `tests/` after every stage.

---

## Stage 0 — Safety net and hygiene (no functional change)

1. Create a branch, e.g. `git checkout -b refactor/data-model`.
2. Remove committed build artifacts and stop tracking them:
   - `src/*.o`, `src/dHRUM.so`, `.Rhistory`; verify `.gitignore` covers them.
3. Golden-output harness (new file, e.g. `tests/golden_run.R`, not shipped in the package):
   - Build the package from `main`, run a representative configuration (a few HRUs, several storage/interception/PET types, ~1 year of daily data, ponds enabled for one HRU), and save **all** `all_ts` series from `get_HbDta()` plus per-HRU series from `getSingleHruTsDta()` to an `.rds` file.
   - After each stage, re-run and compare with `identical()` (valarrays come back as numeric vectors). The harness must be deterministic: single thread for the baseline (`set_num_treads(1)`) so the Stage-1 OMP fix is checked separately (see Stage 1.2).
4. Record a baseline test run: `R CMD INSTALL .` plus the scripts in `tests/` (`test_IntercsptionModels.r`, etc.). Note which ones currently pass so refactors are not blamed for pre-existing failures.

## Stage 1 — Bug fixes

### 1.1 Wrong member copied in `data_HB_1d` copy semantics

`src/data_HB_1d.cpp`, copy constructor (~line 249) and `operator=` (~line 314) both contain:

```cpp
init_CanS = other.init_SteS;   // BUG: should be other.init_CanS
```

Fix both to `init_CanS = other.init_CanS;`. (In Stage 2 these functions are replaced by `= default`, which removes this whole class of bug; still fix it now so the fix is visible in history and covered by its own commit/test.)

**Status: DONE (commit 0015663, 2026-09-06).** The copy-constructor hunk was already fixed in the working tree by the user; the `operator=` hunk was completed alongside.

Add a regression check to the harness: construct an HRU, set distinct init states for CANS vs STES, copy it via `dHRUM` copy/`initHrusVec`, and verify `get_initState(init_Stype::CANS)` survived the copy.

### 1.2 OpenMP data race in `dHRUM::gatherTsFromHrus`

`src/dHRUM.cpp` ~line 458:

```cpp
#pragma omp parallel for num_threads(num_threads)
for(unsigned itHru=0; itHru<dimHM; itHru++) {
  helpValAr += dHruVec[itHru].getSingleHruTsDta(itTS) * Areas[itHru] / basinArea;
}
```

`helpValAr` is a shared accumulator updated from multiple threads — a data race. Basin-aggregated series can be nondeterministically wrong (and `valarray` `+=` on a shared object is not atomic).

Fix: accumulate deterministically. The accumulation is associative in exact arithmetic; to also keep the current floating-point summation order, parallelize over **time steps** instead of HRUs, or sum each HRU into private storage and combine in order:

```cpp
for(unsigned itRts=0; itRts<numTSvars; itRts++) {
  const ts_type itTS = all_ts[itRts];
  hdata helpValAr(0.0, numTs);
  // per-HRU series are hdata temporaries; sum in HRU order
  for(unsigned itHru=0; itHru<dimHM; itHru++) {
    helpValAr += dHruVec[itHru].getSingleHruTsDta(itTS) * Areas[itHru] / basinArea;
  }
  basinDta.s_data(helpValAr, itTS, false);
}
```

If speed matters (many HRUs × long series), parallelize over the time index and loop HRUs inside — each thread owns disjoint `helpValAr[i]` elements and HRU order is preserved. Verify with the harness run at 1 vs. 4 threads: outputs must be identical.

Note: `getSingleHruTsDta()` returns a fresh `hdata` copy per call — 35 series × nHRU × allocation per `gatherTsFromHrus` call. Consider (Stage 3) a `const hdata&` accessor to avoid the copies.

### 1.3 `allParNames` misalignment in `parStructSels.h` — audited: cosmetic only, but exposes a real adjacent bug

`inst/include/parStructSels.h`: `all_pars` has 43 entries; `allParNames` has 44, with `"CAN_ST"` duplicated at positions 10 and 11, shifting every later name off by one relative to the enum.

**Audit result (blast radius: none for `allParNames` itself).** `allParNames` has exactly three consumers, all in `src/dHRU_params_R_Rcpp.cpp` (lines 52, 371, 696), and all three are `std::find` membership checks validating user-supplied name strings. Membership semantics are unaffected by the duplicated `"CAN_ST"` — the 43 valid enum names are all present, so no valid name is rejected and no invalid name accepted. Nothing maps `all_pars[i] -> allParNames[i]` by index. The enum→string output path (`get_param_names` -> `params::par_HRUtype_to_string`, `src/params.cpp:1428`) is a hand-written switch covering all 43 values correctly, and string→enum uses local `std::map<std::string, par_HRUtype>` literals, not `allParNames`.

**Adjacent real bug the audit exposed.** Both `s_mapStringTopar_HRUtype` map literals (`src/dHRU_params_R_Rcpp.cpp:58` and `:377`) contain 42 of 43 entries — `{"SMAXpdm", par_HRUtype::SMAXpdm}` is **missing**, although `allParNames` accepts `"SMAXpdm"` and the downstream switch even has a `case par_HRUtype::SMAXpdm:` (lines 126, 445) intended to receive it. The lookup `s_mapStringTopar_HRUtype[parNameStr[id]]` uses `map::operator[]` on a non-const map, so `"SMAXpdm"` passes validation, inserts a default `par_HRUtype()` (= enum value 0 = `B_SOIL`), and the user's SMAXpdm value is **silently assigned to B_SOIL**. Latent today: no script in `R/`, `tests/`, or `Calibrations/` passes `"SMAXpdm"` (verified by repo-wide grep), and `L_PDM` doesn't request it, which is why it has gone unnoticed — but `setParamsToOnedHru`/`setParamsToAlldHrus` are user-facing exports (`setParsToDistdHRUM` delegates to `setParamsToOnedHru`, so it is affected through the same map).

Fix: add `{"SMAXpdm", par_HRUtype::SMAXpdm}` to both maps, and switch the lookups from `operator[]` to `.at()` (throws on unknown key) so validation/map desync can never silently remap a parameter again. The root-cause fix is Stage 2.3's single name↔enum table shared by validation and mapping; the local per-function map literals (duplicated twice in one file) should be hoisted into one shared `const` map regardless.

## Stage 2 — Enum-indexed storage (data model refactor)

Core idea: the enums are sequential from zero, so they can index `std::array`. Every giant switch collapses to one accessor, and copy/assignment become `= default`.

### 2.1 `numberSel.h`: derive the count instead of hand-maintaining it

Keep `ts_type` values and order unchanged (they may leak into serialized data or R code). Add a sentinel or a constant next to the enum:

```cpp
enum class ts_type {PREC, /* ...unchanged... */ REFR, NUM_TS_TYPES};
constexpr unsigned numTSvars = static_cast<unsigned>(ts_type::NUM_TS_TYPES);
```

and define `all_ts` with `NUM_TS_TYPES` entries. Same pattern for `init_Stype` (12 entries → add `NUM_INIT_STATES`) and `cal_Type`. From now on adding a variable is a two-line change (enum + nothing else) instead of editing 8+ switch statements.

### 2.2 `data_HB_1d`: replace 35 named members with an indexed array

Header (`inst/include/data_HB_1d.h`):

```cpp
private:
  unsigned numTS;
  std::array<hdata, numTSvars> ts;              // was: Prec, Snow, Rain, ...
  std::array<numberSel, static_cast<size_t>(init_Stype::NUM_INIT_STATES)> initState;
  // calendar members (year, month, day, Jday, init_*), Latitude, PETtype,
  // numfastRes, StateFastRes, OutFastRes stay as-is for now
```

Add one private helper:

```cpp
hdata&       at(ts_type t)       { return ts[static_cast<size_t>(t)]; }
const hdata& at(ts_type t) const { return ts[static_cast<size_t>(t)]; }
```

Then each public method shrinks from a ~150-line switch to:

```cpp
void data_HB_1d::s_data(const hdata& dta, const ts_type& t, bool updateNumTS) {
  if (updateNumTS) numTS = dta.size();
  at(t) = dta;
}
numberSel data_HB_1d::g_dta(const unsigned& tst, const ts_type& t) { return at(t)[tst]; }
void data_HB_1d::s_varVal(const numberSel& d, const unsigned& tst, const ts_type& t) { at(t)[tst] = d; }
hdata data_HB_1d::get_HbTsData(const ts_type& t) { return at(t); }
void data_HB_1d::setOneTstoZero(const ts_type& t) { at(t) = 0.0; }   // valarray scalar-assign fills
```

- `s_initStates` / `g_initState` use `initState[static_cast<size_t>(s)]` the same way.
- Delete the hand-written copy ctor / `operator=` and declare both `= default` in the header (all remaining members are value types). This permanently removes the Stage-1.1 bug pattern.
- The default constructor shrinks from a 60-line init list to a loop over `ts` (resize each to `{1}` as today) or, better, empty vectors resized lazily by `s_data` — but **keep the current size-1 initial state** unless callers are audited, since code may read before writing.
- `setAllToZeros` / calendar switches: iterate `all_ts` / index `std::array<caldata,4>`.
- **SRP: extract calendar and PET algorithms into free functions** (decided 2026-09-06): the arrays (`year/month/day/Jday`, `Temp`→`PEt`) stay in `data_HB_1d`; the algorithms move out. New `calendar_util.h/cpp` with a `calendar::` namespace of free functions (`isLeap`, `daysInMonth`, `generateDaily(initY,initM,initD,numTS, y,m,d,jday&)`, `fillJulianDay`), and new `pet_calc.h/cpp` with a `pet::` namespace of pure functions (`oudin/hamon/thornthwaite/blaneycriddle/jensenhaise/mcguinnessbordne`, each `(latitude, temp, year, month, jday) -> hdata`) plus `pet::compute(pet_Type, ...)` dispatching through a `std::array<PetFn,6>` indexed by `static_cast<size_t>(pet_Type)` — the `calc_Pet` switch disappears. `data_HB_1d::s_calender()` / `loadCalData()` / `calc_Pet()` become thin forwarders calling these functions, so `single_HMunit` and the Rcpp glue compile unchanged. Call-order invariants to preserve: `s_initDate` before `s_calender`; `s_Pet_Pars` before `calc_Pet`. `s_calender` currently hard-assumes daily resolution (day-increment logic) — document that; only promote calendar to its own class if sub-daily timesteps are ever planned.
  - Fix while extracting: `BlaneycriddlePET` loops from `tst=1` and never writes `PEt[0]` (stale value). **Confirmed empirically 2026-09-06**: 1-HRU, 10-day June run at lat 50 gives `PET[1] = 0` (the zero fill) vs. formula-expected 5.336; days 2–10 match the formula to full precision. The new `pet::blaneycriddle` must write all steps. **Status: FIXED as standalone commit 0c76cc6 (loop starts at `tst=0`), locked in by tests/pet_regression.r** (fresh-process run: all six methods PASS, BC day 1 = formula value).
  - Note while extracting: the member versions redundantly overwrite `PETtype` at the end of `OudinPET`/`HamonPET`/`ThornthwaitePET` but not in the other three — the pure-function version makes this inconsistency moot (no member to mutate).
  - **PET verification completed for all six methods (2026-09-06)**, 1 HRU, 731 days (2000–2001, leap year), lat 50, seasonal temps: OUDIN, HAMON, JENSENHAISE, MCGUINNESSBORDNE match R transcriptions of the as-coded formulas to machine epsilon over the full series, day 1 included (J-H day-1 = 0 is correct clamping for negative winter temp, not a hole). THORNTHWAITE matches a faithful transcription with 0/731 mismatches, sensible seasonal shape; its loop is 0-based and writes `PEt[0]` explicitly, so only BLAINEYCRIDDLE has the day-1 hole. One latent algorithmic quirk in Thornthwaite (not a boundary bug): the annual heat index for the *final* year in the series sums only Jan–Nov (the loop assigns `I_annual = helpsumI - i_heat[it]` before year-change detection, and the final year never gets its closing update). At lat 50 December mean temps are ≤0 so `tam` clamps to 0 and the impact is nil there, but in warm climates the last year's December is silently dropped from the heat index. **Decision (2026-09-06): KEEP as-is** — the refactor must not change outputs, the current behavior is locked in by `tests/pet_regression.r` (whose transcription deliberately reproduces it) and a comment will mark it in `pet::thornthwaite` during extraction. Revisit only as a deliberate, versioned model change with hydrological evaluation in a warm-climate test basin — not as part of this refactor.
- Optional in the same pass: `dta_dam_1d` has the identical switch pattern (7 switches) — apply the same treatment for consistency.

Risk notes: valarray `operator[]` has no bounds checking, same as today — no regression. `g_dta` currently has no `default` case (returns uninitialized on bad enum); enum-indexing makes out-of-domain enums impossible at compile time for well-formed callers.

### 2.3 `parStructSels.h` + `params`: one table for values and bounds

Replace the three parallel `pdata` valarrays with one record per parameter:

```cpp
// parStructSels.h
enum class par_HRUtype { /* ...unchanged... */, NUM_PARS };
constexpr unsigned numParTypes = static_cast<unsigned>(par_HRUtype::NUM_PARS);

struct ParEntry {
  numberSel val = 0.0;
  numberSel up  = 0.0;
  numberSel low = 0.0;
};

// single source of truth for names — fixes Stage 1.3 by construction
constexpr std::array<const char*, numParTypes> parNames = { "B_SOIL", "C_MAX", /* in enum order */ };
```

`params` then holds `std::array<ParEntry, numParTypes> p;` and:

```cpp
numberSel params::g_par(const par_HRUtype& t)     { return p[idx(t)].val; }
void      params::s_params(const numberSel& d, par_HRUtype t) { p[idx(t)].val = d; }
```

— the four ~180-line switches in `src/params.cpp` disappear. `s_default()` becomes a static table of defaults next to `parNames` (one initializer list instead of 43 switch cases). Keep `PDM_boundary_update()` semantics as-is; document that it mutates bounds as a side effect.

`current_param(...)` keeps its public behavior (building the name/value vectors) but iterates `Current_parameter_list` with `p[idx(...)]` lookups instead of nested switches. The per-structure parameter sets (`L_PDM`, `L_LIN_RES`, ...) are already data-driven lists — keep them untouched; they are the good pattern the rest of the code should imitate.

### 2.4 Checkpoint

Rebuild, run harness (must be bit-identical except the deliberate Stage-1 fixes), run `tests/`. Expect `src/data_HB_1d.cpp` to shrink from ~1,600 to ~600 lines and `src/params.cpp` from ~1,600 to ~700.

## Stage 3 — `dHRUM`: single source of truth

1. **Drop mirrored state.** Remove `dHruVecId`, `gs_STORtypes`, `sw_STORtypes`, `interception_STORtypes`, `surf_STORtypes`, `fast_RESPONSESTypes`, `pondTypes`, `numParsAllHRus`, and `NumFastRes` where they merely mirror `dHruVec[it]` internals.
   - `getHRUIds()` / `getSingleHruId()` → read `dHruVec[it].getIdHru()`.
   - `get_STORtypes()` → collect `dHruVec[it].get_GStype()` on demand.
   - `get_singleHRUnumPars(Id)` → `dHruVec[Id].get_numPars()` directly.
   - Before deleting each vector, verify (grep) that no Rcpp glue reads it as the authoritative copy.
2. **Validate invariants.** `set_numTS()` currently trusts `dHruVec[0]`; check all HRUs report the same `get_numdta()` and throw/`Rcpp::stop` otherwise. Same check for areas length vs. `dimHM` in `setAreasToHrus` (currently a length mismatch silently misbehaves in `gatherTsFromHrus`).
3. **Rename `num_threads`** → `numThreads` (member and `set/get_num_treads` may keep their public names). `num_threads(num_threads)` relying on macro-clause shadowing is fragile.
4. `setParamsToAlldHrus`, `setParamsToOneHru` and the loaders should take `const std::vector<std::pair<numberSel,par_HRUtype>>&` (currently by value — copies the whole vector per call).
5. `basinDta = dHruVec[0].getAllData();` in `initdHRUbasinDTA()` — after Stage 2 this is a defaulted copy and cheap enough; keep, but zero all ts series immediately (as `gatherTsFromHrus` does) so no stale HRU-0 values leak into basin output before gathering.
6. Dead code: delete the large commented blocks (already visible in `calcHbToAllHrus`, `setAreasToHrus`).

## Stage 4 — `single_HMunit`: decompose the god object (incremental, private-only)

Ordered by risk, lowest first:

1. **Privatize public data members** (`Current_par_val`, `Current_uppar_val`, `Current_lowpar_val` — header lines 193–195). Getters already exist (`get_Current_par_values()` etc.); move the vectors private and check the Rcpp glue only uses the getters. If some glue writes them directly, add a setter rather than leaving them public.
2. **Group scratch state.** The 12 `prev_*` members + `et_demand` + `tstRM` become:

```cpp
struct RunState {
  numberSel soil=0, gros=0, canS=0, steS=0, intS=0, snoS=0, surS=0,
            groS1=0, groS2=0, intSnow=0, stemSnow=0, etDemand=0;
  numberDta t=0;
};
RunState st;
```

This makes `set_ZeroStates`/`set_ZeroinitStates` a struct reset and makes the copy constructor trivial: after Stage 2 all members are value types, so delete the hand-written copy ctor / `operator=` in favor of `= default` (audit first that nothing relies on the current copy semantics — e.g. `tstRM` deliberately reset to 0 in the current copy ctor; keep that behavior explicitly if callers depend on it).
3. **Typo'd public API**: add correctly-spelled methods (`update_actualET`, `get_interceptionStorType`, `set_calendar`, `pondMax`) and keep the old names as inline forwarding wrappers marked deprecated in comments. Remove the old names only in a major release — R scripts may call the Rcpp-exposed equivalents.
4. **Extract I/O.** `print_OutputToFile` / `read_InputFromFile` / `print_Pars` move to free functions in a new `single_HMunit_io.{h,cpp}` taking a `const single_HMunit&`. This is the seam that later allows file-formats to change without touching the model class.
5. **Model dispatch stays.** The `switch(_gs_STORtype)`-style dispatchers in `single_HMunit.cpp` (soil, gw, interception, surface, fast response) are fine as variant selection; do **not** convert to virtual polymorphism — HRUs are stored by value in `std::vector<single_HMunit>` and switched per-basin, so enum dispatch is the right design here. Just give every switch a `default:` that throws (silent fall-through today hides configuration errors).
6. Delete the ~commented-out code swaths throughout `single_HMunit.cpp`.

## Stage 5 — SRP: extract `Calendar` and PET from `data_HB_1d`

Motivation: `data_HB_1d` should be a typed container of hydrological series. Calendar is *derived data* (function of start date + length); PET is a *derived computation* (function of calendar + temperature + latitude). Neither needs container internals: all six `*PET()` methods only read `Latitude`, `Temp`, `year`, `month`, `Jday`, `numTS` and write `PEt`. Month-length tables are currently duplicated in three places (`s_calender`, `s_Julian_day`, `get_daysInMonth`), and `dta_dam_1d` holds another calendar implementation.

Noticed defect: `OudinPET()` and `HamonPET()` mutate `PETtype` as a side effect (`src/data_HB_1d.cpp:880, 904`) while the other four PET variants don't. After extraction `s_Pet_Pars` is the single writer of the flag.

1. **New `calendar.h/.cpp`** — pure functions `isLeapYear`, `daysInMonth(month, year)`, `julianDay(y, m, d)`, plus a `Calendar` class owning the year/month/day/Jday arrays and the init date:
   - `generate(startYear, startMonth, startDay, nSteps)` (= today's `s_initDate` + `s_calender` + `s_Julian_day`), `load(yy, mm, dd)` (externally supplied dates, = `loadCalData`), `series(cal_Type)`, `valueAt(cal_Type, ts)`, `size()`.
   - Store the four series in one `std::array<caldata, 4>` indexed by `cal_Type` (add a `NUM_CAL_TYPES` sentinel per Stage 2.1) — collapses the `g_calDta` / `getCalData` / `setOneCalDateToZero` switches.
2. **New `pet.h/.cpp`** — pure functions in `namespace pet`: `oudin`, `hamon`, `thornthwaite`, `blaneycriddle`, `jensenhaise`, `mcguinness`, each `(const Calendar&, const hdata& temp, numberSel latitude) -> hdata`, plus one `compute(pet_Type, ...)` dispatcher. Verified from the implementations that `(Calendar, Temp, Latitude)` covers all six, including Thornthwaite's monthly aggregation. Free functions chosen over strategy classes for the same reason as Stage 4.5 (stateless, enum-dispatched) — and statelessness makes them trivially safe in the OpenMP HRU loop and unit-testable without a `data_HB_1d`.
3. **Rewrite the seams in `data_HB_1d`.** Replace the seven calendar members with `Calendar cal;`; keep the existing public methods (`s_calender`, `s_initDate`, `g_calDta`, `getCalData`, `loadCalData`, `setOneCalDateToZero`, `leap_Check_Year`, `get_daysInMonth`) as one-line forwarders — zero caller changes. Replace the six `*PET()` methods with `at(ts_type::PET) = pet::compute(PETtype, cal, at(ts_type::TEMP), Latitude);` inside `calc_Pet`. `Latitude`/`PETtype` may stay as members at first (config, not data); moving them into `single_HMunit`'s configuration is optional follow-up.
4. **Reuse `Calendar` in `dta_dam_1d`**, deleting its duplicated calendar logic.
5. Optional: move the PET call up to `single_HMunit::calc_Pet` (`hyd_dta.s_data(pet::compute(...), ts_type::PET, false)`) so `data_HB_1d` knows nothing about PET at all. Do this only after the forwarder version is verified.

Verification: calendar and PET are pure, so the golden harness catches any migration error. Expect another ~400 lines removed from `data_HB_1d.cpp` and the first real unit-testable units in the codebase.

## What NOT to change (and why)

- **`std::valarray` → `std::vector`.** Tempting, but `valarray`'s elementwise arithmetic (`helpValAr += ts * area / basinArea`) is used genuinely and widely; a swap would touch every hydrological equation for zero behavioral gain. Revisit only if profiling shows the temporaries matter (then the fix is expression-local, e.g. the Stage-3 accumulator loop, not a container swap).
- **Enum dispatch → inheritance** for model variants: see Stage 4.5.
- **R-facing names:** out of scope by constraint.

## Verification checklist (run after every stage)

0. **Fresh R session mandatory after `R CMD INSTALL`.** Learned the hard way (2026-09-06): re-installing while dHRUM is loaded in the RStudio session corrupts the lazy-load DB, and even after `unloadNamespace()` + `library()` the process keeps executing the **stale DLL image** (dlopen caching + Rcpp XPtr finalizer references), so probes silently test the previous build. Golden-harness and regression runs must go through `Rscript` or Session → Restart R. Snapshot-based comparison (binary-identical golden outputs) makes this workable without trusting in-session reloads.
1. `Rcpp::compileAttributes()` if any Rcpp-visible signature changed (Stages 1–4 should not need it; Stage 4.3/4.4 might).
2. `R CMD INSTALL .` clean, no new compiler warnings (`-Wall` noise in current build is a pre-existing baseline — record it in Stage 0).
3. Golden harness: `identical()` against baseline (Stages 2–4) / documented-diff-only (Stage 1).
4. Determinism check: harness at 1 vs. 4 threads gives identical output.
5. `tests/` scripts run as they did at baseline.
6. One real calibration from `Calibrations/POH/` runs end-to-end and reaches the same objective value.


