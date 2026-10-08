koef =0.125
ee=3
S=seq(0,1000, by=0.1)

r=1/(koef*S^ee+1)

opr=koef*S^ee/(koef*S^ee+1)
plot(r, col="red",type ="l", ylim=c(0,1))
lines(opr, col="blue")


koef =0.125
ee=0
S=seq(0,1000, by=0.1)

r=1/(koef*S^ee+1)

opr=koef*S^ee/(koef*S^ee+1)
plot(r, col="red",type ="l", ylim=c(0,1))
lines(opr, col="blue")


koef =0.125
ee=1
S=seq(0,1000, by=0.1)

r=1/(koef*S^ee+1)

opr=koef*S^ee/(koef*S^ee+1)
plot(r, col="red",type ="l", ylim=c(0,1))
lines(opr, col="blue")

koef =0.125
ee=3
S=seq(0,1000, by=0.1)

r=1/(koef*S^ee+1)

opr=koef*S^ee/(koef*S^ee+1)
plot(r, col="red",type ="l", ylim=c(0,1))
lines(opr, col="blue")

koef =0.19
ee=0.75
S=seq(0,1000, by=0.1)

r=1/(koef*S^ee+1)

opr=koef*S^ee/(koef*S^ee+1)
plot(S,r,col="red",type ="l", ylim=c(0,1))
lines(S, opr, col="blue")

