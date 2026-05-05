# the inteceprtion component
# the vegetation fraction depends on the lai
k=0.4
lai=seq(0.1,8,by=0.1)
c=1-exp(-k*lai)
c
plot(c)


k=seq(0.1,0.8,by=0.01)
lai=5
c=1-exp(-k*lai)
c
plot(c)
