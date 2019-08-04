library(volesti)

d=80
W=30
P = GenZonotope(30,60)
#P=GenRandVpoly(20,40)

#print(exact_vol(P))
di = comp_diam(P)
t1=system.time({retvec1 = volume(P, Algo = "BAN", WalkType = "BilW", Parameters = list("Window"=30, "hpoly"=TRUE),diameter = 0.4*di,diam_zono = diam_hpoly_zono_inter)})


