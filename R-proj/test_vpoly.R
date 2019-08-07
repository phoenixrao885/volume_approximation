library(volesti)

d=15
P = GenRandVpoly(d,2*d, body = 'cube')
W=100+ceiling(6*sqrt(d))
#timiter = system.time({ voliter = volume(P, Algo = "BAN", WalkType = "RDHR", Parameters = list("Window"=30), rounding = TRUE, diam_zono = diam_hpoly_zono_inter) })

voliter = volume(P, Algo = "BAN", WalkType = "BilW", Parameters = list("Window"=W, "N"=120+0.6*sqrt(d)), iii = 0.9, diam_zono = comp_diam_hpoly_zono)



print(voliter)

voliter = volume(P, Algo = "BAN", WalkType = "BilW", rounding = TRUE, Parameters = list("Window"=W, "N"=120+0.6*sqrt(d)), iii = 0.9, diam_zono = comp_diam_hpoly_zono)
print(voliter)
