library(volesti)

P = GenRandVpoly(15,2*15)

timiter = system.time({ voliter = volume(P, Algo = "BAN", WalkType = "RDHR", Parameters = list("Window"=30), rounding = TRUE, diam_zono = diam_hpoly_zono_inter) })

print(as.numeric(timiter[3]))
