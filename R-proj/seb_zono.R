
#d=15
#P=GenZonotope(d,2*d,dist = "uniform")
#voliter = volume(P, WalkType = "RDHR", diam_zono = diam_hpoly_zono_inter)
#zono_seb_unif =c(zono_seb_unif, voliter[3])


zono_seb_gauss = c()
d=5
P=GenZonotope(d,2*d,dist = "gaussian")
voliter = volume(P, WalkType = "RDHR", diam_zono = diam_hpoly_zono_inter)
zono_seb_gauss =c(zono_seb_gauss, voliter[3])

d=10
P=GenZonotope(d,2*d,dist = "gaussian")
voliter = volume(P, WalkType = "RDHR", diam_zono = diam_hpoly_zono_inter)
zono_seb_gauss =c(zono_seb_gauss, voliter[3])

d=15
P=GenZonotope(d,2*d,dist = "gaussian")
voliter = volume(P, WalkType = "RDHR", diam_zono = diam_hpoly_zono_inter)
zono_seb_gauss =c(zono_seb_gauss, voliter[3])



zono_seb_exp = c()
d=5
P=GenZonotope(d,2*d,dist = "exponential")
voliter = volume(P, WalkType = "RDHR", diam_zono = diam_hpoly_zono_inter)
zono_seb_exp =c(zono_seb_exp, voliter[3])

d=10
P=GenZonotope(d,2*d,dist = "exponential")
voliter = volume(P, WalkType = "RDHR", diam_zono = diam_hpoly_zono_inter)
zono_seb_exp =c(zono_seb_exp, voliter[3])

d=15
P=GenZonotope(d,2*d,dist = "exponential")
voliter = volume(P, WalkType = "RDHR", diam_zono = diam_hpoly_zono_inter)
zono_seb_exp =c(zono_seb_exp, voliter[3])
