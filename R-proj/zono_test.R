library(volesti)

d=80
W=30
P = GenCross(d,'V')
#P=GenRandVpoly(20,40)

#print(exact_vol(P))

t1=system.time({retvec1 = volume(P, Algo = "BAN", WalkType = "BilW", Parameters = list("Window"=W), rounding = TRUE)})

print(retvec1)
print(as.numeric(t1[3]))

t2=system.time({retvec2 = volume(P, Algo = "BAN", WalkType = "RDHR", rounding = TRUE)})

print(retvec2)
print(as.numeric(t2[3]))

print(paste0("error1 = ",abs(2^d/prod(1:d)-retvec1[1])/(2^d/prod(1:d))))
print(paste0("error2 = ",abs(2^d/prod(1:d)-retvec2[1])/(2^d/prod(1:d))))