library(volesti)

d=40
P = GenCube(d,'H')

#print(exact_vol(P))

vol = volume(P, Algo = "BAN", WalkType = "BilW", diameter = sqrt(d))

print(vol)

vol = volume(P, Algo = "BAN")

print(vol)
