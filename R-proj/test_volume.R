library(volesti)

P = gen_cube(10,'H')

tim = system.time({ vol = test_volume(P) })

print(vol)
print(tim)
