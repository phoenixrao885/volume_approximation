library(volesti)

#P = gen_cube(100,'H')
path = system.file('extdata', package = 'volesti')
P = file_to_polytope(paste0(path,'/birk10.ine'))

tim = system.time({ vol = test_volume(P, error=0.05, walk_length = 1, parameters = list("Window" = 1000)) })

print(vol)
print(tim)
