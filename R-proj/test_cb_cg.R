library(volesti)

dims = seq(from = 210, to = 500, by =10)

time_cb=c()
time_cg=c()

for (d in dims) {
  
  P = gen_cube(d, 'H')
  
  tims1=0
  tims2=0
  
  for (j in 1:5) {
    tim1 = system.time({ vol = test_volume(P) })
    tim2 = system.time({ vol = volume(P, algo='CG') })
    
    tims1 = tims1 + as.numeric(tim1)[3]
    tims2 = tims2 + as.numeric(tim2)[3]
  
  }
  
  time_cb = c(time_cb, tims1/5)
  time_cg = c(time_cg, tims2/5)
  
  save(time_cb, file = "times_cb_210_500.RData")
  save(time_cg, file = "times_cg_210_500.RData")
  
  print(d)
  
}
