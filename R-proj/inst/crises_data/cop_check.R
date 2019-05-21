library(volesti)
library(tawny)

starting_date = "2007-09-04"
stopping_date = "2009-01-05"

MatReturns <- read.table("https://stanford.edu/class/ee103/data/returns.txt", sep=",")

MatReturns = MatReturns[-c(1,2),]
dates = as.character(MatReturns$V1)
MatReturns = as.matrix(MatReturns[,-c(1,54)])
MatReturns = matrix(as.numeric(MatReturns [,]),nrow = dim(MatReturns )[1],
                    ncol = dim(MatReturns )[2], byrow = FALSE)

row1 = which(dates %in% starting_date)
row2 = which(dates %in% stopping_date)

cp = copulas(MatReturns[row1:row2,])

indicators=c()
numSlices=100
for (i in 1:(dim(cp[1]/52))) {
  
cop = cp[((i-1)*52+1):i*52,]
blue_mass = 0
red_mass = 0

for (row in 1:100) {
  for (col in 1:100) {
    if (row-col<=0.2*numSlices && row-col>=-0.2*numSlices) {
      if (row+col<0.8*numSlices || row+col>1.2*numSlices) {
        red_mass = red_mass + cop[row,col]
      }
    } else {
      if (row+col>=0.8*numSlices && row+col<=1.2*numSlices) {
        blue_mass = blue_mass + cop[row,col]
      }
    }
  }
}
indicators = c(indicators, blue_mass / red_mass)
}
