## Rscript to assess convergence of parameters in a hdf5 file
## output by entropy 

## Sep 2019 -- vshastry
## Downloaded 4/28/22 by Alyssa Phillips

library(rhdf5)

args<-commandArgs(TRUE)

f<-as.character(args[1])

w.q<-h5read(f,"q")

jpeg(file = paste0(f, "convergence.jpeg"))
par(mfrow=c(2,2))
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='')
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='')
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='')
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='')
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='')
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='')
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='')
plot(w.q[,sample(1:dim(w.q)[2],1),sample(1:dim(w.q)[3],1)],type="l",ylab='')
mtext("Trace plots for admixture estimates", outer = TRUE, side=3, line=-2)
dev.off()

#w.p<-h5read(f,"p")

#par(mfrow=c(2,2))
#plot(w.p[,sample(1:dim(w.p)[2],1),sample(1:dim(w.p)[3],1)],type="l",ylab='')
#plot(w.p[,sample(1:dim(w.p)[2],1),sample(1:dim(w.p)[3],1)],type="l",ylab='')
#plot(w.p[,sample(1:dim(w.p)[2],1),sample(1:dim(w.p)[3],1)],type="l",ylab='')
#plot(w.p[,sample(1:dim(w.p)[2],1),sample(1:dim(w.p)[3],1)],type="l",ylab='')
#mtext("Trace plots for allele freq. estimates", outer = TRUE, side=3, line=-2)

