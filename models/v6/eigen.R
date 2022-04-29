

library(popbio)

sp

A <- sp$p
A[1, 4] <- sp$fecundity
A <- A[1:4,]

eigen.analysis(A)
# explanation: http://www.uwyo.edu/dbmcd/popecol/feblects/lect11.html#:~:text=Matrix%20models%20tell%20us%20about,on%20the%20population%20growth%20rate.

