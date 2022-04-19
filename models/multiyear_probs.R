


# number of time steps
t <- 5 

# outcomes of saplings after t years
p <- sample(1:10, 3)
names(p) <- c("s", "a", "m") # sapling, adult, mortality

# annual survival probability
s <- (p[1] / sum(p)) ^ (1 / t) 

# remaining annual probability allocated proportionally
am <- (1 - s) * p[2:3] / sum(p[2:3]) 

# annual transition probabilities
A <- c(s, am)


# iterative simulation with computed transition probabilities
n0 <- c(s = sum(p), a = 0, m = 0) # starting counts (sapling, adult, dead)
n <- n0 # container to accumulate individuals over multiple years
ni <- n0 # temporary container 
for(i in 1:t){
  ni <- n[1] * A # outcomes for remaining saplings
  n[2:3] <- n[2:3] + ni[2:3] # accumulate fates for outgoing saplings
  n[1] <- ni[1]
}

# confirm equivalence
n
p
