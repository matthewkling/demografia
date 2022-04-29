

# test softmax priors



p <- make_priors()$p



## two-level ##

# propagules:
x <- log(p[,"p"])
x <- x[is.finite(x)]
x <- x - min(x)
mu_a = x[1]
mu_b = x[2]
sigma = 2

# adults:
x <- log(p[,"a"])
x <- x[is.finite(x)]
x <- x - min(x)
mu_a = x[1]
mu_b = x[2]
sigma = .75



n <- 10000
d <- data.frame(la = rnorm(n, mu_a, sigma),
                lb = rnorm(n, mu_b, sigma)) %>%
      rowwise() %>%
      mutate(a = softmax(c(la, lb))[1],
             b = softmax(c(la, lb))[2])

d %>%
      ggplot(aes(b)) +
      geom_histogram(position = "identity", alpha = .5, bins = 100, boundary = 1)  +
      scale_x_log10(breaks = c(.00001, .0001, .001, .01, .1, 1))

d %>%
      select(a:b) %>%
      gather(l, p, a:b) %>%
      ggplot(aes(p, color = l, fill = l)) +
      geom_histogram(binwidth = .005, boundary = 1,
                     position = "identity", alpha = .5) +
      xlim(0, 1)



## three-level ##


# seedlings
x <- log(p[,"s"])
x <- x[is.finite(x)]
x <- x - min(x)
mu_a = x[1]
mu_b = x[2]
mu_c = x[3]
sigma = 1

# juveniles
x <- log(p[,"j"])
x <- x[is.finite(x)]
x <- x - min(x)
mu_a = x[1]
mu_b = x[2]
mu_c = x[3]
sigma = 2



n <- 10000
d <- data.frame(la = rnorm(n, mu_a, sigma),
                lb = rnorm(n, mu_b, sigma),
                lc = rnorm(n, mu_c, sigma)) %>%
      rowwise() %>%
      mutate(a = softmax(c(la, lb, lc))[1],
             b = softmax(c(la, lb, lc))[2],
             c = softmax(c(la, lb, lc))[3])


d %>%
      select(a:c) %>%
      gather(l, p, a:c) %>%
      ggplot(aes(p, color = l, fill = l)) +
      geom_histogram(binwidth = .005, boundary = 0,
                     position = "identity", alpha = .5) +
      xlim(0, 1)


