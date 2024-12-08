s1 <- c("A-A-A-A-A-A-A-A-A-A-A-A")
s2 <- c("A-A-A-A-A-A-A-A-A-A-A-B")
s3 <- c("A-A-A-A-B-B-B-B-C-C-C-C")
s4 <- c("A-A-A-A-B-B-B-B-C-C-C-D")
s5 <- c("A-A-B-B-B-B-B-B-C-C-C-D")
s6 <- c("B-B-B-B-C-C-C-C-D-D-D-D")

w <- c(10,1,30,4,15,40)

ex2.weighted <- data.frame(seq=c(s1, s2, s3, s4, s5, s6), weight=w)
ex2.unweighted <- data.frame(seq=c(rep(s1,w[1]), rep(s2,w[2]), rep(s3,w[3]), rep(s4,w[4]), rep(s5,w[5]), rep(s6, w[6])))

rm(s1, s2, s3, s4, s5, s6, w)
