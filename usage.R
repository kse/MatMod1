s1 <- mm_analyzeSums(93107.04, 747.4, 6)
s2 <- mm_analyzeSums(93753.64, 750, 6)

x <- data.frame(Grøn=c(5,13,16), Rød=c(101,213,185), Sort=c(94,174,199))

nord <- c(7.1, 7.2, 7.4, 7.6, 7.6, 7.7, 7.7, 7.9, 8.1, 8.4, 8.5, 8.8)
øst  <- c(6.9, 7.0, 7.1, 7.2, 7.3, 7.3, 7.4, 7.6, 7.8, 8.1, 8.3, 8.5)
syd  <- c(7.8, 7.9, 8.1, 8.3, 8.3, 8.4, 8.4, 8.4, 8.6, 8.9, 9.2, 9.4)
vest <- c(6.4, 6.6, 6.7, 7.1, 7.6, 7.8, 8.2, 8.4, 8.6, 8.7, 8.8, 8.9)

val <- rbind(mm_analyzePoints(nord), mm_analyzePoints(øst), mm_analyzePoints(syd), mm_analyzePoints(vest))
