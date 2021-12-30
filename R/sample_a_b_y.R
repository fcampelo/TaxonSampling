sample_a_b_y <- function(m, m_i, childrenCount, childrenCountSpp, ...){
  x <- table(sample(names(childrenCount), m, replace = TRUE))

  m2 <- 0 * x
  for (nn in names(m_i)){
    m2[nn] <- m_i[nn]
    if(!(nn %in% names(x))) {
      x[nn]  <- 0
    }
  }

  while(any(x - m2 < 0)){
    id1 <- which.max(x - m2)[1]
    id2 <- which.min(x - m2)[1]
    x[id2] <- x[id2] + 1
    x[id1] <- x[id1] - 1
  }

  return(x)
}
