sample_k_b_y <- function(m, m_i, childrenCount, childrenCountSpp, ...){
  x <- table(sample(names(childrenCount),
                    size    = m,
                    replace = TRUE,
                    prob    = childrenCountSpp / sum(childrenCountSpp)))

  while(any(x - m_i < 0)){
    id1 <- which.max(dif)[1]
    id2 <- which.min(dif)[1]
    x[id2] <- x[id2] + 1
    x[id1] <- x[id1] - 1
  }

  return(x)
}
