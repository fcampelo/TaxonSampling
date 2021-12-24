sample_k_b_y <- function(m, m_i, childrenCount, childrenCountSpp, ...){
  m_i <- table(sample(names(childrenCount),
                      size    = m,
                      replace = TRUE,
                      prob    = childrenCountSpp / sum(childrenCountSpp)))

  return(m_i)
}
