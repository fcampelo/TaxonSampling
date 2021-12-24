sample_k_d_a <- function(m, m_i, childrenCount, childrenCountSpp, ...){

  if (m > 0 & sum(childrenCount > m_i) <= m) {
    child       <- names(childrenCount)[childrenCount > m_i]
    m_i[child]  <- m_i[child] + 1
    m           <- m - length(child)
  }
  while (m > 0 & sum(childrenCount > m_i) > 0) {
    child      <- sample(names(childrenCount)[childrenCount > m_i],
                         size = 1,
                         prob = childrenCountSpp[childrenCount > m_i] / sum(childrenCountSpp[childrenCount > m_i]))
    m_i[child] <- m_i[child] + 1
    m          <- m - 1
  }

  return(m_i)
}
