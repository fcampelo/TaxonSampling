sample_a_d_a <- function(m, m_i, childrenCount, childrenCountSpp, ...){
  # First round
  if (m > 0 & sum(childrenCount > m_i) <= m){
    child       <- names(childrenCount)[childrenCount > m_i]
    m_i[child]  <- m_i[child] + 1
    m           <- m - length(child)
  }

  # Subsequent rounds
  m_i <- sample_a_d_y(m, m_i, childrenCount, childrenCountSpp, ...)

  return(m_i)
}
