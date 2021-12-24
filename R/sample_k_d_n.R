sample_k_d_n <- function(m, m_i, childrenCount, childrenCountSpp, ...){

  while (m > 0 & sum(childrenCount > m_i) <= m) {
    child      <- names(childrenCount)[childrenCount > m_i]
    m_i[child] <- m_i[child] + 1
    m          <- m - length(child)
  }
  # If there are still sequences to sample, get one using the taxon
  # diversity as probability values
  if (m > 0) {
    child      <- sample(names(childrenCount)[childrenCount > m_i],
                         size = m,
                         prob = childrenCountSpp[childrenCount > m_i] / sum(childrenCountSpp[childrenCount > m_i]))
    m_i[child] <- m_i[child] + 1
  }

  return(m_i)
}
