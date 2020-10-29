# Adaptación a R del algoritmo de Aguilar & Misztal 2008 

# Assumes individuals renumbered from 1 to n, oldest to youngest 
# n number of individuals
# ped is the pedigree matrix (4 x n) id, sire, dam, selfing generations

# Se podrían guardar elementos calculados en una matriz sparse y queda mucho más eficiente
# No tengo ciclos by default pero no es necesario porque estoy dando el pedigri entero, incluyendo a todos los padres

# Function that returns inbreeding coefficient for individual n considering selfing cycles 
inbSelfing <- function(ped, a) {
  s <- ped[a,2] 
  d <- ped[a,3]
  self <- ped[a,4]
  if (s == 0 | d == 0) { #si alguno de los padres es desconocido 
    inbreedc <- (1 - 0.5 ^ self) # solo generaciones de autofecundación
  }
  else {
    inbreedc <- (1 - 0.5 ^ self * (1 - (0.5 * (cffa(ped, s, d))))) # generaciones de autofecundación y parentesco entre los padres
  }
  return(inbreedc)
}

# Recursive function that returns relationship between a1 and a2 
cffa <- function(ped, a1, a2) {
  #a1 <- s; a2 <- d
  if (a1 == 0 | a2 == 0) { # si padre desconocido
    rel <- 0 }
  else if (a1 == a2) { # si es el individuo consigo mismo
    rel <- 1 + inbSelfing(ped, a1) }
  else if (a1 < a2) {
    rel <- 0.5 * (cffa(ped, a1, ped[a2,2]) + cffa(ped, a1, ped[a2,3])) 
  }
  else {
    rel <- 0.5 * (cffa(ped, a2, ped[a1,2]) + cffa(ped, a2, ped[a1,3]))
  }
  return(rel)
}



