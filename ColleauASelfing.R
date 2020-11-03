source("recursiveASelfing.R")
source("AInvSelfing.R")

#pedigri: id, sire, dam, self
#w vector por el que su multiplica la matriz
#f vector de inbreedings

colleauSelfing <- function(ped, w, f) {
  n <- dim(ped)[1]
  q <- matrix(0, n)
  v <- matrix(0, n)
  for (i in n:1) { #para cada animal del pedigri, del ultimo al primero
    is <- ped[i,2]
    id <- ped[i,3]
    q[i] <- q[i] + w[i]
    if (is > 0) { #si sire conocido 
      q[is] <- q[is] + q[i] / 2 } #recibe 0.5 de lo del hijo
    if (id > 0) { #si dam conocida
      q[id] <- q[id] + q[i] / 2 } #recibe 0.5 de lo del hijo
  }
  for (i in 1:n) { #para cada animal del pedigri, del primero al ultimo (al revés que antes)
    is <- ped[i,2]
    id <- ped[i,3]
    iself <- ped[i,4]
    fs <- 0
    fd <- 0
    if (is > 0) { #si sire conocido
      fs <- f[is] } #asignale su f
    if (id > 0) { #si dam conocida
      fd <- f[id] } #asignarle su f
    di <- 1/findAlphaSelfing(ped, is, id, iself, f) #elementos de la d (no. de padres desconocidos + f de los padres)
    if (is > 0 & id > 0) { # si sire y dam conocidos
      v[i] <- (v[is] + v[id]) / 2} #el del individuo es el promedio de los padres
    if (is > 0 & id == 0) { # si sire conocido y dam desconocida
      v[i] <- v[is] / 2} #el del individuo es la mitad del del sire
    if (id > 0 & is == 0) { # si dam conocida y sire desconocido
      v[i] <- v[id] / 2 } #el del individuo es la mitad del de la dam
    v[i] <- v[i] + di * q[i] #el v[i] + 0.5*q que se había calculado antes
  }
  return(v)
}
