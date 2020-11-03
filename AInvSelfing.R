# Código Ainv con selfing, adaptado del codigo de A con Inbreeding -----------------
# En base a lo de Oaxley a pruebas y errores con distintos ejemplos

# Esta es la clave:
# Alphas considerando generaciones de autofecundación y cantidad de padres conocidos:
findAlphaSelfing <- function(ped, s, d, self, f) {
  if (s > 0 & d > 0) { #si padre y madre conocidos
    x <- 1/(1 + (1 - 0.5^self) * (1 - (0.5 * cffa(ped, s, d))) - 1/4 * ((1 + f[s]) + (1 + f[d]))) #antes acá estaba find_alphaInb
  }
  else if (s > 0) { #solo padre conocido
    x <- 1/(1 + 1 - 0.5^self - 1/4 * (1 + f[s]))
  }
  else if (d > 0) { #solo madre conocida
    x <- 1/(1 + 1 - 0.5^self - 1/4 * (1 + f[d]))
  }
  else { # (s > 0) { #ningun padre conocido
    x <- 1/(1 + 1 - 0.5^self)
  }
  return(x)
}

# (self = NULL, f = NULL)
# if f != NULL una cosa
# estaria bueno que se de cuenta si hay selfing o no y segun eso haga finalphaselfing o haga findalpha comun


AInvSelfing <- function(ped) {
  n <- nrow(ped) #numero de filas del pedigri
  ainv <- matrix(0, n, n) #inicializar la matriz con 0s, nped x nped
  w <- matrix(c(1, -0.5, -0.5), 3, 1) #vector para contribuciones
  f <- c()
  for (i in 1:n) {
    f[i] <- inbSelfing(ped, i)
  }
  for (ii in 1:n) {
    a <- ped[ii,1]
    s <- ped[ii,2]
    d <- ped[ii,3]
    self <- ped[ii,4]
    x <- findAlphaSelfing(ped, s, d, self, f)
    p <- matrix(c(a, s, d), 1, 3) #vector que dice el animal, el sire y el dam
    for (i in 1:3) {
      for (j in 1:3) {
        if (p[i] > 0 & p[j] > 0) {  
          ainv[p[i], p[j]] <- ainv[p[i], p[j]] +  w[i] * w[j] * x
        }
      }
    }
  }
  return(ainv)
}









