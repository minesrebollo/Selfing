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

# AInv <- function(ped) {
#   n <- nrow(ped) #n = número de animales en el pedigree
#   ainv <- matrix(0, n, n) #Crea matriz nxn con 0s
#   w <- matrix(c(1, -0.5, -0.5), 3, 1) #el vector para las contribuciones
#   for (ii in 1:n) { #para cada animal
#     a <- ped[ii,1] #a = el animal
#     s <- ped[ii,2] #s = el padre
#     d <- ped[ii,3] #d = la madre
#     x <- findAlpha(ped, s, d) #función que devuelve la contribución según npadres
#     p <- matrix(c(a, s, d), 1, 3) #matriz 1x3 con el numero de a, s y d
#     for (i in 1:3) {
#       for (j in 1:3) {
#         if (p[i] > 0 & p[j] > 0) {  
#           ainv[p[i], p[j]] <- ainv[p[i], p[j]] +  w[i] * w[j] * x
#         }
#       }
#     }
#   }
#   return(ainv) #la función devuelve a inv una vez que recorrió todos los animales
# }

#Para el seminario
# AInvInb <- function(ped,Fs) {
#   n <- nrow(ped)
#   ainv <- matrix(0, n, n)
#   w <- matrix(c(1, -0.5, -0.5), 3, 1)
#   for (ii in 1:n) {
#     a <- ped[ii,1]
#     s <- ped[ii,2]
#     d <- ped[ii,3]
#     fs <- Fs[s]
#     fd <- Fs[d]
#     x <- findAlpha(ped, s, d, fs, fd)
#     p <- matrix(c(a, s, d), 1, 3)
#     for (i in 1:3) {
#       for (j in 1:3) {
#         if (p[i] > 0 & p[j] > 0) {
#           ainv[p[i], p[j]] <- ainv[p[i], p[j]] +  w[i] * w[j] * x
#         }
#       }
#     }
#   }
#   return(ainv)
# }

# AInvInbUPG <- function(ped,Fs) {
#   n <- nrow(ped)
#   ainv_upg <- matrix(0, n+ngr, n+ngr)
#   w <- matrix(c(1, -0.5, -0.5), 3, 1)
#   for (ii in 1:n) {
#     a <- ped[ii,1]
#     s <- ped[ii,2]
#     d <- ped[ii,3]
#     fs <- Fs[s]
#     fd <- Fs[d]
#     x <- find_upgcode(s, d, fs, fd)
#     if (s < 0) s <- (s * -1) + n
#     if (d < 0) d <- (d * -1) + n
#     p <- matrix(c(a, s, d), 1, 3)
#     for (i in 1:3) {
#       for (j in 1:3) {
#         if (p[i] > 0 & p[j] > 0) { #esta condición siempre se cumple
#         ainv_upg[p[i],p[j]] <- ainv_upg[p[i],p[j]] + w[i] * w[j] * x
#       }
#     }
#   }
#   return(ainv_upg)
# }









