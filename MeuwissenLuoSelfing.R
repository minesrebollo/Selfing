#Meuwissen & Luo 1992 adapted to account for selfing generations
#Requires recursiveASelfing

source("recursiveASelfing.R")

meuwSelfing <- function(ss, dd, nself) {
  n <- length(ss) #antes decía size(ss, dim = 1)
  point <- rep(0, n) #the next oldest ancestor in the link list, = 0 if i is the last ancestor
  l <- rep(0, n) # element ij of matrix L
  d <- c() 
  f <- c() 
  for (i in 1:n) { 
    is <- ss[i] # sire
    id <- dd[i] # dam
    ss[i] <- max(is, id) #el mayor queda en ss sin importar si es realmente sire o si es dam
    dd[i] <- min(is, id) #el menor queda en dd sin importar si es realmente dam o si es sire
    #acá va lo de inbreeding, es lo mismo que findAlphaSelfing pero sin 1/, lo puedo poner como función aparte... 
    self <- nself[i]
    if (is > 0 & id > 0) { #si padre y madre conocidos
      d[i] <- 1 + (1 - 0.5^self) * (1 - (0.5 * cffa(ped, is, id))) - 1/4 * ((1 + f[is]) + (1 + f[id]))
    } else if (is > 0) { #solo padre conocido
      d[i] <- 1 + 1 - 0.5^self - 1/4 * (1 + f[is])
    } else if (id > 0) { #solo madre conocida
      d[i] <- 1 + 1 - 0.5^self - 1/4 * (1 + f[id])
    } else { # (s > 0) { #ningun padre conocido
      d[i] <- 1 + 1 - 0.5^self
    }
    if (is == 0 | id == 0) { #si padre o madre desconocidos
      f[i] <- 1 - 0.5 ^ self #el inbreeding es solo lo de sus generaciones de selfing, sale del loop y entra con el siguiente (antes decía 0)
      next }
    #esto que viene lo hace si padre y madre son conocidos porque sino entró al for con el i sieguiente pon el next
    fi <- -1 #supuestamente F(i) is initially set equal to -1, which is more accurate than calculating F + 1 and then subtracting 1.
    l[i] <- 1.0 #la contribución a si mismo es 1 (?)
    j <- i #oldest ancestor of i al principio es el mismo
    while (j != 0) { #mientras lista de ancestros no esté vacía
      k <- j #variable temporal
      r <- 0.5 * l[k] #calcular contribución a/de los padres como 0.5 la contribución de ese animal
      ks <- ss[k] #padre1 del ancestro
      kd <- dd[k] #padre2 del ancestro
      if (ks > 0) { #si el padre1 de k es conocido
        while (point[k] > ks) { #stop if next ancestor is older(o será younger?) than sire)
          k <- point[k] }
        l[ks] <- l[ks] + r #agrega la contribución del sire
        if (ks != point[k]) { #si ks no está como ancestro de k
          point[ks] <- point[k] #se borra lo que diga en point[ks] y se sustituye por point[k]
          point[k] <- ks } #ks se agrega a la lista como ancestro en la posición k
        if (kd > 0) { #si madre conocida
          while (point[k] > kd) { #si ks > kd (en qué situación esto sería falso? no se definieron ks > kd? si son iguales?)
            k <- point[k] } #k pasa a ser el sire 
          l[kd] <- l[kd] + r #agrega la contribución de la dam
          if (kd != point[k]) {
            point[kd] <- point[k] #el sire pasa a la posición de la dam
            point[k] <- kd } #la dam pasa a la posición del sire
          }
        }
      fi <- fi + l[j] * l[j] * d[j]
      l[j] <- 0 #borra esa contribución
      k <- j #guarda el j en k (k había pasado a ser el sire de j)
      j <- point[j] #nuevo j es el siguiente animal más viejo en la lista que quedó guardado en la posición del animal anterior
      point[k] <- 0 #elimina ese ancestro de la lista
    }
    f[i] <- fi
    ###
    if (f[i] < 0 & abs(f[i]) < 0.001) { ## ?
      f[i] <- 0 }
  }
  return(f)
}

#ped <- data.frame(id = c(1:16),
#                  sire = c(0,0,0,0,0,1,3,5,6,4,8,1,10,8,8,8),
#                  dam = c(0,0,0,0,0,2,2,0,7,7,0,9,9,13,14,14),
#                  nself = c(0,2,0,0,5,0,0,5,0,0,0,0,3,0,0,3)); ped

#ss <- ped$sire
#dd <- ped$dam
#nself <- ped$nself

#fMeu <- meuwSelfing(ss, dd, nself)

#fRec <- c()
#for (i in 1:nrow(ped)) {
#    fRec <- c(fRec, cffa(ped, i, i) - 1) 
#}

#fMeu - fRec # dan iguales!!!
