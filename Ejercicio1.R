#Generador congruencial linear
#nsim= numero de simulaciones 
GCL <- function (nsim, semilla=21849, incremento=1, multiplicador=2456651, M=2^32) {
  seq<- semilla 
  for (i in 1:nsim) {
    new <- (seq [i]*multiplicador + incremento) %% M
    seq <-c(seq, new)
    
  }
  return (seq/M)
}

#Se le pueden poner algunos valores por default para ocupar menos parametros al 
#llamar la funcion
x <-GCL (100)
hist (x)
#Ahora comprobamos la independencia entre valores generados, es decir que en 
#secuencia el siguiente no dependa del anterior 
#El signo -, te dice quitale ese indice a mi vector
anterior <-x[-length (x)] #x[1:(length(x)-1)]
siguiente <- x[-1] #x[2:length(x)]
#Graficas en el eje x el numero anterior y en el y el siguiente para ver si existe
#alguna correlacion entre los valores 
plot(anterior, siguiente)

