#FLIP --------------------------------------------------------------------
FLIP <- function(Datos){
  Datos <- as.matrix(Datos) #Para poder operarlo ahora
  N <- nrow(Datos)
  N
  A <- matrix(data=NA, nrow=N, ncol=1)
  A[1,1] <- Datos[N,1]
  head(A)
  for (i in 1:(N-1)){
    
    A[(i+1),1] <- Datos[(N-i),1]
  }
  return(A)
}

# LAG ---------------------------------------------------------------------
LAG <- function(vector){
  tt <- length(vector)
  matt <- matrix(NA,ncol=2,nrow = (tt-1))
  matt[,1] <- vector[2:tt]
    matt[,2] <- vector[1:(tt-1)]
    return(matt)
}

##DIFF ------------------------------------------------------------------
DIFF <- function(vector){
  N <- length(vector)
  diferencia <- c()
  for (i in 2:N) {
    diferencia[i] <- vector[i] - vector[(i-1)]
  }
  return(diferencia[-1])
}

#MEDIA MOVIL---------------------------------------------------------------
MEDIAM <- function(vector, ventana){
  N <- length(vector)
  N
  mm <- matrix(data = NA,nrow = N,ncol = 1)
  for (i in ventana:N) {
    mm[i] <- mean(vector[(i-ventana+1):i])
    
  }
  return(mm)
}

#ACF ----------------------------------------------------------------
ACF <- function(SERIE,Kmax,alpha=0.05){
  Tt <- length(SERIE)
  Media <- mean(SERIE)
  varianza <- sum((SERIE-Media)^2)/(Tt-1)
  
  gamma <- matrix(data=NA, nrow=(Kmax),ncol=1)
  for(K in 0:Kmax){
    gamma[K] <- (sum((SERIE[1:(Tt-K)]-Media)*(SERIE[(K+1):Tt]-Media)))/(Tt-1)
  }
  
  ACFm <- matrix(data=NA, ncol = 1,nrow = (Kmax))
  for(h in 1:Kmax){
    ACFm[h] <- gamma[h]/varianza
      ACFr<- as.matrix(c(1,ACFm))
  }
  #BANDAS
  ACFr <- ACFr[-1]
  vc <- qnorm(1-(alpha/2))
  Q <- length(ACFr)
  bandas <- matrix(NA,ncol = 3, nrow = Q)
  colnames(bandas) <- c("abajo", "arriba","significancia")
  bandas[1,1] <- ACFr[1]-vc*(1/Tt)
  bandas[1,2] <- ACFr[1]+vc*(1/Tt)
  div <- 1/(Tt)
  
  for (q in 2:Q) {
    bandas[q,1] <- -2*sqrt(div*(1+2*sum(ACFr[1:q]^2)))
    bandas[q,2] <- 2*sqrt(div*(1+2*sum(ACFr[1:q]^2)))
    
    if(ACFr[q]<bandas[q,1]|ACFr[q]>bandas[q,2]){
      bandas[q,3] <- 1
    } else {
      bandas[q,3] <- 0
    }
    
  }
  acfbandas <- cbind(ACFr,bandas)
  # plot(acfbandas[,1], type = "h",ylim = c(min(acfbandas[,2]),max(acfbandas[,3])))
  # lines(acfbandas[,2],lty="dashed")
  # lines(acfbandas[,3],lty="dashed")
  # abline(h=0) #Esta gráfica empieza desde rho1
  
return(ACFr)
}


#PACF -----------------------------------------------------------------------------
PACF <- function(SERIE,Kmax) {
  
  #ACF------------------------------------------------------------
  Tt <- length(SERIE)
  Media <- mean(SERIE)
  varianza <- sum((SERIE-Media)^2)/(Tt-1)
  
  gamma <- matrix(data=NA, nrow=(Kmax),ncol=1)
  for(K in 0:Kmax){
    gamma[K] <- (sum((SERIE[1:(Tt-K)]-Media)*(SERIE[(K+1):Tt]-Media)))/(Tt-1)
  }
  head(gamma)
  dim(gamma)
  
  ACFm <- matrix(data=NA, ncol = 1,nrow = (Kmax))
  for(h in 1:Kmax){
    ACFm[h] <- gamma[h]/varianza
    ACFr<- as.matrix(c(1,ACFm))
  }
  
    #PACF--------------------------------------------------------
  PACFm <- matrix(data=NA, nrow=Kmax, ncol=1)
  for (h in 1:Kmax) {
    denominador <- matrix(data=NA, ncol=h, nrow=h)
    for (i in 1:h) {
      for (j in 1:h) {
        denominador[i,j] <- ACFr[abs(i-j)+1]
      }
    }
    numerador <- denominador[,-j]
    numerador <- cbind(numerador,ACFr[2:(i+1)])
    
    phik <- det(numerador)/det(denominador)
    PACFm[h] <- as.numeric(phik)
  }
  
  #BANDAS
  plot(x=c(1:Kmax),y=PACFm,type="h",
       xlab="k",ylab="phikk",main="pacf", ylim = c((-qnorm(0.975)/sqrt(Tt)),(qnorm(0.975)/sqrt(Tt))))
  abline(a=0,b=0)
  abline(qnorm(0.975)/sqrt(Tt),0,col="blue",lty="dashed")
  abline(-qnorm(0.975)/sqrt(Tt),0,col="blue",lty="dashed") #Esta gráfica empieza desde rho1
  
  return(PACFm)
}


# DESCOMPOSICIÓN DE SERIE -------------------------------------------------
DESCOMPONER <- function(SERIE,K){
  library(dplyr)
  #Tendencia de la serie----------------------------------------------
  datos <- data.frame(serie= SERIE,MM= MEDIAM(SERIE,K))
  datos <-  datos%>% mutate(sin_tendencia= (serie-MM))
  datos <-  datos[-c(1:K),]
  #Componente estacional----------------------------------------------
  n <- trunc(nrow(datos)/K)
  n # Número de Años
  
  #Creo matriz que me contenga los datos sin tendendia por mes
  anualxmes <- matrix(datos$sin_tendencia,nrow = n,ncol = K,byrow = T) #Tengo n datos, y K=12 meses
  
  #Vector que va a tener la media de cada mes
  mediames <- c()
  for (j in 1:K) {
    mediames[j] <- (sum(anualxmes[,j]))/(n) #lo indexo en j porque voy a tener K número de medias
  }
  
  #Sj= Ej_barra - E_barra
  ss <- c()
  for (j in 1:K) {
    ss[j] <- mediames[j]-mean(datos$sin_tendencia)
  }
  #Meterlo en el data frame
  estacc <- rep(x=ss,times= n) #Para que quede la estacionalidad para cada año
  datos <- datos%>% mutate(estacionalidad=estacc)
  
  #Componente Aleatorio-----------------------------------------------
  datos <- datos%>% mutate(innovaciones= (sin_tendencia-estacionalidad))
  datos <- datos%>% mutate(desestacionalizada= (serie-estacionalidad))
  par(mfrow= c(5,1),mar = c(1, 1, 1, 1))
  plot(datos$serie,type="l", main= "serie")
  plot(datos$sin_tendencia,type="l", main= "sin tendencia")
  plot(datos$estacionalidad,type="l", main= "estacionalidad")
  plot(datos$innovaciones,type="l" , main= "innovaciones")
  plot(datos$desestacionalizada,type="l", main= "desestacionalizada")
  
  return(datos)
  
}


# ESTABILIZACIÓN VARIANZA -------------------------------------------------
E_VARIANZA <- function(datos,obsxgrupo){
  #Parámetros para acomodar serie ---------------------------------------
  N <- length(datos) #Número de observaciones total 
  R <- obsxgrupo   #Número de observaciones por grupo 
  H <- trunc(N/R)  #Número de grupos 
  n <- N- H*R      #Número de observaciones excluidas
  
  datos <- datos[1:(N-n)]
  #Medias para cada grupo -----------------------------------------------
  datosmatriz <- matrix(data=datos, ncol = H, nrow = R, byrow = FALSE)
  
  medias <- matrix(data = NA, nrow = 1, ncol = H)
  for (h in 1:H) {
    medias[h] <- mean(datosmatriz[,h])
  }
  
  #Desviación estándar para cada grupo ------------------------------------
  desviaciones <- matrix(data = NA, nrow=1, ncol= H)
  for (h in 1:H) {
    desviaciones[h] <- 
      sqrt(sum((datosmatriz[,h]- medias[h])^2)/(R-1))
  }
  
  #Coeficientes de variación dependiendo de cada lambda ----------------------
  lambdas <- c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)
  L <- length(lambdas)
  Mlambda <- matrix(data=NA, nrow = 1, ncol = L)
  SDlambda <-  matrix(data = NA, nrow = 1, ncol = L)
  for (l in 1:L) {
    cvpotencia <- matrix(data = NA, nrow = 1, ncol = H)
    for (h in 1:H) {
      cvpotencia[h] <- 
        desviaciones[h]/ (medias[h]^(1-lambdas[l]))
    }
    Mlambda[l] <- sum(cvpotencia)/H
    SDlambda[l] <- sqrt(sum((cvpotencia- Mlambda[l])^2)/(H-1))
  }
  
  #En el siguiente código vamos a tener los CV para los lambdas y no para los grupos
  CVlambda <- matrix(data = NA, ncol = L, nrow = 1)
  for (l in 1:L) {
    CVlambda[l] <- SDlambda[l]/Mlambda[l]  
  }
  
  #Escoger lambda con menor coeficiente de variación -------------------------
  indice <- which.min(CVlambda)
  lambdacv <- lambdas[indice]
  texto <- print(paste("el lambda que debería usar es:",lambdacv))
  return(lambdacv)
}

# ESTABILIZACIÓN NIVEL  ---------------------------------------------------
 #funcioncita para casa S
SDNIVEL <- function(transformada,j=0){
  N <- length(transformada)
  
  primero <- 1/(N-1)
  suma2 <- sum(transformada/(N))
  sumsum <- sum((transformada - suma2)^2)
  ss <- primero*sumsum
  so <- sqrt(ss)
  return(so)
}
# Iteración para cada S y generación de gráficos ACF
NIVEL <- function(datos,J=3){
  matricita <- matrix(NA, nrow= 4, ncol= 1)
  data <- datos
  par(mfrow= c((J+1),2),mar = c(1, 1, 1, 1))
  for (j in 0:J) {
    if(j>0){
      data <- DIFF(data)
    }
    plot(data, type="l", main= "Serie con diferencia 0")
    plot(c(ACF(data,24)), type = "h", main= paste("serie con dif de", j))
    desves <- SDNIVEL(data,j)
    matricita[(j+1)] <- desves
  }
  nv <- which.min(matricita)-1
  texto <- print(paste("el nivel más apropiado es:",nv))
  return(nv)
}



 ##### COMPROBACIÓN DE SUPUESTOS------------------------------
# PRUEBA DE MEDIA CERO ----------------------------------------------------
MEDIA0 <- function(choques,d,p,q){
  Tt <- length(choques)
  tprima <- d+p+1
  
  choqnuevos <- choques[tprima:Tt]
  tnuevo <- length(choqnuevos)
  media <- sum(choqnuevos)/(tnuevo)
  
  numerador <- sum((choqnuevos- media)^2)
  desves <- sqrt(numerador/(tnuevo-q))
  
  decision <- abs((sqrt(tnuevo-d)* media)/desves)
  
  if (decision<2) {
    print("no existe evidencia para rechazar H0, media=0 estadísticamente")
  } else  {
    print("La media no es 0")
  }
  return(decision)
}


#RUIDO BLANCO
# LJUNG BOX ---------------------------------------------------------------
LJUNGBOX <- function(choques,d,p,q,alpha){
  tt <- length(choques)
  ka <- tt/4
  
  mult1 <- (tt-d-p)
  mult2 <- (tt-d-p+2)
  denom <- (tt-d-p-ka)
  
  rhos <- ACF(choques,ka,0.05)
  rhoss <- rhos[-1]
  
  LB <- mult1*mult2* sum((rhoss^2)/denom)
  
  VC <- qchisq((1-alpha),(ka-p-q))
  #if(LB<VC){
    #print("No se rechaza ruido blanco")
  #}else{print("Se rechaza ruido blanco")
  #}
  ll <- cbind(LB,VC)
  return(ll)
  }
# BOX-PIERCE --------------------------------------------------------------
BOXPIERCE <- function(choques,Tt,d,p,ka,alpha=0.05){
  per <- Tt-d-p
  acfchoq <- ACF(choques,(ka+1))
  BP <- per* sum((acfchoq[2:(ka+1)])^2)
  
  vc <- qchisq((1-alpha),(ka-p-q))
  if(BP<vc){
    print("No se rechaza ruido blanco")
  }else{print("Se rechaza ruido blanco")
  }
  return(BP)
  
}

# PRUEBA DE NORMALIDAD ----------------------------------------------------
NORMALIDAD <- function(choques,d,p){
  Tt <- length(choques)
  tprima <- d+p+1
  
  choqnuevos <- choques[tprima:Tt]
  tnuevo <- length(choqnuevos)
  media <- sum(choqnuevos)/(tnuevo)
  
  numerador <- sum((choqnuevos- media)^2)
  desves <- sqrt(numerador/(tnuevo-q))

maxfuera <- ceiling((Tt-d-p)/20)
limsup <- 2*desves
liminf <- -2*desves
menoresa<- sum(choques<liminf)
mayoresa <- sum(choques>limsup)
porfueratotal <- menoresa + mayoresa

if(porfueratotal> maxfuera){
  print("No hay normalidad")
} else {
  print("Hay normalidad")
}
compara <- matrix(cbind(maxfuera, porfueratotal),ncol = 2)
colnames(compara) <- c("maxfuera","porfueratotal")
return(compara)
}

# JARQUE BERA -------------------------------------------------------------
JARQUEBERA <- function(vector,d){
#CÁLCULO SIMETRÍA
Tt <- length(vector)
Tt
#Sacamos la media para ajustar las medidas 
media <- mean(vector)
#Primero vamos a sacar sigma cuadrado y los momentos 3 y 4 
sigma2 <- sum((vector-media)^2)/(Tt-d)
Desviacion<- sqrt(sigma2)

#momento 3
m3<- sum((vector-media)^3)/(Tt-d)
#momento 4 
m4<- sum((vector-media)^4)/(Tt-d)
#Simetría 
Simetria<- m3/(Desviacion^3)
Simetria

#CALCULO DE EXCESO DE CURTOSIS
Ex_curtosis<- m4/(Desviacion^4)
Ex_curtosis

#CÁLCULO DEL TEST DE JARQUE BERA
Jarque_Bera <- (Tt-d)*(((Simetria^2)/6)+((Ex_curtosis^2)/24))
Jarque_Bera

VC <- qchisq(1-0.05,2)
if(Jarque_Bera< VC){
  print("no se rechaza normalidad")
} else { 
  print("Se rechaza normalidad")
  }

return(Jarque_Bera)
}

# DATOS ABERRANTES --------------------------------------------------------
D_ABERRANTES <- function(choques,d,p,q){
  Tt <- length(choques)
  tprima <- d+p+1
  
  choqnuevos <- choques[tprima:Tt]
  tnuevo <- length(choqnuevos)
  media <- sum(choqnuevos)/(tnuevo)
  
  numerador <- sum((choqnuevos- media)^2)
  desves <- sqrt(numerador/(tnuevo-q))

  limsup <- 3*desves
  liminf <- -3*desves
  
  dummies <- matrix(NA, ncol = 1, nrow = Tt)
    for (t in 1:Tt) {
      if(
      choques[t]<= liminf | choques[t]>= limsup){
        dummies[t] <- 1
      } else{
        dummies[t] <- 0
      }
    }
  plot(dummies, type= "h")
  return(dummies)
  }

# AR(P)  -------------------------------------------------------------------
ARP <- function(N,p,phis,const=0){
  N <- N+600
  phis <- as.vector(phis)
  p <- length(phis)
  
  vo <- c()
  ar <- matrix(NA,ncol=1,nrow=N)
  e <- rnorm(N,mean = 0,sd = 1)
  for (i in 1:p) {
    vo[i] <- runif(1,min = 1,max = 20)
    ar[i,] <- vo[i]
    
    for (t in (p+1):N) {
      ar[t] <- sum(phis[i]*ar[(t-i)]) + e[t]
    }
  }
  arfin <- ar[601:N]
  return(arfin)
}

#----------------------------------SEGUNDO CORTE-----------------------------------
#DGPS ----------------- 
DGPS <- function(tt,nsimul,valor0,miu,beta,phi,sdproceso){
  
  tusar <- tt+500
  #Matrices vacías de choques y dgps con columnas iguales al # de simulaciones
  choques <- matrix(NA, nrow = tusar,ncol=1)
  dgp <- matrix(NA, nrow = tusar,ncol = 1)
  dgp[1,] <- valor0
  dgpsimul <- matrix(NA, nrow = tusar,ncol = nsimul)
  
  #Tendencia
  contador <- c(1:tusar)
  
  for (s in 1:nsimul) {
    #choques 
    choques[,1] <- rnorm(n = tusar,mean = 0,sd = sdproceso)
    #for para número de simulaciones que quiero
    #dgps
    
    for (t in 2:tusar) {
      dgp[t,1] <- miu + beta*contador[t] + phi*dgp[(t-1),1] + choques[t,1]
    }
    #Llenar por # de simulaciones
    dgpsimul[,s] <- dgp
  }
  dgpsimul <- dgpsimul[-(1:500),]
  return(dgpsimul)
}

# MATRIZ AUMENTADA SIN INTERCEPTO NI PENDIENTE ----------------------------
MXDELT <- function(y,p){
  y <- as.vector(y)
  N <- length(y)
  deltay <- DIFF(y)
  yn <- y[((p+1):N-1)] #Primer fila de la matriz de diseño
  if(p==1){
    x <- matrix(yn,ncol=1, nrow=(N-p))
  } else {
  maumentada <- matrix(data=NA, ncol=(p-1),nrow= (N-p)) #Para luego hacerle un cbind con yn
  #Si tomamos N como el tamaño de la función sin diferenciar: 
  #En la fila de gamma(fila con y_t-1) pierdo la primera posición porque no tengo momento 0, luego en la primera posición de la parte aumentada pierdo la diferencia de la posicion t-1 y la de la posición t, así voy a tener N-1 en la parte sin estar aumentada y N-2 en la parte con P=1, o sea que mi dimensión, en general, va a ser N-p-1
  
  for (i in 1:(p-1)) {
    maumentada[,i] <- deltay[((p-i):(N-i-1))]
    #Yo tengo que ponerle el -1 porque mi vector dif viene sin el NA
  }
  
  x <- cbind(yn,maumentada)
  }
  return(x)
}

# REGRESIÓN LINEAR --------------------------------------------------------
REGRESION <- function(x,y){
  N <- nrow(x)
  k <- ncol(x)
  coef <- solve(t(x)%*%x)%*%t(x)%*%y
  ystim <- x%*%coef
  e <- y-ystim
  sigma <- (t(e)%*%e)/(N-k)
  varcov <- as.numeric(sigma)*solve(t(x)%*%x)
  
  lista <- list(ndatos=N,k=k,coef=coef,ystim=ystim,e=e,sigma=sigma,varcov=varcov)
  return(lista)
}

# TAU ADF -----------------------------------------------------------------
TAUADF <- function(y,p,esq){
  y <- as.matrix(y)
  N <- length(y)
  deltay <- DIFF(y)
  dim(y)
  #Inicio de condicionales
  if (esq==1){
    if(p==1){
      yusar <- deltay[p:(N-1)]
      yn <- y[((p+1):N-1)]
      x <- cbind(yn)
      yusar <- deltay[p:(N-1)]
      
      coef <- REGRESION(x,yusar)$coef
      varcov <- REGRESION(x,yusar)$varcov
      desves <- sqrt(varcov[1,1])
      tau <- coef[1]/desves
      
    } else{
      maumentada <- matrix(data=NA, ncol=(p-1), nrow=N-p) #Para luego hacerle un cbind con yn
      yn <- y[((p+1):N-1)] #Primer fila de la matriz de diseño
      #Si tomamos N como el tamaño de la función sin diferenciar: 
      #En la fila de gamma(fila con y_t-1) pierdo la primera posición porque no tengo momento 0, luego en la primera posición de la parte aumentada pierdo la diferencia de la posicion t-1 y la de la posición t, así voy a tener N-1 en la parte sin estar aumentada y N-2 en la parte con P=1, o sea que mi dimensión, en general, va a ser N-p-1
      
      for (i in 1:(p-1)) {
        maumentada[,i] <- deltay[((p-i):(N-i-1))]
        #Yo tengo que ponerle el -1 porque mi vector dif viene sin el NA
      }
      x <- cbind(yn,maumentada)
      yusar <- deltay[p:(N-1)]
      
      coef <- REGRESION(x,yusar)$coef
      varcov <- REGRESION(x,yusar)$varcov
      desves <- sqrt(varcov[1,1])
      tau <- coef[1]/desves
    }
    
  } else {
    if(esq==2){
      if(p==1){
        yn <- y[((p+1):N-1)]
        x <- cbind(1,yn)
        yusar <- deltay[p:(N-1)]
        
        coef <- REGRESION(x,yusar)$coef
        varcov <- REGRESION(x,yusar)$varcov
        desves <- sqrt(varcov[2,2])
        tau <- coef[1]/desves
        
      } else{
        maumentada <- matrix(data=NA, ncol=(p-1), nrow=N-p) 
        yn <- y[((p+1):N-1)] #Primer fila de la matriz de diseño
        
        for (i in 1:(p-1)) {
          maumentada[,i] <- deltay[((p-i):(N-i-1))]
        }
        
        x <- cbind(1,yn,maumentada)
        yusar <- deltay[p:(N-1)] 
        
        coef <- REGRESION(x,yusar)$coef
        varcov <- REGRESION(x,yusar)$varcov
        desves <- sqrt(varcov[2,2])
        tau <- coef[2]/desves
      }
      
    } else {
      if (esq==3){
        if(p==1){
          contador <- c((p+1):N)
          yusar <- deltay[p:(N-1)]
          yn <- y[((p+1):N-1)]
          x <- cbind(1,contador,yn)
          
          coef <- REGRESION(x,yusar)$coef
          varcov <- REGRESION(x,yusar)$varcov
          desves <- sqrt(varcov[3,3])
          tau <- coef[3]/desves
          
        } else{
          
          maumentada <- matrix(data=NA, ncol=(p-1), nrow=N-p) 
          yn <- y[((p+1):N-1)] #Primer fila de la matriz de diseño
          
          for (i in 1:(p-1)) {
            maumentada[,i] <- deltay[((p-i):(N-i-1))]
          }
          contador <- c((p+1):N)
          x <- cbind(1,contador,y[(p:(N-1))],maumentada)
          yusar <- deltay[p:(N-1)]
          
          coef <- REGRESION(x,yusar)$coef
          varcov <- REGRESION(x,yusar)$varcov
          desves <- sqrt(varcov[3,3])
          tau <- coef[3]/desves
        }
        
      } else{
        tau <- print("Este modelo no está dentro de los esquemas usados")
      }
    }
  } #Fin primer else 
  return(tau)
}

# PHIS ADF ----------------------------------------------------------------

PHISADF <- function(y,p,esq){
  y <- as.matrix(y)
  n <- length(y)
  deltay <- DIFF(y)
  #Inicio de condicionales
  #PHI1
  if (esq==2){
    
    # SIN PARTE AUMENTADA -----------------------------------------------------
    if(p==1){
      yn <- y[(p:(n-1))]
      x <- cbind(1,yn)
      yusar <- deltay[p:(n-1)]
      regre <- REGRESION(x,yusar)
      coef <- regre$coef
      yestim <- regre$ystim
      e <- regre$e
      #Suma de residuales no restringidos
      SRCNR <- (t(e)%*%e)
      sigma <- SRCNR/(n-p-1)
      #Suma de residuales restringida 
      SRCR <- (t(yusar)%*%yusar)
      #Estadístico
      phi <- (SRCR-SRCNR)/(sigma*2)
    }  
    # CON PARTE AUMENTADA -----------------------------------------------------
    else{
      yn <- y[(p:(n-1))] #Primer fila de la matriz de diseño
      mdeltas <- matrix(data=NA, ncol=p-1, nrow=n-p) #Para luego hacerle un cbind con yn
      for (c in 1:(p-1)) {
        mdeltas[,c] <- deltay[(p-c):(n-c-1)] #Segunda parte de la 
      }
      x <- cbind(1,yn,mdeltas)
      yusar <- deltay[p:(n-1)]
      regre <- REGRESION(x,yusar)
      coef <- regre$coef
      yestim <- regre$ystim
      e <- regre$e
      #Suma de residuales no restringidos
      SRCNR <- (t(e)%*%e)
      sigma <- SRCNR/(n-p-1)
      #Suma de residuales restringida 
      regre2 <- REGRESION(mdeltas,yusar)
      coefr <- regre2$coef
      yestimr <- regre2$ystim
      er <- regre2$e
      SRCR <- (t(er)%*%er)
      #Estadístico
      phi <- (SRCR-SRCNR)/(sigma*2)
      
    }
  } 
  
  # PHI 2 y 3 -------------------------------------------------------------------
  else {
    #PHI2
    if(esq==3){
      
      # SIN PARTE AUMENTADA -----------------------------------------------------
      if(p==1){
        contador <- c(p:(n-1))
        
        x <- cbind(1,contador,y[(p:(n-1))])
        yusar <- deltay[p:(n-1)]
        
        regre <- REGRESION(x,yusar)
        coef <- regre$coef
        yestim <- regre$ystim
        e <- regre$e
        #Suma de residuales no restringidos
        SRCNR <- (t(e)%*%e)
        sigma <- SRCNR/(n-p-3)
        #Suma de residuales restringida
        SRCR <- (t(yusar)%*%yusar)
        #Estadístico
        phi2 <- (SRCR-SRCNR)/(sigma*3)
        
        #PHI3 
        #Suma de residuales restringida 
        erphi3 <- yusar-coef[3]
        SRCR3 <- (t(erphi3)%*%erphi3)
        #Estadístico
        phi3 <- (SRCR3-SRCNR)/(sigma*2)
        phi3
        phi <- cbind(phi2,phi3)
        colnames(phi) <- c("phi2","phi3")
      }
      
      # CON PARTE AUMENTADA -----------------------------------------------------
      
      else{
        yn <- y[(p:(n-1))] #Primer fila de la matriz de diseño
        mdeltas <- matrix(data=NA, ncol=p-1, nrow=n-p) #Para luego hacerle un cbind con yn
        for (c in 1:(p-1)) {
          mdeltas[,c] <- deltay[(p-c):(n-c-1)] 
        }
        contador <- c(p:(n-1))
        
        x <- cbind(1,contador,y[(p:(n-1))],mdeltas)
        yusar <- deltay[p:(n-1)]
        
        regre <- REGRESION(x,yusar)
        coef <- regre$coef
        yestim <- regre$ystim
        e <- regre$e
        #Suma de residuales no restringidos
        SRCNR <- (t(e)%*%e)
        sigma <- SRCNR/(n-p-3)
        #Suma de residuales restringida
        regre2 <- REGRESION(mdeltas,yusar)
        coefr <- regre2$coef
        yestimr <- regre2$ystim
        er <- regre2$e
        SRCR <- (t(er)%*%er)
        #Estadístico
        phi2 <- (SRCR-SRCNR)/(sigma*3)
        
        #PHI3 
        #Suma de residuales restringida 
        xr <- cbind(1,mdeltas)
        coeff3 <- solve(t(xr)%*%xr)%*%t(xr)%*%yusar
        yestif3 <- xr%*%coeff3
        er3 <- yusar-yestif3
        SRCR3 <- (t(er3)%*%er3)
        #Estadístico
        phi3 <- (SRCR3-SRCNR)/(sigma*2)
        phi3
        phi <- cbind(phi2,phi3)
        colnames(phi) <- c("phi2","phi3")
      }
      
    } else{
      phi <- print("Este modelo no está dentro de los esquemas usados")
    }
  }
  return(phi)
}  

## DICKEY FULLER----------------------------------

## P OPTIMO -------------------------------------------------------------
P_ESQUEMA3 <- function(y,alpha=0.05){
  
  N <- length(y)
  p <- N/4
  
  #Esquema sin parte aumentada 
  #Matriz diseño esquema 3 
  delt1 <- MXDELT(y,1)
  x1 <- cbind(1,c((2):N),delt1)
  #Coeficientes
  deltay1 <- DIFF(y)
  yusar1 <- deltay1[1:(N-1)]
  
  regre1 <- REGRESION(x1,yusar1)
  coef1 <- regre1$coef
  yestim1 <- regre1$ystim
  e1 <- regre1$e
  sigma1 <- regre1$sigma
  varcov1 <- regre1$varcov
  desves1 <- sqrt(varcov1[(3),(3)])
  ##Ruido Blanco
  lb1 <- LJUNGBOX(e1,0,1,0,alpha)
  lb11 <- lb1[1]
  vclb1 <- lb1[2]
  
  if(lb11<vclb1){
    pusar <- 1
    listap <- list(p=pusar, tabla=1, ystim=yestim1,e=e1,varcov=varcov1,desves=desves1, coef=coef1)
  } else{
    
    # CREACIÓN DE LA TABLA PARA PARTE AUMENTADA
    tablita <- matrix(NA,ncol = 3,nrow = p)
    row.names(tablita) <- c(1:p)
    colnames(tablita) <- c("Signif","RB","Choose")
    #SIGNIFICANCIA Y RUIDO BLANCO -----------------------------------------------------
    ##Significancia ----------------------------------
    for (i in 2:p) {
      
      #Matriz diseño esquema 3 
      delt <- MXDELT(y,i)
      x <- cbind(1,c((i+1):N),delt)
      #Coeficientes
      deltay <- DIFF(y)
      yusar <- deltay[i:(N-1)]
      
      regre <- REGRESION(x,yusar)
      coef <- regre$coef
      yestim <- regre$ystim
      e <- regre$e
      sigma <- regre$sigma
      varcov <- regre$varcov
      desves <- sqrt(varcov[(i+2),(i+2)])
      tsig <- coef[(i+2)]/desves
      vct <- qt((1-(alpha/2)),(N-(i+2)-i))
      ##Ruido Blanco
      lb <- LJUNGBOX(e,0,i,0,alpha)
      lb1 <- lb[1]
      vclb <- lb[2]
      #Tabla de decisión ------------ 
      
      #Condicional para significancia
      if(tsig>vct|tsig<(-vct)){
        tablita[i,1] <- 1
      } else{
        tablita[i,1] <- 0
      }
      
      #Condicional para ruido blanco
      
      if(lb1<vclb){
        tablita[i,2] <- 1
      } else{
        tablita[i,2] <- 0
      }   
      
      #Condicional decisión
      if((tablita[i,2]+ tablita[i,1])==2){
        tablita[i,3] <- 1
      } else{
        tablita[i,3] <- 0
      }
      
      #Condicional para terminar iteración
      if(tablita[i,3]==1){
        break
      }
      
    } #Cierre primer for
    
    tablita <- na.omit(tablita)
    tablita
    if(p>1){
      pusar <- nrow(tablita)+1
    }
    
    listap <- list(p=pusar, tabla=tablita, ystim=yestim,e=e,varcov=varcov,desves=desves, coef=coef)
  }
  return(listap)
  
}


## PRUEBA DE RAIZ UNIRTARIA SOBRE ESQUEMA 3 ---------------
RAIZESQUEMA3 <- function(y,p){
  N <-  length(y)
  #Para tau rechazo si el estadístico es menor al vc 
  #Para phis si phi es mayor al vc rechazo 
  tautau <- TAUADF(y,p,3)
  t_a2 <- coef[2]/(sqrt(varcov[2,2]))
  vct_a2 <- qt((1-(0.05/2)),(N-p-nrow(coef)))
  
  if(tautau<vctautau){ #Rechazo al 5%
    # Comprobamos si es DGP 1 O DGP2
    if(abs(t_a2)<abs(vct_a2)){ #No rechazo al 5%
      resultado <- print(paste("Con una prueba tautau de", tautau,"y con un estadístico individual de",t_a2, "se rechaza la hipótesis nula de gamma= 0, pero no se rechaza que y a2= 0, se debe usar un DGP1"))
    }else{ 
      resultado <- print(paste("Con una prueba tautau de", tautau,"y con un estadístico individual de",t_a2, "existe evidencia para rechazar que gamma y a2 son iguales a 0, se debe usar un DGP2"))
    }
  } else {
    resultado <- print(paste("A un nivel de significancia de 5%, un estadístico tautau de", tautau,"Hay evidencia para aceptar H0, por lo que se dice que hay raíz unitaria. Pasar a la prueba phi3"))
  }
  
  milista <- list(resultado= resultado, tautau= tautau)
  return(milista)
}

## PHI3 ----------------
PHI3 <- function(y,p){
  N <-  length(y)
  #Para phis si phi es mayor al vc rechazo 
  phi3 <- PHISADF(y,p,3)[2]
  
  #Parte individual para a2
  #Nueva regresión 
  if(p==1){
    x3 <- cbind(1,c((2):N))
    #Coeficientes
    deltyn3 <- DIFF(y)
    yn3 <- deltyn3[1:(N-1)]
    
    regre3 <- REGRESION(x3,yn3)
    coef3 <- regre3$coef[2]
    varcov3 <- regre3$varcov
    t_a2 <- coef3/(sqrt(varcov3[2,2]))
    vct_a2 <- qt((1-(0.05/2)),(N-p-ncol(x3)))
  }else{
    delx3 <- MXDELT(y,p)
    delx3 <- delx3[,-1]
    x3 <- cbind(1,c((p+1):N),delx3)
    #Coeficientes
    deltay3 <- DIFF(y)
    yusar3 <- deltay3[p:(N-1)]
    
    regre3 <- REGRESION(x3,yusar3)
    coef3 <- regre3$coef
    varcov3 <- regre3$varcov
    t_a2 <- coef3[2]/(sqrt(varcov[2,2]))
    vct_a2 <- qt((1-(0.05/2)),(N-p-ncol(x3)))
  }
  
  if(phi3<6.25){ #No Rechazo al 5%
    resultado <- print(paste("Con un estadístico Phi3 de", phi3,"No hay evidencia para rechazar la hipótesis nula, gamma y a2 son iguales a 0, bajar al esquema 2"))
  } else{
    #Prueba de significancia individual
    
    if(abs(t_a2)<abs(vct_a2)){ #No rechazo al 5%
      resultado <-  print(paste("Con un estadístico Phi3 de", phi3," y una prueba individual de", t_a2, "Hay evidencia para aceptar que gamma=0, pero no hay evidencia par rechazar la hipótesis del estadístico individual y a2=0, bajar al esquema 2"))
    }else{ 
      resultado <- print(paste("Con un estadístico Phi3 de", phi3,"y una prueba individual de",t_a2, "Hay evidencia para aceptar que gamma=0, pero a2 es diferente de 0, entonces se debe usar un DGP2"))
    }
  }
  milista <- list(resultado= resultado, phi3= phi3)
  return(milista)
}

# Esquema 2 --------------------------------

#P optimo-----------------
P_ESQUEMA2 <- function(y,alpha=0.05){
  
  N <- length(y)
  p <- N/4
  
  #Esquema sin parte aumentada 
  #Matriz diseño esquema 3 
  delt1 <- MXDELT(y,1)
  x1 <- cbind(1,delt1)
  #Coeficientes
  deltay1 <- DIFF(y)
  yusar1 <- deltay1[1:(N-1)]
  
  regre1 <- REGRESION(x1,yusar1)
  coef1 <- regre1$coef
  yestim1 <- regre1$ystim
  e1 <- regre1$e
  sigma1 <- regre1$sigma
  varcov1 <- regre1$varcov
  ##Ruido Blanco
  lb1 <- LJUNGBOX(e1,0,1,0,alpha)
  lb11 <- lb1[1]
  vclb1 <- lb1[2]
  
  if(lb11<vclb1){
    pusar <- 1
    listap <- list(p=pusar, tabla=1, ystim=yestim1,e=e1,varcov=varcov1, coef=coef1)
  } else{
    tablita <- matrix(NA,ncol = 3,nrow = p)
    row.names(tablita) <- c(1:p)
    colnames(tablita) <- c("Signif","RB","Choose")
    #SIGNIFICANCIA Y RUIDO BLANCO -----------------------------------------------------
    ##Significancia ----------------------------------
    for (i in 2:p) {
      
      #Matriz diseño esquema 3 
      delt <- MXDELT(y,i)
      x <- cbind(1,delt)
      #Coeficientes
      deltay <- DIFF(y)
      yusar <- deltay[i:(N-1)]
      
      regre <- REGRESION(x,yusar)
      coef <- regre$coef
      yestim <- regre$ystim
      e <- regre$e
      sigma <- regre$sigma
      varcov <- regre$varcov
      desves <- sqrt(varcov[(i+1),(i+1)])
      tsig <- coef[(i+1)]/desves
      vct <- qt((1-(alpha/2)),(N-(i+1)-i))
      ##Ruido Blanco
      lb <- LJUNGBOX(e,0,i,0,alpha)
      lb1 <- lb[1]
      vclb <- lb[2]
      #Tabla de decisión ------------ 
      
      #Condicional para significancia
      if(tsig>vct|tsig<(-vct)){
        tablita[i,1] <- 1
      } else{
        tablita[i,1] <- 0
      }
      
      #Condicional para ruido blanco
      
      if(lb1<vclb){
        tablita[i,2] <- 1
      } else{
        tablita[i,2] <- 0
      }   
      
      #Condicional decisión
      if((tablita[i,2]+ tablita[i,1])==2){
        tablita[i,3] <- 1
      } else{
        tablita[i,3] <- 0
      }
      
      #Condicional para terminar iteración
      if(tablita[i,3]==1){
        break
      }
      
    } #Cierre primer for
    
    tablita <- na.omit(tablita)
    tablita
    
    if(p>1){
      pusar <- nrow(tablita)+1
    }
    
    listap <- list(p=pusar, tabla=tablita, ystim=yestim,e=e,varcov=varcov,desves=desves, coef=coef)
  }
  return(listap)
  
}




#Raiz unitaria esquema 2 --------------------
RAIZESQUEMA2<- function(y,p) {
  N <- length(y)
  #Para tau rechazo si el estadístico es menor al vc 
  taumiu <- TAUADF(y,pusar2,2)
  t_a0 <- coef2[1]/sqrt((varcov2[1,1]))
  vca0 <- qt(1-(0.05/2),df = (N-pusar2-nrow(coef2)),lower.tail = F)
  
  if(taumiu<vctaumiu){ #Rechazo al 5%
    
    if(abs(t_a0) < abs(vca0)) {# No Rechazo
      resultado <- print(paste("Con un estadístico t de", t_a0,"Probar el esquema 1"))} else{
        resultado <- print(paste("con un estadístico individual de",t_a0, "existe evidencia para rechazar que a2 es diferente de 0 y se debe usar un DGP1"))
      }
  } else{
    resultado <- print("Hay existencia de raíz unitaria, hacer prueba phi1")
  }
  
  milista <- list(resultado= resultado, taumiu= taumiu, indiv= t_a0,vc= vca0)
  return(milista)
}

#### PRUEBA PHI 1 --------------------
PHI1<- function(y,p){
  N <- length(y)
  #Para tau rechazo si el estadístico es menor al vc 
  #Para phis si phi es mayor al vc rechazo 
  phi1 <- PHISADF(y,p,2)
  
  if(p==1){
    x3 <- matrix(rep(1,(N-1)),ncol = 1)
    #Coeficientes
    deltyn3 <- DIFF(y)
    yn3 <- deltyn3[1:(N-1)] 
    
    regre3 <- REGRESION(x3,yn3)
    coef3 <- regre3$coef
    varcov3 <- regre3$varcov
    t_a0 <- coef3[1]/(sqrt(varcov3[1]))
    vct_a0 <- qt((1-(0.05/2)),(N-p-ncol(x3)))
  }else{
    delx3 <- MXDELT(y,p)
    delx3 <- delx3[,-1]
    x3 <- cbind(1,delx3)
    #Coeficientes
    deltay3 <- DIFF(y)
    yusar3 <- deltay3[p:(N-1)]
    
    regre3 <- REGRESION(x3,yusar3)
    coef3 <- regre3$coef
    varcov3 <- regre3$varcov
    t_a0 <- coef3[1]/(sqrt(varcov3[1,1]))
    vct_a0 <- qt((1-(0.05/2)),(N-p-ncol(x3)))
  }
  
  if(phi1<4.59){ #No Rechazo al 5%
    resultado <- print(paste("Con un estadístico Phi1 de", phi1,"No hay evidencia para rechazar la hipótesis nula, gamma y a0 son iguales a 0, bajar al esquema 1"))
  } else{
    #Prueba de significancia individual
    
    if(abs(t_a0)<abs(vct_a0)) { #No rechazo al 5%
      resultado <- print("No existe evidencia para rechazar que a0=0, bajar al esquema 1")
    }else { 
      resultado <- print(paste("con un phi1 de" ,phi1, "y con un estadístico individual de",t_a0, "no existe evidencia para rechazar la existencia de rapiz unitaria y que a0=0, por tanto se debe usar un DGP4"))
    }
  }
  
  milista <- list(resultado= resultado,phi1=phi1, indiv= t_a0,v=varcov3,vc=vct_a0)
  return(milista)
}

### ESQUEMA 1 

##P optimo 
P_ESQUEMA1 <- function(y,alpha=0.05){
  
  N <- length(y)
  p <- N/4
  
  #Esquema sin parte aumentada 
  #Matriz diseño esquema 3 
  x1 <- MXDELT(y,1)
  #Coeficientes
  deltay1 <- DIFF(y)
  yusar1 <- deltay1[1:(N-1)]
  
  regre1 <- REGRESION(x1,yusar1)
  coef1 <- regre1$coef
  yestim1 <- regre1$ystim
  e1 <- regre1$e
  sigma1 <- regre1$sigma
  varcov1 <- regre1$varcov
  desves1 <- sqrt(varcov1[1,1])
  tsig1 <- coef1[1]/desves1
  vct1 <- qt((1-(alpha/2)),(N-(3)-1))
  ##Ruido Blanco
  lb1 <- LJUNGBOX(e1,0,1,0,alpha)
  lb11 <- lb1[1]
  vclb1 <- lb1[2]
  
  if(lb11<vclb1){
    pusar <- 1
    listap <- list(p=pusar, tabla=1, ystim=yestim1,e=e1,varcov=varcov1,desves=desves1, coef=coef1)
  } else{
    tablita <- matrix(NA,ncol = 3,nrow = p)
    row.names(tablita) <- c(1:p)
    colnames(tablita) <- c("Signif","RB","Choose")
    #SIGNIFICANCIA Y RUIDO BLANCO -----------------------------------------------------
    ##Significancia ----------------------------------
    for (i in 1:p) {
      
      #Matriz diseño esquema 3 
      x <- MXDELT(y,i)
      #Coeficientes
      deltay <- DIFF(y)
      yusar <- deltay[i:(N-1)]
      
      regre <- REGRESION(x,yusar)
      coef <- regre$coef
      yestim <- regre$ystim
      e <- regre$e
      sigma <- regre$sigma
      varcov <- regre$varcov
      desves <- sqrt(varcov[(i),(i)])
      tsig <- coef[(i)]/desves
      vct <- qt((1-(alpha/2)),(N-(i)-i))
      ##Ruido Blanco
      lb <- LJUNGBOX(e,0,i,0,alpha)
      lb1 <- lb[1]
      vclb <- lb[2]
      #Tabla de decisión ------------ 
      
      #Condicional para significancia
      if(tsig>vct|tsig<(-vct)){
        tablita[i,1] <- 1
      } else{
        tablita[i,1] <- 0
      }
      
      #Condicional para ruido blanco
      
      if(lb1<vclb){
        tablita[i,2] <- 1
      } else{
        tablita[i,2] <- 0
      }   
      
      #Condicional decisión
      if((tablita[i,2]+ tablita[i,1])==2){
        tablita[i,3] <- 1
      } else{
        tablita[i,3] <- 0
      }
      
      #Condicional para terminar iteración
      if(tablita[i,3]==1){
        break
      }
      
    } #Cierre primer for
    
    tablita <- na.omit(tablita)
    tablita
    
    if(p>1){
      pusar <- nrow(tablita)+1
    }
    
    listap <- list(p=pusar, tabla=tablita, ystim=yestim,e=e,varcov=varcov,desves=desves, coef=coef)
  }
  
  return(listap)
  
}


#prueba de raiz uniaria-----------------
RAIZESQUEMA1<- function(y,p){
  N <- length(y)
  #Para tau rechazo si el estadístico es menor al vc 
  #Para phis si phi es mayor al vc rechazo 
  tau <- TAUADF(y,p,1)
  t_a0 <- coef3[1]/sqrt((varcov3[1,1]))
  vca0 <- qt(1-(0.05/2),df = (N-p-1),lower.tail = F)
  
  if(tau<vctau){ #Rechazo al 5%
    resultado <- print(paste("Con un estadístico Tau de", tau,"Hay evidencia para rechazar la hipótesis nula y no hay raiz unitaria, el mejor modelo es un DGP1 con media nula y ", (p-1),"rezagos en parte aumentada"))
  } else{ #Hay raíz unitaria
    resultado <- print(paste("Existe evidencia para decir que hay raíz unitaria, modelar con un DGP3 con", p, "rezagos"))
  } 
  
  milista <- list(resultado= resultado, tau= tau, indiv= t_a0)
  return(milista)
}
