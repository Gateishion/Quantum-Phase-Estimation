
arb_sampling <- function(f_d,sz,l_s, u_s){
#This function returns a sample of an
#arbitrary PDF
#f is a probability distribution
#size of sample

part <- 1000
x<-seq(l_s,u_s, by = abs(l_s-u_s)/part)

Smp <-  function(f){
  force(f)

  storage.vector <- rep(NA,length(x))
  for (i in 1:length(storage.vector)){
    storage.vector[i] <-  f(x[i])}
  return(storage.vector)
}

rsamp <- function(f,n){
  force(f)
  y1<-cumsum(Smp(f))*diff(seq(l_s,u_s, by = abs(l_s-u_s)/part))[1]
  pf <- approxfun(x,y1)
  qf <- approxfun(y1,x)
  return(qf(runif(n)))
}

return(rsamp(f_d,sz))

}

maxims <- function(a,tol,low,up) {
#This functions returns all local maximums of a
#function
#tol is the subdivitions

  force(a)
  x <- seq(low, up, length= tol)
  y <- vector("numeric", length = length(x))
  for (j in 1:length(x)) {
    y[j] <- a(x[j])
  }
  #print(x)
  m <- rep(NA, length(x)-1)
  for (i in 2:(length(x)-1)) {
    #print(c(y[i-1],y[i],y[i+1]))
    if( y[i-1] < y[i] & y[i] > y[i+1]){
    m[i] <- optimize(a, lower = x[i-1], x[i+1], maximum = T, tol=1e-16)$maximum }
  }

  #Delete de NAs
  m <- m[!is.na(m)]
  #m <- c(m,a(x[1]), a(x[length(x)]))
  m <- c(m,x[1], x[length(x)])
  #if (length(m)==0){m <-optimize(a, lower = low, up, maximum = T,tol=1e-16)$maximum }
  return(m)

}


choosemax <- function(atol,vec,f){
#choose the global maximum of a vector of maximums
#atol is the digits to compare
    force(f)
  #if( is.na(vec) == TRUE){
  #  return(NA)
  #}
  #v1 <- f(vec)
  v1 <- vector("numeric", length(vec))
  for (j in 1:length(vec)) {
    v1[j] <- round(f(vec[j]), digits = atol)
  }

  v2 <- vec[which(v1 == max(v1))]
  prop <- length(v2)
  #ss <- sample(1:prop, size = 1, prob = rep((1/(prop)), prop))
  ss <- sample(1:prop, size = 1)
  return(v2[ss])
}


logfunc2 <- function(f){
  #returns the log of a function
  force(f)
  function(...){log(f(...))}
}


sumfunct2<-function(a,b){
  #returns the sum of two functions
  force(a)
  force(b)
  function(...){a(...)+b(...)}
}


sumfctb2 <- function(vec,n){
  #sum vector of functions
  if (n==1) {vec[[1]]}
  else{
    for(j in 1:n){
      force(vec[[j]])
    }
    out <- vec[[1]]
    for(i in 2:n){
      out <- sumfunct2(out,vec[[i]])
    }
    return(out)
  }
}


pdfunct2<-function(a,b){
  #Product of 2 functs
  force(a)
  force(b)
  function(...){a(...)*b(...)}
}


pdfctb2 <- function(vec,n){
  #Product of list of functions
  if (n==1) {vec[[1]]}
  else{
    for(j in 1:n){
      force(vec[[j]])
    }
    out <- vec[[1]]
    for(i in 2:n){
      out <- pdfunct2(out,vec[[i]])
    }
    return(out)
  }
}


Holevo_Var <- function(est){
  #Returns the Holevo variance
  #est, vector of estimations
  v <- (CircStats::circ.disp(est)$rbar^(-2)-1)
  return(v)
}


#-----------------------------------
#Pauli Basis
#-----------------------------------
Pauli_M <- function(i){
#i in 0 to 3
S <- vector(mode="list", length = 4)
sigma0 <- matrix(nrow=2, ncol = 2, c(1,0,0,1))
sigma1 <- matrix(nrow=2, ncol = 2, c(0,1,1,0))
sigma2 <- matrix(nrow=2, ncol = 2, c(0,complex(real=0, imaginary = 1),complex(real = 0, imaginary = -1),0))
sigma3 <-  matrix(nrow=2, ncol = 2, c(1,0,0,-1))

S[[1]] <- sigma0
S[[2]] <- sigma1
S[[3]] <- sigma2
S[[4]] <- sigma3

if(i > 4 | i< 0){
  return( "i must be 0,1,2,3 or 4")
}

return(S[[i+1]])

}

P_vector <- function(v) {
  #returns the "dot product" of paulis with a vector
  sigma0 <- matrix(nrow=2, ncol = 2, c(1,0,0,1))
  sigma1 <- matrix(nrow=2, ncol = 2, c(0,1,1,0))
  sigma2 <- matrix(nrow=2, ncol = 2, c(0,complex(real=0, imaginary = 1),complex(real = 0, imaginary = -1),0))
  sigma3 <-  matrix(nrow=2, ncol = 2, c(1,0,0,-1))
  (v[1]*sigma1)+(v[2]*sigma2)+(v[3]*sigma3)
  }

Qubit <- function(a){
  #Returns a qubit
  #Input Bloch vector
  if(length(a)>4){
    return("the length of a must be 3")
  }
  asig <- P_vector(a)
  rho <- (1/2)*(Pauli_M(0)+asig)
  return(rho)
}


#--------------------------------------
#Covariant Estimation
#--------------------------------------

P_covariante <- function(qbit_vec,axis,phase,x){
#return p(x | theta) de la medida covariante M_star

n <- axis
J <- (1/2)*(P_vector(n))
rho <- Qubit(qbit_vec)
theta <- phase

SpecJ <- eigen(J)
Pj1 <- SpecJ$vectors[,1]%*%Conj(t(SpecJ$vectors[,1]))
Pjm1 <- SpecJ$vectors[,2]%*%Conj(t(SpecJ$vectors[,2]))
Pb1 <- Re(sum(mgcv::sdiag(Pj1%*%rho)))
Pbm1 <- Re(sum(mgcv::sdiag(Pjm1%*%rho)))
SuppJ <- c(SpecJ$values[1],SpecJ$values[2])
ProbJ <- c(Pb1, Pbm1)
MatofpbJ <- cbind(SuppJ, ProbJ)


eilam <- function(x,t,m) {complex(real = cos(m*(x-t)), imaginary = -sin(m*(x-t)))}
PSn <- MatofpbJ

create_eilam <- function(i){
  ei_out <- function(x,t) {eilam(x=x, t=t,m=i)*sqrt(PSn[,2][i])}
}
Sumandos <- lapply(1:length(PSn[,2]), create_eilam)

SP <- sumfctb2(Sumandos, length(Sumandos))
ProbSN <- function(x) {pracma::Real ((1/(2*pi))*(SP(x,t=theta)*Conj(SP(x,t=theta))))}
vfpbMS <- Vectorize(ProbSN)
return(vfpbMS(x))

}

#---------------------------------------
#Entangled
#---------------------------------------

Conv2 <- function(M,N) {
  kSamples::conv(M[,1],M[,2],N[,1],N[,2])
}
#Convolution of a vector of matrices of probabilities
ConvNb <- function(M,n){
  if (n==1) {
    M}
  else{
    out <- M
    for(i in 2:n){
      out <- Conv2(out,M)
    }
    return(out)
  }
}


Matrix_Prob_J <- function(axis,qbit_vec){

  n <- axis
  J <- (1/2)*(P_vector(n))
  rho <- Qubit(qbit_vec)

  SpecJ <- eigen(J)
  Pj1 <- SpecJ$vectors[,1]%*%Conj(t(SpecJ$vectors[,1]))
  Pjm1 <- SpecJ$vectors[,2]%*%Conj(t(SpecJ$vectors[,2]))
  Pb1 <- Re(sum(mgcv::sdiag(Pj1%*%rho)))
  Pbm1 <- Re(sum(mgcv::sdiag(Pjm1%*%rho)))
  SuppJ <- c(SpecJ$values[1],SpecJ$values[2])
  ProbJ <- c(Pb1, Pbm1)
  MatofpbJ <- cbind(SuppJ, ProbJ)

  return(MatofpbJ)

}




