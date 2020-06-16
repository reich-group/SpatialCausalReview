



################################################
# Spatial regression with A and X as covariates
################################################

spatial_prop <- function(Y,A,ADJ,prop_score=NULL,X=NULL,
                        iters=20000,burn=1000,thin=1,
                        plot=FALSE,debug=FALSE,
                        filename="bayesmod.txt"){

 library("R2OpenBUGS")
 tick <- proc.time()[3]

 sink(filename)
 cat("
   model{
     for(i in 1:N){
       Y[i]     ~ dnorm(muY[i],taue)
       muY[i]  <- beta*A[i] + base[i]
       base[i] <- int + X[i,1]*gamma[1] + X[i,2]*gamma[2] + X[i,3]*gamma[3] + 
                        X[i,4]*gamma[4] + X[i,5]*gamma[5] + U[i]
     }
     U[1:N] ~ car.normal(adj[], weights[], num[], tauu)  

     int   ~ dnorm(0,0.001)
     beta  ~ dnorm(0,0.001)
     taue  ~ dgamma(0.5,0.005)
     tauu  ~ dgamma(0.5,0.005)
     for(j in 1:5){
       gamma[j] ~ dnorm(0,0.001)
     }  
   }
 ", fill=TRUE)
 sink()

 N           <- ncol(ADJ)
 num         <- apply(ADJ,2,sum)
 sumNumNeigh <- sum(num)
 adj         <- NULL
 for(j1 in 1:N){for(j2 in 1:N){
   if(ADJ[j1,j2]){adj<-c(adj,j2)}
 }}
 weights <- rep(1,sumNumNeigh)
 if(is.null(X)){X <- matrix(0,N,5)}  

 dat       <- list("N","Y","A","X","adj","num","weights")
 vars2keep <- list("beta","taue","tauu","gamma")
 if(!is.null(prop_score)){
   vars2keep <- list("beta","taue","tauu","gamma","base")
 }
 inits     <- function(){list(taue=1/var(Y), tauu=10,int=mean(Y),beta=0,U=0*Y,gamma=rep(0,5))}

 output<-bugs(
   model.file=filename,
   data=dat,
   inits = inits,
   parameters.to.save = vars2keep,
   n.chains=2,
   n.iter=iters,
   n.burnin=burn,
   n.thin=thin,
   debug=debug,
   DIC=FALSE
 )

 beta  <- output$sims.list$beta

 if(plot){plot(output)}
 tock <- proc.time()[3]

 out  <- list(beta=beta,seconds=tock-tick)
return(out)}



