################################################
# Multiple linear regression
################################################

BLR_prop <- function(Y,A,X,
                     iters=20000,burn=1000,thin=1,
                     plot=FALSE,debug=FALSE,filename="bayesmod.txt"){

 library("R2OpenBUGS")
 tick <- proc.time()[3]

 sink(filename)
 cat("
   model{
     for(i in 1:N){
       Y[i]     ~ dnorm(muY[i],taue)
       muY[i]  <- int + beta*A[i] + X[i,1]*gamma[1] + X[i,2]*gamma[2] + X[i,3]*gamma[3] + 
                                    X[i,4]*gamma[4] + X[i,5]*gamma[5]
     }

     int   ~ dnorm(0,0.001)
     beta  ~ dnorm(0,0.001)
     taue  ~ dgamma(0.5,0.005)
     for(j in 1:5){
       gamma[j] ~ dnorm(0,0.001)
     }  
   }
 ", fill=TRUE)
 sink()

 N           <- length(Y)
 if(is.null(X)){X <- matrix(0,N,5)}  
 p         <- ncol(X)
 dat       <- list("N","Y","A","X")
 vars2keep <- list("beta","taue","gamma")
 inits     <- function(){list(sige=sd(Y),int=mean(Y),beta=0,gamma=rep(0,5))}

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




