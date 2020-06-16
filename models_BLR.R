################################################
# Multiple linear regression
################################################

BLR <- function(Y,A,
                iters=20000,burn=1000,thin=1,
                plot=FALSE,debug=FALSE,filename="bayesmod.txt"){

 library("R2OpenBUGS")
 tick <- proc.time()[3]

 sink(filename)
 cat("
   model{
     for(i in 1:N){
       Y[i]     ~ dnorm(muY[i],taue)
       muY[i]  <- int + beta*A[i] 
     }

     int   ~ dnorm(0,0.001)
     beta  ~ dnorm(0,0.001)
     taue  ~ dgamma(0.5,0.005)
   }
 ", fill=TRUE)
 sink()

 N         <- length(Y)
 dat       <- list("N","Y","A")
 vars2keep <- list("beta","taue")
 inits     <- function(){list(taue=1/var(Y),int=mean(Y),beta=0)}

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

