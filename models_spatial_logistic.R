################################################
# Spatial logistic regression
################################################

spatial_logit <- function(A,ADJ,
                          iters=5000,burn=1000,thin=1,
                          plot=FALSE,debug=FALSE,keep_samples=FALSE,
                          filename="bayesmod.txt"){

 library("R2OpenBUGS")
 tick <- proc.time()[3]

 sink(filename)
 cat("
   model{
     for(i in 1:N){
       A[i]         ~ dbin(q[i],1)
       logit(q[i]) <- int + V[i]
     }
     V[1:N] ~ car.normal(adj[], weights[], num[], tauv)  

     int   ~ dnorm(0,0.1)
     tauv  ~ dgamma(0.05,0.005)  
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

 dat       <- list("N","A","adj","num","weights")
 vars2keep <- list("q")
 inits     <- function(){list(int=0,tauv=10,V=0*A)}

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

 if(plot){plot(output)}
 tock <- proc.time()[3]
 q    <- output$sims.list$q 
 if(!keep_samples){q<-colMeans(q)}
 out  <- list(prop_score=q,seconds=tock-tick)
return(out)}

