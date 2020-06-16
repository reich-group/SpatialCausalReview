
################################################
# Spatial regression with A as covariate
################################################

spatial_gaussian <- function(Y,A,ADJ,prop_score=NULL,
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
       base[i] <- int + U[i]
     }
     U[1:N] ~ car.normal(adj[], weights[], num[], tauu)  

     int   ~ dnorm(0,0.001)
     beta  ~ dnorm(0,0.001)
     taue  ~ dgamma(0.5,0.005)
     tauu  ~ dgamma(0.5,0.005)
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

 dat       <- list("N","Y","A","adj","num","weights")
 vars2keep <- list("beta","taue","tauu")
 if(!is.null(prop_score)){
   vars2keep <- list("beta","taue","tauu","base")
 }
 inits     <- function(){list(taue=1/var(Y), tauu=10,int=mean(Y),beta=0,U=0*Y)}

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
  # Post-hoc IDW debiasing
  if(!is.null(prop_score)){
     base  <- output$sims.list$base
     sigma <- 1/sqrt(output$sims.list$taue)
    
     b     <- NULL
     for(i in 1:length(beta)){
       Y0   <- base[i,]
       Y1   <- base[i,] + beta[i]
       temp <- (A*Y - (A-prop_score)*Y1)/prop_score - 
               ((1-A)*Y - (prop_score-A)*Y0)/(1-prop_score)
       b    <- c(b,mean(temp))
     }     
     beta <- b
  }


 if(plot){plot(output)}
 tock <- proc.time()[3]

 out  <- list(beta=beta,seconds=tock-tick)
return(out)}


