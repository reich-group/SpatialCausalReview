

################################################
# Joint spatial model for Y and A
################################################

joint_nofeedback <- function(Y,A,ADJ,
                             iters=20000,burn=1000,thin=1,
                             plot=FALSE,debug=FALSE,filename="bayesmod.txt"){

 library("R2OpenBUGS")
 tick <- proc.time()[3]

 sink(filename)
 cat("
   model{
     for(i in 1:N){
       Y[i]         ~ dnorm(muY[i],taue)
       A[i]         ~ dbin(q[i],1)
       muY[i]      <- int1 + beta*A[i] + U[i] + phi*V_cut[i]
       V_cut[i]    <- cut(V[i])
       logit(q[i]) <- int2 + V[i]
     }
     U[1:N] ~ car.normal(adj[], weights[], num[], tauu)  
     V[1:N] ~ car.normal(adj[], weights[], num[], tauv)  

     int1  ~ dnorm(0,0.1)
     int2  ~ dnorm(0,0.1)
     beta  ~ dnorm(0,0.1)
     phi   ~ dnorm(0,0.1)
     taue  ~ dgamma(0.5,0.005)
     tauu  ~ dgamma(0.5,0.005)
     tauv  ~ dgamma(0.5,0.005)
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
 vars2keep <- list("beta","phi","taue","tauu","tauv")
 inits     <- function(){list(taue=1/var(Y), tauu=0.01,int1=mean(Y),beta=0,phi=0,tauv=0.01,U=0*Y,V=0*Y)}

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
 out  <- list(beta=output$sims.list$beta,seconds=tock-tick)
return(out)}


