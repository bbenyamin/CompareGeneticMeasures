#Nolan App Plot Functions#
##########################################
# Function 1: Produces graph. First calculates measures of gene effects using lambdas
# Function 2: CProduces graph. First calculates measures of gene effects using h2L

####################################################################################
####################################################################################
# Function 1:
get_plot_lambdas <- function(variant,p,RRBb,RRBB,K,lambdas){
  
  #########################
  # Falconer: heritability of liability given K and lambdas.
  T0=qnorm(1-K)       # Threshold among population
  z<-dnorm(T0)      	# Height of the normal distribution at threshold     
  i<-z/K        			# Mean liability of affected group (selection intensity)
  T1 <- qnorm(1 - lambdas * K)
  h2L <- 2*(T0-T1*sqrt(1-(T0^2-T1^2)*(1-T0/i)))/(i+T1^2*(i-T0))
  
  kbb <- K/((1-p)^2+2*p*(1-p)*RRBb+p*p*RRBB)
  PARH=(K-kbb)/K
  
  # Check for impossible values (pDGeno3>1); then constrain 
  n = length(RRBb)
  flag = c(rep(0,n))
  flag[kbb*RRBB >1]=1
  pGeno1 <- (1-p)^2
  pGeno2 <- 2*p*(1-p)
  pGeno3 <- p*2
  kbb[kbb*RRBB>1]=(K-pGeno3)/(pGeno1+pGeno2*RRBb)
  
  ##########################
  #functions
  #follows Falconer & Mackay Table 7.3
  get_a=function(wbb,wBb,wBB){(wBB-(wbb+wBB)/2)}
  get_d=function(wbb,wBb,wBB){wBb-(wbb+wBB)/2}
  get_VA=function(p,a,d){2*p*(1-p)*(a+d*((1-p)-p))^2}
  get_VD=function(p,a,d){(2*p*(1-p)*d)^2}
  
  ######################
  # Observed Risk scale
  # 1) Risch 1990 Nature; Hemminiki 2008 (Typo in Risch eqn; fixed below)
  wbbO=kbb
  wBbO=kbb*RRBb
  wBBO=kbb*RRBB
  aO=get_a(wbbO,wBbO,wBBO)
  dO=get_d(wbbO,wBbO,wBBO)
  VAO=get_VA(p,aO,dO)
  VDO=get_VD(p,aO,dO)
  VGO = VAO+VDO
  # check
  Vao = kbb^2*2*p*(1-p)*(p*(RRBB-RRBb)+(1-p)*(RRBb-1))^2 #Additive genetic variance, observed scale =VAO 
  Vdo = kbb^2*p^2*(1-p)^2*(RRBB+1-2*RRBb)^2			 # Dominance, observed scale =VDO 
  VGOsib=(0.5*VAO+0.25*VDO)
  Vgsib = 1 + (0.5*VAO+0.25*VDO)/(K^2) #Equivalent to lambdasib
  PARH <- (K-kbb)/K 
  lambdasib <- 1+(0.5*VAO+0.25*VDO)/(kbb^2)*((1-PARH)^2);
  proplambdas <- log(lambdasib)/log(lambdas) # Proportion of sibling risk explained, on log scale
  
  #########################
  # Falconer: liability threshold h^2_SNP
  # Falconer & Mackay; Sham Human Genetics text book
  wbbL <- -qnorm(1-wbbO)  
  wBbL <- -qnorm(1-wBbO)
  wBBL <- -qnorm(1-wBBO)
  aL <- get_a(wbbL,wBbL,wBBL)
  dL <- get_d(wbbL,wBbL,wBBL)
  VAL <- get_VA(p,aL,dL)
  VDL <- get_VD(p,aL,dL)
  VGL <- VAL+VDL
  H2Li<- VGL/(1+VGL)
  h2Li<- VAL/(1+VGL)
  h2Liprop <- h2Li/h2L
  #approximation from Falconer & Mackay as used in Purcell et al ISC, Nature 2009
  VGLapprox <- 2*p*(1-p)*(RRBb-1)^2/(i^2) 
  h2Lapprox <- VGLapprox/(1+VGLapprox)
  
  ##########################
  # Log Risk
  wbbl<-log(1)
  wBbl<-log(RRBb)
  wBBl<-log(RRBB)
  al <- get_a(wbbl,wBbl,wBBl)
  dl <- get_d(wbbl,wBbl,wBBl)
  VAl <- get_VA(p,al,dl)
  VDl <- get_VD(p,al,dl)
  VGl <- VAl+VDl
  #check  Pharoah 2008 NEJM and Pharaoh et al, 2002
  #m=2*p*(1-p)*log(RRBb)+ p*p*log(RRBB)
  #VGl=(1-p)*(1-p)*m*m+2*p*(1-p)*(log(RRBb)-m)^2+p*p*(log(RRBB)-m)^2
  vgEaston <- exp(VGl)
  lambda2S <- 2*lambdas-1
  lambdaS2 <- lambdas^2  #MZ lambda
  proplamb2S <-VGl/log(lambda2S)
  proplambS2 <- VGl/log(lambdaS2)
  
  #########################  
  # Population attributable risk (PAR) Witte
  PARnum<-2*p*(1-p)*(RRBb-1)+p*p*(RRBB-1)
  PAR <- PARnum/(1+PARnum)  
  
  #########################
  # AUC Calculations; get AUCmax for h2L and AUCSNP for Vgliab
  v  <- -i * (K / (1-K))
  dn <- function(rho,i,v,T0){pnorm((i-v)*rho/sqrt(rho*(1-rho*i*(i-T0)+1-rho*v*(v-T0))))}
  AUCSNP <- dn(h2Li,i,v,T0)
  AUCMAX <- dn(h2L,i,v,T0)
  AUCprop <- ((AUCSNP-0.5)/(AUCMAX-0.5))^2
  
  ##############
  #Output for plot
  d <- data.frame(raf=p, 
                       or=RRBb, 
                       h2=h2Liprop*100,
                       H2approx=h2Lapprox*100,
                       L_S=proplambdas*100,
                       L_S2=proplambS2*100,
                       AUC_Prop=AUCprop*100) 
  
  #Defining variables to use in plot:
  raf<-p
  or<-RRBb
  h2<-h2Liprop*100
  H2approx<-h2Lapprox*100
  L_S<-proplambdas*100
  L_S2<-proplambS2*100
  AUC_Prop<-AUCprop*100
  
  #Defining variables for x axis labels:
  sumh2=sum(d$h2)
  sumh2=format(sumh2,digits=3)
  x1=paste("(",sumh2,"%)", sep="", collapse="")
  x1=paste("Heritability", x1, sep="\n")
  
  sumH2approx=sum(d$H2approx)
  sumH2approx=format(sumH2approx,digits=3)
  x2=paste("(",sumH2approx,"%)", sep="", collapse="")
  x2=paste("Approx. Herit.", x2, sep="\n")
  
  sumL_S=sum(d$L_S)
  sumL_S=format(sumL_S,digits=3)
  x3=paste("(",sumL_S,"%)", sep="", collapse="")
  x3=paste("Sibling RR", x3, sep="\n")
  
  sumL_S2=sum(d$L_S2)
  sumL_S2=format(sumL_S2,digits=3)
  x4=paste("(",sumL_S2,"%)", sep="", collapse="")
  x4=paste("Family RR", x4, sep="\n")
  
  sumAUC_Prop=sum(d$AUC_Prop)
  sumAUC_Prop=format(sumAUC_Prop,digits=3)
  x5=paste("(",sumAUC_Prop,"%)", sep="", collapse="")
  x5=paste("AUC", x5, sep="\n")
  
##NOW PLOTTING:

mmax = max(L_S) + 2

mtitle = "" #Add plot title if desired
ch_l=c(x1, x2, x3, x4, x5)

ch = c('h2', 'H2approx', 'L_S', 'L_S2', 'AUC_Prop')
ch_m = c(1,1,1,1,1)

tr = function(t)
  return(sqrt(t))
itr = function(t)
  return(t^2)

## Start up the plotting

FUDGE = 0.15

#mmax = max(d[, ch], na.rm=TRUE)
plot(c(1+FUDGE, length(ch)), tr(c(0, mmax)), type='n', axes=FALSE, ylab='Percentage', xlab='Measure', font.lab=2, main=mtitle)
alpha=0.5

for(i in 1:nrow(d)){ #d is reading in the data from infile (rows=different SNPs). This could be done still as for loop, calling data from previous function (which would be included above)
  or = d$or[i] #or=RR
  raf = d$raf[i]
  
  col = 'yellow' ## safety precaution, should never show
  # Common, low penetrance
  if(or<=1.3) col=rgb(0,1,0,alpha) #'green'
  if(or>1.3) col=rgb(0,0,1,alpha) #'blue'
  if(or>2) col=rgb(1,0,0,alpha) #red
  if(or>15) col=rgb(0,0,0,alpha) #black
  
  lines(1:5, tr(d[i, ch]), col=col)
}

axis(2, at=axTicks(2), labels=itr(axTicks(2)), lwd=0, lwd.ticks=1)
axis(1, at=1:length(ch), ch_l, lwd.ticks=0)
axis(3, labels=NA, lwd.ticks=0)
for(x in 1:length(ch))
  abline(v=x, col='black')
  
  legend("top", legend=c("RR<1.3", "1.3<RR<2", "2<RR<15", "RR>15"),col=c("green","blue","red","black"), horiz=T,lwd=1.5, xpd=T)
  
}

####################################################################################
####################################################################################
# Function 2: (only 2 lines are different to function 1)
get_plot_h2L <- function(variant,p,RRBb,RRBB,K,h2L){
  
  #########################
  # Falconer: heritability of liability given K and lambdas.
  T0=qnorm(1-K)       # Threshold among population
  z<-dnorm(T0)        # Height of the normal distribution at threshold     
  i<-z/K        			# Mean liability of affected group (selection intensity)
  T1 <- (2 * T0 - i * h2L) / sqrt(4 - h2L^2 * i * (i - T0)) #DIFFERENT to function 1
  lambdas <- (1 - pnorm(T1)) / K #DIFFERENT to function 1
  
  kbb <- K/((1-p)^2+2*p*(1-p)*RRBb+p*p*RRBB)
  PARH=(K-kbb)/K
  
  # Check for impossible values (pDGeno3>1); then constrain 
  n = length(RRBb)
  flag = c(rep(0,n))
  flag[kbb*RRBB >1]=1
  pGeno1 <- (1-p)^2
  pGeno2 <- 2*p*(1-p)
  pGeno3 <- p*2
  kbb[kbb*RRBB>1]=(K-pGeno3)/(pGeno1+pGeno2*RRBb)
  
  ##########################
  #functions
  #follows Falconer & Mackay Table 7.3
  get_a=function(wbb,wBb,wBB){(wBB-(wbb+wBB)/2)}
  get_d=function(wbb,wBb,wBB){wBb-(wbb+wBB)/2}
  get_VA=function(p,a,d){2*p*(1-p)*(a+d*((1-p)-p))^2}
  get_VD=function(p,a,d){(2*p*(1-p)*d)^2}
  
  ######################
  # Observed Risk scale
  # 1) Risch 1990 Nature; Hemminiki 2008 (Typo in Risch eqn; fixed below)
  wbbO=kbb
  wBbO=kbb*RRBb
  wBBO=kbb*RRBB
  aO=get_a(wbbO,wBbO,wBBO)
  dO=get_d(wbbO,wBbO,wBBO)
  VAO=get_VA(p,aO,dO)
  VDO=get_VD(p,aO,dO)
  VGO = VAO+VDO
  # check
  Vao = kbb^2*2*p*(1-p)*(p*(RRBB-RRBb)+(1-p)*(RRBb-1))^2 #Additive genetic variance, observed scale =VAO 
  Vdo = kbb^2*p^2*(1-p)^2*(RRBB+1-2*RRBb)^2			 # Dominance, observed scale =VDO 
  VGOsib=(0.5*VAO+0.25*VDO)
  Vgsib = 1 + (0.5*VAO+0.25*VDO)/(K^2) #Equivalent to lambdasib
  PARH <- (K-kbb)/K 
  lambdasib <- 1+(0.5*VAO+0.25*VDO)/(kbb^2)*((1-PARH)^2);
  proplambdas <- log(lambdasib)/log(lambdas) # Proportion of sibling risk explained, on log scale
  
  #########################
  # Falconer: liability threshold h^2_SNP
  # Falconer & Mackay; Sham Human Genetics text book
  wbbL <- -qnorm(1-wbbO)  
  wBbL <- -qnorm(1-wBbO)
  wBBL <- -qnorm(1-wBBO)
  aL <- get_a(wbbL,wBbL,wBBL)
  dL <- get_d(wbbL,wBbL,wBBL)
  VAL <- get_VA(p,aL,dL)
  VDL <- get_VD(p,aL,dL)
  VGL <- VAL+VDL
  H2Li<- VGL/(1+VGL)
  h2Li<- VAL/(1+VGL)
  h2Liprop <- h2Li/h2L
  #approximation from Falconer & Mackay as used in Purcell et al ISC, Nature 2009
  VGLapprox <- 2*p*(1-p)*(RRBb-1)^2/(i^2) 
  h2Lapprox <- VGLapprox/(1+VGLapprox)
  
  ##########################
  # Log Risk
  wbbl<-log(1)
  wBbl<-log(RRBb)
  wBBl<-log(RRBB)
  al <- get_a(wbbl,wBbl,wBBl)
  dl <- get_d(wbbl,wBbl,wBBl)
  VAl <- get_VA(p,al,dl)
  VDl <- get_VD(p,al,dl)
  VGl <- VAl+VDl
  #check  Pharoah 2008 NEJM and Pharaoh et al, 2002
  #m=2*p*(1-p)*log(RRBb)+ p*p*log(RRBB)
  #VGl=(1-p)*(1-p)*m*m+2*p*(1-p)*(log(RRBb)-m)^2+p*p*(log(RRBB)-m)^2
  vgEaston <- exp(VGl)
  lambda2S <- 2*lambdas-1
  lambdaS2 <- lambdas^2  #MZ lambda
  proplamb2S <-VGl/log(lambda2S)
  proplambS2 <- VGl/log(lambdaS2)
  
  #########################  
  # Population attributable risk (PAR) Witte
  PARnum<-2*p*(1-p)*(RRBb-1)+p*p*(RRBB-1)
  PAR <- PARnum/(1+PARnum)  
  
  #########################
  # AUC Calculations; get AUCmax for h2L and AUCSNP for Vgliab
  v  <- -i * (K / (1-K))
  dn <- function(rho,i,v,T0){pnorm((i-v)*rho/sqrt(rho*(1-rho*i*(i-T0)+1-rho*v*(v-T0))))}
  AUCSNP <- dn(h2Li,i,v,T0)
  AUCMAX <- dn(h2L,i,v,T0)
  AUCprop <- ((AUCSNP-0.5)/(AUCMAX-0.5))^2
  
  ##############
  #Output for plot
  d <- data.frame(raf=p, 
                  or=RRBb, 
                  h2=h2Liprop*100,
                  H2approx=h2Lapprox*100,
                  L_S=proplambdas*100,
                  L_S2=proplambS2*100,
                  AUC_Prop=AUCprop*100) 
  
  #Defining variables to use in plot:
  raf<-p
  or<-RRBb
  h2<-h2Liprop*100
  H2approx<-h2Lapprox*100
  L_S<-proplambdas*100
  L_S2<-proplambS2*100
  AUC_Prop<-AUCprop*100
  
  #Defining variables for x axis labels:
  sumh2=sum(d$h2)
  sumh2=format(sumh2,digits=3)
  x1=paste("(",sumh2,"%)", sep="", collapse="")
  x1=paste("Heritability", x1, sep="\n")
  
  sumH2approx=sum(d$H2approx)
  sumH2approx=format(sumH2approx,digits=3)
  x2=paste("(",sumH2approx,"%)", sep="", collapse="")
  x2=paste("Approx. Herit.", x2, sep="\n")
  
  sumL_S=sum(d$L_S)
  sumL_S=format(sumL_S,digits=3)
  x3=paste("(",sumL_S,"%)", sep="", collapse="")
  x3=paste("Sibling RR", x3, sep="\n")
  
  sumL_S2=sum(d$L_S2)
  sumL_S2=format(sumL_S2,digits=3)
  x4=paste("(",sumL_S2,"%)", sep="", collapse="")
  x4=paste("Family RR", x4, sep="\n")
  
  sumAUC_Prop=sum(d$AUC_Prop)
  sumAUC_Prop=format(sumAUC_Prop,digits=3)
  x5=paste("(",sumAUC_Prop,"%)", sep="", collapse="")
  x5=paste("AUC", x5, sep="\n")
  
  ##NOW PLOTTING:
  
  mmax = max(L_S) + 2
  
  mtitle = "" #Add plot title if desired
  ch_l=c(x1, x2, x3, x4, x5)
  
  ch = c('h2', 'H2approx', 'L_S', 'L_S2', 'AUC_Prop')
  ch_m = c(1,1,1,1,1)
  
  tr = function(t)
    return(sqrt(t))
  itr = function(t)
    return(t^2)
  
  ## Start up the plotting
  
  FUDGE = 0.15
  
  #mmax = max(d[, ch], na.rm=TRUE)
  plot(c(1+FUDGE, length(ch)), tr(c(0, mmax)), type='n', axes=FALSE, ylab='Percentage', xlab='Measure', font.lab=2, main=mtitle)
  alpha=0.5
  
  for(i in 1:nrow(d)){ #d is reading in the data from infile (rows=different SNPs). This could be done still as for loop, calling data from previous function (which would be included above)
    or = d$or[i] #or=RR
    raf = d$raf[i]
    
    col = 'yellow' ## safety precaution, should never show
    # Common, low penetrance
    if(or<=1.3) col=rgb(0,1,0,alpha) #'green'
    if(or>1.3) col=rgb(0,0,1,alpha) #'blue'
    if(or>2) col=rgb(1,0,0,alpha) #red
    if(or>15) col=rgb(0,0,0,alpha) #black
    
    lines(1:5, tr(d[i, ch]), col=col)
  }
  
  axis(2, at=axTicks(2), labels=itr(axTicks(2)), lwd=0, lwd.ticks=1)
  axis(1, at=1:length(ch), ch_l, lwd.ticks=0)
  axis(3, labels=NA, lwd.ticks=0)
  for(x in 1:length(ch))
    abline(v=x, col='black')
  
  legend("top", legend=c("RR<1.3", "1.3<RR<2", "2<RR<15", "RR>15"),col=c("green","blue","red","black"), horiz=T,lwd=1.5, xpd=T)
  
}

