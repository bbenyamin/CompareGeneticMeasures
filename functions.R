#Nolan App Functions#
##########################################
# Function 1: Calculates measures of gene effects using lambdas
# Function 2: Calculates measures of gene effects using h2L
# Function 3: Reads parameters.txt for explanation in R shiny

#####################################################################################################################
####################################################################################################################
#Function 1: Calculation with lambdas
# Calculates measures of gene effects:
# liability threshold variance, risk scale lambdas, AUC and PAR for single biallelic loci
# p: allele frequency of the risk allele (RAF)
# RRBb: relative risk of one risk allele compared to no risk allele
# RRBB: relative risk of two risk alleles compared to no risk allele = RRBb^2
# K: baseline population risk of disease (assumed)
# lambdas: sibling relative risk (assumed) or h2L (heritability of liability)

get_Vgs_lambdas <- function(variant,p,RRBb,RRBB,K,lambdas){
  
    #########################
  # Falconer: heritability of liability given K and lambdas.
  T0=qnorm(1-K)       # Threshold among population
  z<-dnorm(T0)  			# Height of the normal distribution at threshold     
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
  result_simple=list(h2Li=h2Li*100,h2Liprop=h2Liprop*100,h2Lipropapprox=h2Lapprox*100,lambdasib=lambdasib,proploglambdasib=proplambdas*100,proplambS2=proplambS2*100,AUCi=AUCSNP,AUCprop=AUCprop*100,PAR=PAR*100)
  
 }

#################################################################################
#################################################################################
# Function 2: Calculates measures of gene effects using h2L (only 2 lines are different)

get_Vgs_h2L <- function(variant,p,RRBb,RRBB,K,h2L){
  
  #########################
  # Falconer: heritability of liability given K and lambdas.
  T0=qnorm(1-K)       # Threshold among population
  z<-dnorm(T0)    		# Height of the normal distribution at threshold     
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
  result_simple=list(h2Li=h2Li*100,h2Liprop=h2Liprop*100,h2Lipropapprox=h2Lapprox*100,lambdasib=lambdasib,proploglambdasib=proplambdas*100,proplambS2=proplambS2*100,AUCi=AUCSNP,AUCprop=AUCprop*100,PAR=PAR*100)
  
}

#################################################################################
################################################################################
# Function 3: Reads parameters.txt for explanation in R shiny
## Explanations of variables###
parameters <- read.table("parameters.txt", sep = ",", header = T)
inputs <- read.table("inputs.txt", sep=",", header=F)



