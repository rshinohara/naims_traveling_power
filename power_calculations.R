################################################
# Power calculations for NAIMS Traveling Study Proposal
# Russell Shinohara
# February 23, 2018
################################################

# This function gets the ANOVA F statistic for site.
get.F<-function(y,site) {
  suppressWarnings(full.model<-lm(y~site))
  suppressWarnings(af <- anova(full.model))
  af$"F value"[1]
}

# This function does a permutation test for the null of no site differences (by permuting sites within subjects).
perm.test.F<-function(y,site,id,B.perm=100) {
  y.demean<-resid(lm(y~id))
  obsd.F<-get.F(y.demean,site)
  perm.F<-rep(NA,B.perm)
  for (bb in 1:B.perm) {
    perm.F[bb]<-get.F(y.demean[sample.int(length(y))],site)[1]
  }
  mean(obsd.F[1]<perm.F)
}

# Parametric bootstrap for random slope https://stats.stackexchange.com/questions/4858/how-to-test-random-effects-in-a-multilevel-model-in-r

bootstrapAnova <- function(mA, m0, B=1000){
  oneBootstrap <- function(m0, mA){
    d <- drop(simulate(m0))
    try(m2 <-refit(mA, newresp=d))
    try(m1 <-refit(m0, newresp=d))
    xx<-try(out<-anova(m2,m1)$Chisq[2])
    if(!inherits(xx, "try-error")) {
      return(xx) 
    } else {
      return(NA)
    }
  }  
  nulldist <-  if(!require(multicore)){
    replicate(B, oneBootstrap(m0, mA))
  } else {
    unlist(mclapply(1:B, function(x) oneBootstrap(m0, mA)))
  }   
  ret <- anova(mA, m0)
  ret$"Pr(>Chisq)"[2] <- mean(ret$Chisq[2] < nulldist,na.rm=TRUE)
  names(ret)[7] <- "Pr_boot(>Chisq)"
  attr(ret, "heading") <- c(attr(ret, "heading")[1], 
                            paste("Parametric bootstrap with", B,"samples."),
                            attr(ret, "heading")[-1])
  attr(ret, "nulldist") <- nulldist
  return(ret)
}

#This function interleaves two vectors
merge.vecs<-function(x,y) {
  if (length(x)!=length(y)) { 
    print('error - vectors of differing lengths') 
    return(NULL)
  } else {
    z<-rep(NA,2*length(x))
    z[2*(1:(length(x)))-1]<-x
    z[2*(1:(length(x)))]<-y
    return(z)
  }
}

################################################


set.seed(529349245)
library(nlme)
library(pbkrtest)
library(dplyr)
library(parallel)

################################################



B<-1000 # Number of simulations
B.perm<-1000 # Number of permutations to use for each permutation test
J<-8 # number of sites
n.1<-3 # number of subjects in Phase I
n.2<-32 # number of subjects in Phase II

m<-3 # number of subjects recruited at each site

t.vol<-13.5 # mean volume of thalamus
l.vol<-18 # mean volume of T2 lesions
sd.t.vol<-2 # sd of volume of thalamus
sd.l.vol<-2 # sd of volume of T2 lesions

min.dd<-1 # min disease duration
max.dd<-20  # max disease duration

min.age<-18 # min disease duration
max.age<-70 # max disease duration

p.female<-2/3 #Based on 2:1 ratio reported at http://www.sciencemag.org/site/feature/data/983519.xhtml

t.age.beta<-(12-14)/40 #based on figure 2 of http://www.sciencedirect.com/science/article/pii/S0926641001000106
t.sex.beta<--1 # based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2262922/
t.dd.beta<--0.17*t.vol/10 # based on Thalamic atrophy and cognition in multiple sclerosis, Neurology, 2007

l.dd.beta<-4/10 #based on table 1 of http://jnnp.bmj.com/content/jnnp/74/11/1551.full.pdf
l.sex.beta<-0 # based on https://www.ncbi.nlm.nih.gov/pubmed/16043800

sd.t.meas<-0.4 ## estimated from pilot data FIRST thalamic volumes - using sample sd as sqrt of average Euclidean distance /2 between observations within site (based on U-statistic sample variance representation)
sd.l.meas<-0.5 ## estimated from pilot data manual T2 lesion volumes - using sample sd as sqrt of average Euclidean distance /2 between observations within site (based on U-statistic sample variance representation)

sd.t.site<-0.4 ## estimated from pilot data FIRST thalamic volumes - using sample sd of first acquisitions on Skyra scanners. 
sd.l.site<-0.6 # estimated from pilot data manual T2 lesion volumes - using sample sd of first acquisitions on Skyra scanners. 

scanner.models<-c(1,2,0,2,1,1,2,1) # based on scanner model (2=Prisma, 1=Skyra, 0=Trio)
coil.channels<-factor(c(20,20,64,64,32,32,32,32,32,32,20,20,64,64,20,20,32,32)) 
## NOTE - there is an extra 2 entries in coil.channels since we're scanning twice on the last (NIH) scanner
## NOTE - scanner order is BWH, CHOP, JHU, MGH, Sinai, Cornell, Penn, NIH

################################################
################################################
############## CROSS-SECTIONAL ANALYSES
################################################
################################################

### NOTE - for simplicity, these power calculations do not assume disease-duration random slopes. Although this indicates that effects of scanner differences in Aim 2 will not be present, this is less relevant for Aim 1 power calculations and allows us to use pilot data directly.

### CROSS-SECTIONAL ANALYSIS - PHASE I POWER CALCULATIONS FOR DETECTING SITE EFFECTS
get.power.phase1<-function(bias.par,B=500) {
  id<-factor(rep(101:103,2*(J+1)))  
  site<-factor(rep(c(1:J,J),each=2*n.1))
  l.p<-rep(NA,B);t.p<-rep(NA,B) # tests for proportion of variation explained by site
  
  for (b in 1:B) {
    #Sample demographics
    age<-runif(3,min.age,max.age)
  	onset<-pmax(age-runif(3,min.dd,max.dd),10)
  	dd<-age-onset
  	sex<-c(1,1,0) # ASSUME WE WILL RECRUIT 2 FEMALES FOR PHASE I

    #Sample imaging measures
  	t<-rnorm(n.1,t.vol,sd.t.vol) + (age-50)*t.age.beta + sex*t.sex.beta + dd*t.dd.beta
  	l<-rnorm(n.1,l.vol,sd.l.vol) + sex*l.sex.beta + dd*l.dd.beta

  	#Sample site effects - add NIH site again for second scan session
  	site.eff.t<-rnorm(J,0,sd.t.site);site.eff.t<-c(site.eff.t,site.eff.t[length(site.eff.t)])
    site.eff.l<-rnorm(J,0,sd.l.site);site.eff.l<-c(site.eff.l,site.eff.l[length(site.eff.l)])

    #Specify acquisition biases
  	scanner.bias<-bias.par*c(scanner.models,1) # based on scanner model (2=Prisma, 1=Skyra, 0=Trio)
    coil.bias<-bias.par/2*as.numeric(coil.channels) # conservatively, assume receive coils are responsible for an additional 50% of variation

    #Sample observed data with both scanner and coil effects
  	t.obs<-rep(t,2*(J+1))+rep(scanner.bias,each=2*n.1)+rep(coil.bias,each=n.1)+rep(site.eff.t,each=2*n.1)+rnorm(2*(J+1)*n.1,0,sd.t.meas)
  	l.obs<-rep(l,2*(J+1))+rep(scanner.bias,each=2*n.1)+rep(coil.bias,each=n.1)+rep(site.eff.l,each=2*n.1)+rnorm(2*(J+1)*n.1,0,sd.l.meas)

    #Do permutation test for site effects (any)
    t.p[b]<-perm.test.F(t.obs,site,id)
    l.p[b]<-perm.test.F(l.obs,site,id)
  }

  c(mean(t.p<0.05),mean(l.p<0.05))
}

set.seed(529349245)
mclapply(seq(0,1,by=0.1),get.power.phase1,mc.cores=20)

### CROSS-SECTIONAL ANALYSIS - PHASE II POWER CALCULATIONS FOR DETECTING SCANNER AND COIL EFFECTS
get.power.phase2<-function(bias.par,B=500) {  

	l.p.scanner<-rep(NA,B);t.p.scanner<-rep(NA,B) # tests for variation explained by scanner
  l.p.coil<-rep(NA,B);t.p.coil<-rep(NA,B) # tests for variation explained by coil
	pb <- txtProgressBar(min = 0, max = B, style = 3)

	for (b in 1:B) {
  
		#Sample patient-level variables
		id.phase2<-1:n.2
		age<-runif(n.2,min.age,max.age)
		onset<-pmax(age-runif(n.2,min.dd,max.dd),10)
		dd<-age-onset
		sex<-rbinom(n.2,1,p.female)
		t<-rnorm(n.2,t.vol,sd.t.vol) + (age-50)*t.age.beta + sex*t.sex.beta + dd*t.dd.beta
		l<-rnorm(n.2,l.vol,sd.l.vol) + sex*l.sex.beta + dd*l.dd.beta
	
    #Sample site variable
    site<-NULL
    home.sites<-rep(1:J,each=n.2/J)
    for (j in 1:n.2) {
      site<-c(site,home.sites[j],sample(1:J,2))
    }

    #Sample imaging biases
    site.eff.t<-rnorm(J,0,sd.t.site)
    site.eff.l<-rnorm(J,0,sd.l.site)

    t.obs.phase2<-NULL
    l.obs.phase2<-NULL
    age.vec.phase2<-NULL
    sex.vec.phase2<-NULL
    dd.vec.phase2<-NULL
    scanner.bias.vec.phase2<-NULL
    id.vec.phase2<-NULL
    site.fac.phase2<-NULL
    coil.fac.phase2<-NULL

    #Sample imaging variables
    for (j in 1:n.2) { # for each subject
      for (k in 1:m) { # for each site that subject is imaged at

        this.visit.site<-site[rep(id.phase2,each=3)==j][k]
        # if NIH, who will scan subjects on two coils
        if (this.visit.site==J) {
          t.obs.phase2<-c(t.obs.phase2,rep(t[j],4)+bias.par*rep(scanner.models[J],4)+bias.par/2*as.numeric(coil.channels)[(2*J-1):length(coil.channels)]+rep(site.eff.t[J],4)+rnorm(4,0,sd.t.meas))
          l.obs.phase2<-c(l.obs.phase2,rep(l[j],4)+bias.par*rep(scanner.models[J],4)+bias.par/2*as.numeric(coil.channels)[(2*J-1):length(coil.channels)]+rep(site.eff.l[J],4)+rnorm(4,0,sd.l.meas))

          age.vec.phase2<-c(age.vec.phase2,rep(age[j],4))
          sex.vec.phase2<-c(sex.vec.phase2,rep(sex[j],4))
          dd.vec.phase2<-c(dd.vec.phase2,rep(dd[j],4))
          scanner.bias.vec.phase2<-c(scanner.bias.vec.phase2,rep(scanner.models[J],4))
          id.vec.phase2<-c(id.vec.phase2,rep(id.phase2[j],4))
          site.fac.phase2<-c(site.fac.phase2,rep(J,4))
          coil.fac.phase2<-c(coil.fac.phase2,coil.channels[(2*J-1):length(coil.channels)])

        } else {
          t.obs.phase2<-c(t.obs.phase2,rep(t[j],2)+bias.par*rep(scanner.models[this.visit.site],2)+bias.par/2*as.numeric(coil.channels)[(2*(this.visit.site)-1):(2*(this.visit.site))]+rep(site.eff.t[this.visit.site],2)+rnorm(2,0,sd.t.meas))
          l.obs.phase2<-c(l.obs.phase2,rep(l[j],2)+bias.par*rep(scanner.models[this.visit.site],2)+bias.par/2*as.numeric(coil.channels)[(2*(this.visit.site)-1):(2*(this.visit.site))]+rep(site.eff.l[this.visit.site],2)+rnorm(2,0,sd.l.meas))

          age.vec.phase2<-c(age.vec.phase2,rep(age[j],2))
          sex.vec.phase2<-c(sex.vec.phase2,rep(sex[j],2))
          dd.vec.phase2<-c(dd.vec.phase2,rep(dd[j],2))
          scanner.bias.vec.phase2<-c(scanner.bias.vec.phase2,rep(scanner.models[this.visit.site],2))
          id.vec.phase2<-c(id.vec.phase2,rep(id.phase2[j],2))
          site.fac.phase2<-c(site.fac.phase2,rep(this.visit.site,2))
          coil.fac.phase2<-c(coil.fac.phase2,coil.channels[(2*(this.visit.site)-1):(2*(this.visit.site))])

        }
      }
    }

    scanner.bias.vec.phase2<-factor(scanner.bias.vec.phase2)
    coil.fac.phase2<-factor(coil.fac.phase2)

    # Fit models and compare with KR-approximation based test for scanner effects
		m.small <- lmer(t.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2))
		m.large <- lmer(t.obs.phase2 ~ scanner.bias.vec.phase2 + sex.vec.phase2 + dd.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2))
		t.p.scanner[b]<-KRmodcomp(m.large, m.small)$stats$p.value
		
		m.small <- lmer(l.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2))
		m.large <- lmer(l.obs.phase2 ~ scanner.bias.vec.phase2 +  sex.vec.phase2 + dd.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2))
		l.p.scanner[b]<-KRmodcomp(m.large, m.small)$stats$p.value

    # Fit models and compare with KR-approximation based test for coil effects in the presence of scanner effects
    m.small <- lmer(t.obs.phase2 ~ scanner.bias.vec.phase2 + sex.vec.phase2 + dd.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2))
    m.large <- lmer(t.obs.phase2 ~ scanner.bias.vec.phase2 + coil.fac.phase2 + sex.vec.phase2 + dd.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2))
    t.p.coil[b]<-KRmodcomp(m.large, m.small)$stats$p.value
    
    m.small <- lmer(l.obs.phase2 ~ scanner.bias.vec.phase2 + sex.vec.phase2 + dd.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2))
    m.large <- lmer(l.obs.phase2 ~ scanner.bias.vec.phase2 + coil.fac.phase2 +  sex.vec.phase2 + dd.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2))
    l.p.coil[b]<-KRmodcomp(m.large, m.small)$stats$p.value

    setTxtProgressBar(pb, b) 
    }
	c(mean(t.p.scanner<0.05),mean(l.p.scanner<0.05),mean(t.p.coil<0.05),mean(l.p.coil<0.05))
}

# LOOK AT HOW SMALL THE MANUFACTURER BIASES CAN BE
set.seed(529349245)
mclapply(seq(0,1,by=0.1),get.power.phase2,mc.cores=20)
mclapply(seq(1,2,by=0.1),get.power.phase2,mc.cores=20)

# we can detect:
# for thalamus - scanner effects of 0.8mL or larger (under assumptions)
# for lesion - scanner effects of 1.2mL or larger (under assumptions)
# for thalamus - coil effects of 0.9 mL or larger (under assumptions)
# for lesion - coil effects of 1.1 mL or larger (under assumptions)

################################################
################################################
############## AIM 2 - LONGITUDINAL ANALYSES WITH DROPOUT
################################################
################################################

### AIM 2 - PHASE I POWER CALCULATIONS FOR SITE EFFECTS
get.power.phase1.longit.dropout<-function(sd.site.longit,longit.bias.par,B=500) {
  id<-factor(rep(101:102,2*(J+1)))  
  n.1<-2
  site<-factor(c(rep(1:J,each=2*n.1),rep(J,each=2*n.1)))
  l.p<-rep(NA,B);t.p<-rep(NA,B) # tests for proportion of variation explained by site
  l.comp.p<-rep(NA,B);t.comp.p<-rep(NA,B) # tests for proportion of variation explained by site
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in 1:B) {
    
    ## Here, we analyze the difference in volumes and assess site effects therein
    t<-2*t.age.beta + 2*t.dd.beta
    l<-2*l.dd.beta

    #Sample site effects - add NIH site again for second scan session
    site.eff.t<-rnorm(J,0,sd.site.longit);site.eff.t<-c(site.eff.t,site.eff.t[length(site.eff.t)])
    site.eff.l<-rnorm(J,0,sd.site.longit);site.eff.l<-c(site.eff.l,site.eff.l[length(site.eff.l)])

    #Specify acquisition biases
    scanner.bias<-longit.bias.par*c(scanner.models,1) # based on scanner model (2=Prisma, 1=Skyra, 0=Trio)
    coil.bias<-longit.bias.par/2*as.numeric(coil.channels) # conservatively, assume receive coils are responsible for an additional 50% of variation

    #Sample observed data with both scanner and coil effects
    t.obs<-rep(t,2*(J+1))+rep(scanner.bias,each=2*n.1)+rep(coil.bias,each=n.1)+rep(site.eff.t,each=2*n.1)+rnorm(2*(J+1)*n.1,0,sd.t.meas)
    l.obs<-rep(l,2*(J+1))+rep(scanner.bias,each=2*n.1)+rep(coil.bias,each=n.1)+rep(site.eff.l,each=2*n.1)+rnorm(2*(J+1)*n.1,0,sd.l.meas)

    #Do permutation tests
    t.p[b]<-perm.test.F(t.obs,site,id)
    l.p[b]<-perm.test.F(l.obs,site,id)

    setTxtProgressBar(pb, b) 
  }
  c(mean(t.p<0.05),mean(l.p<0.05))
}

# LOOK AT HOW SMALL THE MANUFACTURER BIASES CAN BE
set.seed(529349245)
mclapply(seq(0,1,by=0.1),get.power.phase1.longit.dropout,0,B=500,mc.cores=20)
switch.args<-function(longit.bias.par,sd.site.longit,B) get.power.phase1.longit.dropout(sd.site.longit,longit.bias.par,B)
mclapply(seq(0,1,by=0.1),switch.args,0,B=500,mc.cores=20)

# we can detect:
# for thalamus - longitudinal scanner platform effects of 0.4mL or larger (under assumptions)
# for lesion - longitudinal scanner platform effects of 0.5mL or larger (under assumptions)
# for thalamus - intersite differences with sd of 0.4mL or larger (under assumptions)
# for lesion - intersite differences with sd of 0.4mL or larger (under assumptions)


### AIM 2 - PHASE II POWER CALCULATIONS FOR SCANNER AND COIL EFFECTS
get.power.phase2.longit.dropout<-function(longit.bias.par,bias.par,B=500) {  

  l.p.scanner<-rep(NA,B);t.p.scanner<-rep(NA,B) # tests for variation explained by scanner
  l.p.coil<-rep(NA,B);t.p.coil<-rep(NA,B) # tests for variation explained by coil
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  for (b in 1:B) {
  
    #Sample patient-level variables
    id.phase2<-1:n.2
    age<-runif(n.2,min.age,max.age)
    onset<-pmax(age-runif(n.2,min.dd,max.dd),10)
    dd<-age-onset
    sex<-rbinom(n.2,1,p.female)
    time.in.study<-rnorm(n.2,2,0.2)
    t.baseline<-rnorm(n.2,t.vol,sd.t.vol) + (age-50)*t.age.beta + sex*t.sex.beta + dd*t.dd.beta
    l.baseline<-rnorm(n.2,l.vol,sd.l.vol) + sex*l.sex.beta + dd*l.dd.beta
    t.followup<-t.baseline + time.in.study*t.age.beta + time.in.study*t.dd.beta
    l.followup<-l.baseline + time.in.study*l.dd.beta

    #Sample site variable
    site<-NULL
    home.sites<-rep(1:J,each=n.2/J)
    for (j in 1:n.2) {
      site<-c(site,home.sites[j],sample(1:J,2))
    }

    #Sample imaging biases
    site.eff.t<-rnorm(J,0,sd.t.site)
    site.eff.l<-rnorm(J,0,sd.l.site)
    bias.fac<-factor(scanner.models)
    coil.bias.vec<-rep(NA,2*length(site))
    for (j in 1:length(site)) coil.bias.vec[(2*j-1):(2*j)]<-coil.channels[which(rep(1:J,each=2)==site[j])]    

    t.obs.phase2.baseline<-NULL
    l.obs.phase2.baseline<-NULL
    t.obs.phase2.followup<-NULL
    l.obs.phase2.followup<-NULL
    age.vec.phase2<-NULL
    sex.vec.phase2<-NULL
    dd.vec.phase2<-NULL
    scanner.vec.phase2<-NULL
    id.vec.phase2<-NULL
    site.fac.phase2<-NULL
    coil.fac.phase2<-NULL
    time.in.study.vec<-NULL

    
    for (j in 1:n.2) { # for each subject
      for (k in 1:m) { # for each site that subject is imaged at

        this.visit.site<-site[rep(id.phase2,each=3)==j][k]
        # if NIH
        if (this.visit.site==J) {

          #Sample observed data
          t.obs.phase2.baseline<-c(t.obs.phase2.baseline,rep(t.baseline[j],4)+bias.par*rep(scanner.models[J],4)+bias.par/2*as.numeric(coil.channels)[(2*J-1):length(coil.channels)]+rep(site.eff.t[J],4)+rnorm(4,0,sd.t.meas))
          l.obs.phase2.baseline<-c(l.obs.phase2.baseline,rep(l.baseline[j],4)+bias.par*rep(scanner.models[J],4)+bias.par/2*as.numeric(coil.channels)[(2*J-1):length(coil.channels)]+rep(site.eff.l[J],4)+rnorm(4,0,sd.l.meas))
          t.obs.phase2.followup<-c(t.obs.phase2.followup,rep(t.followup[j],4)+bias.par*rep(scanner.models[J],4)+longit.bias.par*rep(time.in.study[j],4)*rep(scanner.models[J],4)+bias.par/2*as.numeric(coil.channels)[(2*J-1):length(coil.channels)]+longit.bias.par/2*rep(time.in.study[j],4)*as.numeric(coil.channels)[(2*J-1):length(coil.channels)]+rep(site.eff.t[J],4)+rnorm(4,0,sd.t.meas))
          l.obs.phase2.followup<-c(l.obs.phase2.followup,rep(l.followup[j],4)+bias.par*rep(scanner.models[J],4)+longit.bias.par*rep(time.in.study[j],4)*rep(scanner.models[J],4)+bias.par/2*as.numeric(coil.channels)[(2*J-1):length(coil.channels)]+longit.bias.par/2*rep(time.in.study[j],4)*as.numeric(coil.channels)[(2*J-1):length(coil.channels)]+rep(site.eff.l[J],4)+rnorm(4,0,sd.l.meas))

          age.vec.phase2<-c(age.vec.phase2,rep(age[j],4))
          sex.vec.phase2<-c(sex.vec.phase2,rep(sex[j],4))
          dd.vec.phase2<-c(dd.vec.phase2,rep(dd[j],4))
          scanner.vec.phase2<-c(scanner.vec.phase2,rep(scanner.models[this.visit.site],4))
          id.vec.phase2<-c(id.vec.phase2,rep(id.phase2[j],4))
          site.fac.phase2<-c(site.fac.phase2,rep(J,4))
          coil.fac.phase2<-c(coil.fac.phase2,coil.channels[(2*J-1):length(coil.channels)])
          time.in.study.vec<-c(time.in.study.vec,rep(time.in.study[j],4))

        } else {

          #Sample observed data
          t.obs.phase2.baseline<-c(t.obs.phase2.baseline,rep(t.baseline[j],2)+bias.par*rep(scanner.models[this.visit.site],2)+bias.par/2*as.numeric(coil.channels)[(2*(this.visit.site)-1):(2*(this.visit.site))]+rep(site.eff.t[this.visit.site],2)+rnorm(2,0,sd.t.meas))
          l.obs.phase2.baseline<-c(l.obs.phase2.baseline,rep(l.baseline[j],2)+bias.par*rep(scanner.models[this.visit.site],2)+bias.par/2*as.numeric(coil.channels)[(2*(this.visit.site)-1):(2*(this.visit.site))]+rep(site.eff.l[this.visit.site],2)+rnorm(2,0,sd.l.meas))
          t.obs.phase2.followup<-c(t.obs.phase2.followup,rep(t.followup[j],2)+bias.par*rep(scanner.models[this.visit.site],2)+longit.bias.par*rep(time.in.study[j],2)*rep(scanner.models[this.visit.site],2)+bias.par/2*as.numeric(coil.channels)[(2*(this.visit.site)-1):(2*(this.visit.site))]+longit.bias.par/2*rep(time.in.study[j],2)*as.numeric(coil.channels)[(2*(this.visit.site)-1):(2*(this.visit.site))]+rep(site.eff.t[this.visit.site],2)+rnorm(2,0,sd.t.meas))
          l.obs.phase2.followup<-c(l.obs.phase2.followup,rep(l.followup[j],2)+bias.par*rep(scanner.models[this.visit.site],2)+longit.bias.par*rep(time.in.study[j],2)*rep(scanner.models[this.visit.site],2)+bias.par/2*as.numeric(coil.channels)[(2*(this.visit.site)-1):(2*(this.visit.site))]+longit.bias.par/2*rep(time.in.study[j],2)*as.numeric(coil.channels)[(2*(this.visit.site)-1):(2*(this.visit.site))]+rep(site.eff.l[this.visit.site],2)+rnorm(2,0,sd.l.meas))

          age.vec.phase2<-c(age.vec.phase2,rep(age[j],2))
          sex.vec.phase2<-c(sex.vec.phase2,rep(sex[j],2))
          dd.vec.phase2<-c(dd.vec.phase2,rep(dd[j],2))
          scanner.vec.phase2<-c(scanner.vec.phase2,rep(scanner.models[this.visit.site],2))
          id.vec.phase2<-c(id.vec.phase2,rep(id.phase2[j],2))
          site.fac.phase2<-c(site.fac.phase2,rep(this.visit.site,2))
          coil.fac.phase2<-c(coil.fac.phase2,coil.channels[(2*(this.visit.site)-1):(2*(this.visit.site))])
          time.in.study.vec<-c(time.in.study.vec,rep(time.in.study[j],2))

        }
      }
    }

    #Re-organize
    t.obs.phase2<-merge.vecs(t.obs.phase2.baseline,t.obs.phase2.followup); l.obs.phase2<-merge.vecs(l.obs.phase2.baseline,l.obs.phase2.followup)
    age.vec.phase2<-rep(age.vec.phase2,each=2);sex.vec.phase2<-rep(sex.vec.phase2,each=2); dd.vec.phase2<-rep(dd.vec.phase2,each=2); scanner.vec.phase2<-rep(scanner.vec.phase2,each=2); id.vec.phase2<-factor(rep(id.vec.phase2,each=2)); site.fac.phase2<-factor(rep(site.fac.phase2,each=2)); coil.fac.phase2<-factor(rep(coil.fac.phase2,each=2))
    time.in.study.vec<-merge.vecs(rep(0,length(time.in.study.vec)),time.in.study.vec)
    visit.id<-merge.vecs(rep(0,length(t.obs.phase2.baseline)),rep(1,length(t.obs.phase2.followup)))

    #Sample dropout
    full.df<-data.frame(l.obs.phase2,t.obs.phase2,scanner.vec.phase2,sex.vec.phase2,dd.vec.phase2,time.in.study.vec,site.fac.phase2,id.vec.phase2,coil.fac.phase2)
    which.lost.to.fu<-sample(1:n.2,size=0.3*n.2)
    which.fu<-which(visit.id==1)
    which.fu.obs<-which.fu[!id.vec.phase2[which.fu]%in%which.lost.to.fu]
    dropout.df<-full.df[c(which(visit.id==0),which.fu.obs),]

    # Fit models and compare with KR-approximation based test for scanner effects
    m.small <- lmer(t.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + time.in.study.vec + scanner.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2),data=dropout.df)
    m.large <- lmer(t.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + time.in.study.vec*scanner.vec.phase2  + (1| id.vec.phase2) + (1| site.fac.phase2),data=dropout.df)
    t.p.scanner[b]<-KRmodcomp(m.large, m.small)$stats$p.value
    
    m.small <- lmer(l.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + time.in.study.vec + scanner.vec.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2),data=dropout.df)
    m.large <- lmer(l.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + time.in.study.vec*scanner.vec.phase2+ (1| id.vec.phase2) + (1| site.fac.phase2),data=dropout.df)
    l.p.scanner[b]<-KRmodcomp(m.large, m.small)$stats$p.value

    # Fit models and compare with KR-approximation based test for coil effects in the presence of scanner effects
    m.small <- lmer(t.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + time.in.study.vec*scanner.vec.phase2  + coil.fac.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2),data=dropout.df)
    m.large <- lmer(t.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + time.in.study.vec*scanner.vec.phase2  + time.in.study.vec*coil.fac.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2),data=dropout.df)
    t.p.coil[b]<-KRmodcomp(m.large, m.small)$stats$p.value
    
    m.small <- lmer(l.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + time.in.study.vec*scanner.vec.phase2  + coil.fac.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2),data=dropout.df)
    m.large <- lmer(l.obs.phase2 ~ sex.vec.phase2 + dd.vec.phase2 + time.in.study.vec*scanner.vec.phase2  + time.in.study.vec*coil.fac.phase2 + (1| id.vec.phase2) + (1| site.fac.phase2),data=dropout.df)
    l.p.coil[b]<-KRmodcomp(m.large, m.small)$stats$p.value

      setTxtProgressBar(pb, b) 
    }
  c(mean(t.p.scanner<0.05),mean(l.p.scanner<0.05),mean(t.p.coil<0.05),mean(l.p.coil<0.05))
}

set.seed(529349245)
mclapply(seq(0,0.5,by=0.05),get.power.phase2.longit.dropout,bias.par=0,B=500,mc.cores=20)
mclapply(seq(0,0.5,by=0.05),get.power.phase2.longit.dropout,bias.par=1,B=500,mc.cores=20)
mclapply(seq(0.5,0.9,by=0.05),get.power.phase2.longit.dropout,bias.par=0,B=500,mc.cores=20)
mclapply(seq(0.5,0.9,by=0.05),get.power.phase2.longit.dropout,bias.par=1,B=500,mc.cores=20)
