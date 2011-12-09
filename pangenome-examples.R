## This is an example script for calculating
## properties of pangenomes and core genomes
## using the Infinitely Many Genes model of
## Pfaffelhuber and Baumdicker 2010
## with modifications by Collins and Higgs (2011)
##
## written by R. Eric Collins, 2011
## released under the GNU Public License
## The most recent version can be found at
## http://github.com/rec3141/pangenome
##
## Usage:
## If you have R properly installed you
## should be able to run the examples from the
## command line via: 
## R CMD BATCH pangenome-examples.R
## then look at the plots produced in Rout.pdf.
##
## From within R, run the examples via:
## source("pangenome-examples.R")
##
## To work with your own dataset, make sure
## to realize that the fitting routines are
## sensitive to initial parameters, so you 
## will likely need to change them
## to find the best fit

rm(list=ls())

#load file containing functions
source("f-pangenome.R")

# read in matrix of all Bacilli gene clusters
# less strict clustering
infile <- "Bacilli-clusters.1e-10.6.csv"
# more strict clustering
# infile <- "Bacilli-clusters.1e-30.7.csv"
mat.all <- read.csv(infile,sep="\t")

mat.read <- mat.all[,3+1:(dim(mat.all)[2]-7)]
keys <- read.table("key-Bacillaceae",sep="\t",row.names=1)
mycols <- unique(unlist(lapply(keys[,2],function(x) grep(x,colnames(mat.read)))))
mat <- mat.read[,mycols]
mat <- mat[rowSums(mat)>0,]

genomesize <- mean(colSums(mat>0)) # mean genome size measured in gene families
ng <- dim(mat)[2] #number of genomes

# we're left with a matrix of clusters 
# containing gene families in the genomes of interest

##---------------------------------
## Gene Family Frequency Spectrum
##---------------------------------

#calculate gene family frequency spectrum
Gk <- f.getspectrum(mat)

#plot gene family frequency spectrum
bp.x <- barplot(Gk, ylab="Frequency of Gene Families",xlab="Number of Genomes",names.arg=1:ng,main="Bacillaceae 1\nGene Family Frequency Spectrum",log='y')

## Predict the gene family frequency spectrum
## using model 1D+E on a coalescent

mytree<-"coalescent.spec" #use the coalescent tree w/G(k)
myfitting<-"chi2" #fit it using the Chi^2
constr<- 1 # G0 constrained to the mean genome size during fitting
mymethod <- "Nelder" # alternative fitting routine: "BFGS"

params.cde <- c(1,100) #starting parameters theta1, ness

opt.cde <- optim(
  params.cde,
  function(x) f.fit(x,constr=constr,data=Gk,treetype=mytree,fitting=myfitting,genomesize=genomesize,ng=ng),
  method=mymethod,
  control=list(trace=1,parscale=params.cde+1e-6,maxit=10000)
)

# best fit for non-strict dataset:
# chi2 = 21589.43
# params = 1.037437 842.326605

# best fit for strict dataset:
# chi2 = 9088.489462
# params = 0.9686929 854.3339099

# get optimized parameters, calculate constrained parameter
params.cde <- c(opt.cde$par[1],opt.cde$par[1]*genomesize,opt.cde$par[2])
print(params.cde)
spec.cde <- f.coalescent.spec(params.cde,ng)
lines(bp.x,spec.cde,col='red',lwd=2,lty=2)


## Predict the gene family frequency spectrum
## using model 2D+E on a coalescent

mytree<-"coalescent.spec" #use the coalescent tree w/G(k)
myfitting<-"chi2" #fit it using the Chi^2
constr<- 1 # G0 constrained to the mean genome size during fitting
mymethod <- "Nelder" # alternative fitting routine: "BFGS"

params.c2de <- c(1,1000,1000,100) #rho1, theta1, ness, rho2, (theta2 constrained)

opt.c2de <- optim(
  params.c2de,
  function(x) f.fit(x,constr=constr,data=Gk,treetype=mytree,fitting=myfitting,genomesize=genomesize,ng=ng),
  method=mymethod,
  control=list(trace=1,parscale=params.c2de+1e-6,maxit=10000)
)

# best fit for non-strict dataset:
# chi2 = 1007.792
# params = 1.047101 2064.168648 1103.227536  200.462348

# best fit for strict dataset:
# chi2 = 593.373075
# params = 0.952658  1488.599551  1037.741085    75.78778

# get optimized parameters, calculate constrained parameter
params.c2de <- c(opt.c2de$par,opt.c2de$par[4]*(genomesize-opt.c2de$par[2]/opt.c2de$par[1]-opt.c2de$par[3]))
print(params.c2de)
spec.c2de <- f.coalescent.spec(params.c2de,ng)
lines(bp.x,spec.c2de,col='red',lwd=2)


## Predict the gene family frequency spectrum
## using model 2D+E on a fixed tree

# load tree table
treetable<-read.table("treetable-Bacillaceae",sep="\t",row.names=1)
colnames(treetable) <- c("desc_a","desc_b","dist")

mytree<-"fixed.spec" #use a fixed tree w/G(k)
myfitting<-"chi2" #fit it using chi2 or sumsq
constr<- 1 # G0 constrained to the mean genome size during fitting
mymethod <- "BFGS" # alternative fitting routines: "BFGS","Nelder"

params.f2de <- opt.c2de$par #rho1, theta1, ness, rho2, (theta2 constrained)

opt.f2de <- optim(
  params.f2de,
  function(x) f.fit(x,constr=constr,data=Gk,treetype=mytree,fitting=myfitting,genomesize=genomesize,ng=ng,treetable=treetable),
  method=mymethod,
  control=list(trace=1,parscale=params.f2de+1e-6,maxit=10000)
)

# best fit for non-strict dataset:
# chi2 = 1104.552
# params = 1.194561 2076.469124 1062.999820  236.381120

# best fit for strict dataset:
# chi2 = 495.020973
# params = 1.271341  1735.432410  1013.196002    70.250007


# get optimized parameters, calculate constrained parameter
params.f2de <- c(opt.f2de$par,opt.f2de$par[4]*(genomesize-opt.f2de$par[2]/opt.f2de$par[1]-opt.f2de$par[3]))
print(params.f2de)
spec.f2de <- f.fixed.spec(params.f2de,treetable)
lines(bp.x,spec.f2de,col='blue',lwd=2)

legend("top",c("1D+E on coalescent","2D+E on coalescent","2D+E on fixed tree"),lty=c(2,1,1),col=c('red','red','blue'),lwd=2)
dev.copy2pdf(file="genespectrum.pdf",width=6,height=6)

##---------------------------------
## Pangenome and core genome curves
##---------------------------------

# Calculate 100 permutations of the pangenome and core genome
perm.pangenome <- f.pangenome(mat,100)
perm.core <- f.core(mat,100)

# Calculate the exact mean pan and core genome curves
# from the gene frequency spectrum G(k)
mean.pangenome <- f.meanpancore(Gk)$pan
mean.core <- f.meanpancore(Gk)$core
pancore <- c(mean.pangenome,mean.core)

# Prepare a new plot window
plot(1:ng,xlim=c(1,ng),ylim=c(0.9*min(mean.core), 1.1*max(mean.pangenome)),log='y',xlab="Genomes added", ylab="Clusters of Gene Families",main="Bacillaceae 1\nPangenome and Core genome",pch='')

# Plot polygons outlining permutations
polygon(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))), col="gray88",border=NA)
polygon(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.pangenome,1,max))), col="gray88",border=NA)

# Add the mean pan and core genome curves to the plot
points(1:ng,mean.pangenome)
points(1:ng,mean.core)

## Predict the core and pangenome curves
## using model 2D+E on a coalescent tree

mytree<-"coalescent" #use the coalescent tree w/G(k)
myfitting<-"chi2" #fit it using the Chi^2
constr<- 1 # G0 constrained to the mean genome size during fitting
mymethod <- "Nelder" # alternative fitting routine: "Nelder"

params.c2de <- c(1,1000,1000,10) #rho1, theta1, ness, rho2, (theta2 constrained)
params.c2de <- c(1,1000,1000,5) #rho1, theta1, ness, rho2, (theta2 constrained)

opt.c2de <- optim(
  params.c2de,
  function(x) f.fit(x,constr=constr,data=pancore,treetype=mytree,fitting=myfitting,genomesize=genomesize,ng=ng),
  method=mymethod,
  control=list(trace=1,parscale=params.c2de+1e-6,maxit=10000)
)

# best fit for non-strict dataset:
# chi2 = 6.216309
# params = 1.011401 2216.125953  896.717547  204.852939

# best fit for strict dataset:
# chi2 = 7.569648
# params = 0.8522393 1432.3099454  892.6032901   56.2204857

params.c2de <- c(opt.c2de$par,opt.c2de$par[4]*(genomesize-opt.c2de$par[2]/opt.c2de$par[1]-opt.c2de$par[3]))
print(params.c2de)
pancore.c2de <- f.coalescent(params.c2de,ng)

# Add the predicted curves to the plot
lines(1:ng,pancore.c2de$pan,col='red',lwd=2)
lines(1:ng,pancore.c2de$core,col='red',lwd=2)

## Predict the core and pangenome curves
## using model 2D+E on a star tree

mytree<-"star" #use the coalescent tree w/G(k)
myfitting<-"chi2" #fit it using the Chi^2
constr<- 1 # G0 constrained to the mean genome size during fitting
mymethod <- "BFGS" # alternative fitting routine: "Nelder"

params.s2de <- c(0.01,10,1000,1) #rho1, theta1, ness, rho2, (theta2 constrained)

opt.s2de <- optim(
  params.s2de,
  function(x) f.fit(x,constr=constr,data=pancore,treetype=mytree,fitting=myfitting,genomesize=genomesize,ng=ng),
  method=mymethod,
  control=list(trace=1,parscale=params.s2de+1e-6,maxit=10000)
)

# best fit for non-strict dataset: 
# chi2 = 3640.319110
# params = 0.002859687 2.996654750 5.029033100 0.365977668

# best fit for strict dataset:
# chi2 = 4306.126837
# params = 7.746796e-04 7.737486e-01 3.812043e-01 2.636258e-01

params.s2de <- c(opt.s2de$par,opt.s2de$par[4]*(genomesize-opt.s2de$par[2]/opt.s2de$par[1]-opt.s2de$par[3]))
print(params.s2de)
pancore.s2de <- f.star(params.s2de,ng)

# Add the predicted curves to the plot
lines(1:ng,pancore.s2de$pan,col='blue',lwd=2,lty=3)
lines(1:ng,pancore.s2de$core,col='blue',lwd=2,lty=3)


## Predict the core and pangenome curves
## using model 2D+E on a fixed tree
## (calculated from the gene frequency spectrum)

mytree<-"fixed" #use the coalescent tree w/G(k)
myfitting<-"chi2" #fit it using the Chi^2
constr<- 1 # G0 constrained to the mean genome size during fitting
mymethod <- "BFGS" # alternative fitting routine: "BFGS"

params.f2de <- opt.c2de$par #rho1, theta1, ness, rho2, (theta2 constrained)

opt.f2de <- optim(
  params.f2de,
  function(x) f.fit(x,constr=constr,data=pancore,treetype=mytree,fitting=myfitting,treetable=treetable,genomesize=genomesize,ng=ng),
  method=mymethod,
  control=list(trace=1,parscale=params.f2de+1e-6,maxit=10000)
)

# best fit for non-strict dataset:
# chi2 = 8.818609
# params = 1.564745 3336.309 914.9655 947.0657

# best fit for strict dataset
# chi2 = 2.795111
# params = 1.078041  1493.368589   867.686765    37.312469

params.f2de <- c(opt.f2de$par,opt.f2de$par[4]*(genomesize-opt.f2de$par[2]/opt.f2de$par[1]-opt.f2de$par[3]))
print(params.f2de)
pancore.f2de <- f.meanpancore(f.fixed.spec(params.f2de,treetable))

lines(1:ng,pancore.f2de$pan,col='orange',lwd=2,lty=1)
lines(1:ng,pancore.f2de$core,col='orange',lwd=2,lty=1)

legend("right",c("2D+E on coalescent","2D+E on star tree","2D+E on a fixed tree"),lty=c(1,3,1),col=c('red','blue','orange'),lwd=2)
dev.copy2pdf(file="pancore.pdf",width=6,height=6)
