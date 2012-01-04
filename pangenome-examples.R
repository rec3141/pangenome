## This is an example script for calculating
## properties of pangenomes and core genomes
## and for fitting to user generated data
## using the Infinitely Many Genes model of
## Baumdicker and Pfaffelhuber (2010)
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
## then look at the plots produced in Rplots.pdf.
##
## From within R, run the example script via:
## source("pangenome-examples.R")
##
## To work with your own dataset, make sure
## to realize that the fitting routines are
## sensitive to initial parameters, so you 
## will likely need to adjust them
## to find the best fit
##
## you will also need to run tre2table.pl
## on your phylogenetic trees to get the treetable
## needed for fitting fixed trees

rm(list=ls())

#load file containing functions
source("f-pangenome.R")

# read in matrix of all Bacilli gene clusters
# less strict clustering (default):
# infile <- "Bacilli-clusters.1e-10.6.csv"
# more strict clustering:
infile <- "Bacilli-clusters.1e-30.7.csv"
mat.all <- read.csv(infile,sep="\t")
mat.read <- mat.all[,3+1:(dim(mat.all)[2]-7)] #remove annotations

# get list of key files to process
allkeys <- list.files(path="./",pattern="key-.+")

table1 <- NULL
table2 <- NULL
table3 <- NULL

for (mykey in allkeys) {
  # the key matches the names to RefSeq IDs
  keys <- read.table(mykey,sep="\t",row.names=1)
  # the taxaname is the name of the taxonomic group
  taxaname <- unlist(strsplit(mykey,'-'))[2]
  print(taxaname)
  # get the columns which are of interest
  mycols <- unique(unlist(lapply(keys[,2],function(x) grep(x,colnames(mat.read)))))
  # get a matrix of clusters containing gene families in the genomes of interest
  mat <- mat.read[rowSums(mat.read)>0,mycols]

  genomesize <- mean(colSums(mat>0)) # mean genome size measured in gene families
  ng <- dim(mat)[2] #number of genomes

  # load phylogenetic tree table for this taxonomic group
  treename <- paste("treetable",taxaname,sep="-")
  treetable <- read.table(treename,sep="\t",row.names=1)
  colnames(treetable) <- c("desc_a","desc_b","dist")

  ##---------------------------------
  ## Table 1 -- genome properties
  ##---------------------------------

  Ngenes <- mean(colSums(mat)) # mean number of genes
  G0 <- mean(colSums(mat>0)) # mean number of gene families
  Gcore <- length(which(rowSums(mat>0)==ng)) # mean core size
  Gpan <- nrow(mat) #pangenome size
  dprot <- f.dprot(treetable,ng) # mean distance from LCA

  # add results to Table1
  table1 <- rbind(table1,cbind(
    taxaname,
    ng,
    Ngenes,
    G0,
    Ngenes/G0,
    Gcore,
    Gcore/G0,
    Gpan,
    Gpan/G0,
    dprot))

  ##-------------------------------------------
  ## Fit to the Gene Family Frequency Spectrum
  ##-------------------------------------------

  #calculate gene family frequency spectrum
  Gk <- f.getspectrum(mat)

  #plot gene family frequency spectrum
  bp.x <- barplot(Gk, ylab="Frequency of Gene Families",xlab="Number of Genomes",names.arg=1:ng,main=paste(taxaname,"Gene Family Frequency Spectrum",sep="\n"),log='y')
  legend("top",c("1D+E on coalescent","2D+E on coalescent","2D+E on fixed tree"),lty=c(2,1,1),col=c('red','red','blue'),lwd=2)

  ## Predict the gene family frequency spectrum
  ## using model 1D+E on a coalescent

  mymodel <- "coalescent.spec" #use the coalescent tree w/G(k)
  myfitting <- "chi2" #fit it using the Chi^2
  constr <- 1 # G0 constrained to the mean genome size during fitting
  # !! if constr changes, need to change input parameters as well
  mymethod <- "BFGS" # alternative fitting routine: "BFGS"

  #set initial parameters and recursively optimize
  opt.spec.cde <- f.recurse(c(1,100),r.data=Gk)

  # get optimized parameters, calculate constrained parameter
  params.spec.cde <- c(opt.spec.cde$par[1],opt.spec.cde$par[1]*(genomesize-opt.spec.cde$par[2]),opt.spec.cde$par[2])

  print(params.spec.cde)
  spec.cde <- f.coalescent.spec(params.spec.cde,ng)
  lines(bp.x,spec.cde,col='red',lwd=2,lty=2)

  ## Predict the gene family frequency spectrum
  ## using model 2D+E on a coalescent

  mymodel <- "coalescent.spec" #use the coalescent tree w/G(k)
  myfitting <- "chi2" #fit it using the Chi^2
  constr <- 1 # G0 constrained to the mean genome size during fitting
  # !! if constr changes, need to change input parameters as well
  mymethod <- "BFGS" # alternative fitting routine: "BFGS"

  #set initial parameters and recursively optimize
  pinitial.spec.c2de <- c(params.spec.cde[1]/10,params.spec.cde[2]/100,params.spec.cde[3],params.spec.cde[1])
  opt.spec.c2de <- f.recurse(pinitial.spec.c2de,r.data=Gk)

  # get optimized parameters, calculate constrained parameter
  params.spec.c2de <- c(opt.spec.c2de$par,opt.spec.c2de$par[4]*(genomesize-opt.spec.c2de$par[2]/opt.spec.c2de$par[1]-opt.spec.c2de$par[3]))
  print(params.spec.c2de)
  spec.c2de <- f.coalescent.spec(params.spec.c2de,ng)
  lines(bp.x,spec.c2de,col='red',lwd=2)

  ## Predict the gene family frequency spectrum
  ## using model 2D+E on a fixed tree

  mymodel <- "fixed.spec" #use a fixed tree w/G(k)
  myfitting <- "chi2" #fit it using chi2 or sumsq
  constr <- 1 # G0 constrained to the mean genome size during fitting;
  # !! if constr changes, need to change input parameters as well
  mymethod <- "BFGS" # alternative fitting routines: "BFGS","Nelder"

  #set initial parameters and recursively optimize
  pinitial.spec.f2de <- opt.spec.c2de$par/c(10,10,1,1)
  opt.spec.f2de <- f.recurse(pinitial.spec.f2de,r.data=Gk)

  # get optimized parameters, calculate constrained parameter
  params.spec.f2de <- c(opt.spec.f2de$par,opt.spec.f2de$par[4]*(genomesize-opt.spec.f2de$par[2]/opt.spec.f2de$par[1]-opt.spec.f2de$par[3]))
  print(params.spec.f2de)
  spec.f2de <- f.fixed.spec(params.spec.f2de,treetable)
  lines(bp.x,spec.f2de,col='blue',lwd=2)

  # print to PDF if interactive
  if(.Internal(dev.displaylist())) {
    dev.copy2pdf(file=paste("genespectrum-",taxaname,".pdf",sep=""),width=6,height=6)
  }

  ##-----------------------------------------------------
  ## Fit to the Pangenome and core genome curves
  ##-----------------------------------------------------

  # Calculate 100 permutations each of the pangenome and core genome
  perm.pangenome <- f.pangenome(mat,100)
  perm.core <- f.core(mat,100)

  # Calculate the exact mean pan and core genome curves
  # from the gene frequency spectrum G(k)
  mean.pangenome <- f.meanpancore(Gk)$pan
  mean.core <- f.meanpancore(Gk)$core
  pancore <- c(mean.pangenome,mean.core)

  # Calculate the RMS value for the permutations
  rms.perm <- mean(f.rms(c(mean.pangenome,mean.core),rbind(perm.pangenome,perm.core)))

  # Prepare a new plot window
  plot(1:ng,xlim=c(1,ng),ylim=c(0.9*min(mean.core), 1.1*max(mean.pangenome)),log='y',xlab="Genomes added", ylab="Clusters of Gene Families",main=paste(taxaname,"Pangenome and Core genome",sep="\n"),pch='')
  legend("right",c("2D+E on coalescent","2D+E on star tree","2D+E on a fixed tree"),lty=c(1,3,1),col=c('red','blue','orange'),lwd=2)

  # Plot polygons outlining permutations
  polygon(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))), col="gray88",border=NA)
  polygon(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.pangenome,1,max))), col="gray88",border=NA)

  # Add the mean pan and core genome curves to the plot
  points(1:ng,mean.pangenome)
  points(1:ng,mean.core)

  ## Predict the core and pangenome curves
  ## using model 1D+E on a coalescent tree

  mymodel <- "coalescent" #use the coalescent tree w/pancore
  myfitting <- "chi2" #fit it using the Chi^2
  constr <- 1 # G0 constrained to the mean genome size during fitting
  mymethod <- "Nelder" # alternative fitting routine: "Nelder"

  # recursively optimize using initial parameters and dataset to fit
  opt.cde <- f.recurse(opt.spec.cde$par,r.data=pancore)

  # get optimized parameters and calculate constrained param
  params.cde <- c(opt.cde$par[1], opt.cde$par[1]*(genomesize-opt.cde$par[2]), opt.cde$par[2])
  print(params.cde)
  pancore.cde <- f.coalescent(params.cde,ng)
  rms.cde <- f.rms(c(mean.pangenome,mean.core),c(pancore.cde$pan,pancore.cde$core))

  # Add the predicted curves to the plot
#   lines(1:ng,pancore.cde$pan,col='red',lty=2,lwd=2)
#   lines(1:ng,pancore.cde$core,col='red',lty=2,lwd=2)

  ## Predict the core and pangenome curves
  ## using model 2D+E on a coalescent tree

  mymodel <- "coalescent" #use the coalescent tree w/G(k)
  myfitting <- "chi2" #fit it using the Chi^2
  constr <- 1 # G0 constrained to the mean genome size during fitting
  mymethod <- "Nelder" # alternative fitting routine: "Nelder"

  # recursively optimize using initial parameters and dataset to fit
  opt.c2de <- f.recurse(opt.spec.c2de$par,r.data=pancore)

  # get optimized parameters and calculate constrained param
  params.c2de <- c(opt.c2de$par,opt.c2de$par[4]*(genomesize-opt.c2de$par[2]/opt.c2de$par[1]-opt.c2de$par[3]))
  print(params.c2de)
  pancore.c2de <- f.coalescent(params.c2de,ng)
  rms.c2de <- f.rms(c(mean.pangenome,mean.core),c(pancore.c2de$pan,pancore.c2de$core))

  # Add the predicted curves to the plot
  lines(1:ng,pancore.c2de$pan,col='red',lwd=2)
  lines(1:ng,pancore.c2de$core,col='red',lwd=2)

  ## Predict the core and pangenome curves
  ## using model 2D+E on a star tree

  mymodel <- "star" #use the coalescent tree w/G(k)
  myfitting <- "chi2" #fit it using the Chi^2
  constr <- 1 # G0 constrained to the mean genome size during fitting
  mymethod <- "Nelder" # alternative fitting routine: "BFGS","Nelder"

  # recursively optimize using initial parameters and dataset to fit
  pinitial.s2de <- opt.c2de$par/c(10,10,1,1)
  opt.s2de <- f.recurse(pinitial.s2de,r.data=pancore)

  # get optimized parameters and calculate constrained param
  params.s2de <- c(opt.s2de$par,opt.s2de$par[4]*(genomesize-opt.s2de$par[2]/opt.s2de$par[1]-opt.s2de$par[3]))
  print(params.s2de)
  pancore.s2de <- f.star(params.s2de,ng)
  rms.s2de <- f.rms(c(mean.pangenome,mean.core),c(pancore.s2de$pan,pancore.s2de$core))

  # Add the predicted curves to the plot
  lines(1:ng,pancore.s2de$pan,col='blue',lwd=2,lty=3)
  lines(1:ng,pancore.s2de$core,col='blue',lwd=2,lty=3)

  ## Predict the core and pangenome curves
  ## using model 2D+E on a fixed tree
  ## (calculated via the gene frequency spectrum)

  mymodel <- "fixed" #use the coalescent tree w/pancore
  myfitting <- "chi2" #fit it using the Chi^2
  constr <- 1 # G0 constrained to the mean genome size during fitting
  mymethod <- "Nelder" # alternative fitting routine: "BFGS","Nelder"

  # recursively optimize using initial parameters and dataset to fit
  opt.f2de <- f.recurse(opt.c2de$par,r.data=pancore)

  # get optimized parameters and calculate constrained param
  params.f2de <- c(opt.f2de$par,opt.f2de$par[4]*(genomesize-opt.f2de$par[2]/opt.f2de$par[1]-opt.f2de$par[3]))
  print(params.f2de)
  pancore.f2de <- f.meanpancore(f.fixed.spec(params.f2de,treetable))
  rms.f2de <- f.rms(c(mean.pangenome,mean.core),c(pancore.f2de$pan,pancore.f2de$core))

  lines(1:ng,pancore.f2de$pan,col='orange',lwd=2,lty=1)
  lines(1:ng,pancore.f2de$core,col='orange',lwd=2,lty=1)

  # print to PDF if interactive
  if(.Internal(dev.displaylist())) {
    dev.copy2pdf(file=paste("pancore-",taxaname,".pdf",sep=""),width=6,height=6)
  }

  ##----------------------------------------
  ## Table 2 -- RMS values for 2D+E fits
  ##----------------------------------------
  # add results to Table 2
  table2 <- rbind(table2,cbind(
    taxaname,
    rms.perm,
    rms.cde/rms.perm,
    rms.c2de/rms.perm,
    rms.s2de/rms.perm,
    rms.f2de/rms.perm))

  ##--------------------------------------------------------
  ## Table 3 -- coalescent 2D+E best fits and predictions
  ##--------------------------------------------------------

  aa <- theta1 <- params.c2de[2]
  bb <- rho1 <- params.c2de[1]
  cc <- theta2 <- params.c2de[5]
  dd <- rho2 <- params.c2de[4]

  fslow <- theta1/rho1
  ffast <- theta2/rho2
  Gess <- params.c2de[3]

  y <- 100
  Gnew100 <- 1/2/y * (sqrt(aa^2-2*aa*bb*y+2*aa*cc+2*aa*dd*y+bb^2*y^2+2*bb*cc*y-2*bb*dd*y^2+cc^2-2*cc*dd*y+dd^2*y^2)+aa-bb*y+cc-dd*y+2*y)
  y <- 1000
  Gnew1000<- 1/2/y * (sqrt(aa^2-2*aa*bb*y+2*aa*cc+2*aa*dd*y+bb^2*y^2+2*bb*cc*y-2*bb*dd*y^2+cc^2-2*cc*dd*y+dd^2*y^2)+aa-bb*y+cc-dd*y+2*y)
  Gcore100 <- f.coalescent(params.c2de,100)$core[100]
  Gcore1000 <- f.coalescent(params.c2de,1000)$core[1000]

  # add results to Table 3
  table3 <- rbind(table3,cbind(
    taxaname,
    theta1,
    rho1,
    theta2,
    rho2,
    Gess,
    fslow,
    ffast,
    Gnew100,
    Gnew1000,
    Gcore100,
    Gcore1000))

} #end allkeys

colnames(table1) <- c("taxaname","ng","Ngenes","G0","Ngenes/G0","Gcore","Gcore/G0","Gpan","Gpan/G0","dprot")
colnames(table2) <- c("taxaname","RMS(data)","Coalescent 1D+E","Coalescent 2D+E","Star 2D+E","Fixed 2D+E")
colnames(table3) <- c("taxaname","theta1","rho1","theta2","rho2","Gess","fslow","ffast","Gnew(100)","Gnew(1000)","Gcore(100)","Gcore(1000)")
write.table(table1,f="table1.csv",sep="\t")
write.table(table2,f="table2.csv",sep="\t")
write.table(table3,f="table3.csv",sep="\t")

print("done")
