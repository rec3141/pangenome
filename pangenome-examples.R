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
## USAGE:
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
datadir <- "./data/"
# datadir <- ""

# read in matrix of all Bacilli gene clusters
# less strict clustering (default):
# infile <- paste(datadir,"Bacilli-clusters.1e-10.6.csv",sep="")
# more strict clustering:
infile <- paste(datadir,"Bacilli-clusters.1e-30.7.csv",sep="")
# infile <- "mbgd.Bacilli.csv"

mat.all <- read.csv(infile,sep="\t")
#mat.read <- mat.all[,3+1:(dim(mat.all)[2]-7)] #remove annotations

# initialize tables
table1 <- NULL
table2 <- NULL
table3 <- NULL
table4 <- NULL

# get list of all key files to process
# allkeys <- list.files(path="./",pattern="key-.+")
allkeys <- c("key-Staph_aureus", "key-Strep_pyogenes", "key-Strep_pneumoniae", "key-B_cereus", "key-Listeria", "key-Staphylococcus", "key-Streptococcus", "key-Bacillaceae","key-Lactobacillales", "key-Bacilli")

for (mykey in allkeys) {
  # the key matches the names to RefSeq IDs
  keys <- read.table(paste(datadir,mykey,sep=""),sep="\t",row.names=1)
  # the taxaname is the name of the taxonomic group
  taxaname <- unlist(strsplit(mykey,'-'))[2]
  print(taxaname)
  # get the columns which are of interest
  mycols <- unique(unlist(lapply(keys[,2],function(x) grep(x,colnames(mat.all)))))
  # get a matrix of clusters containing gene families in the genomes of interest
  mat <- mat.all[rowSums(mat.all[,mycols]>0)>0,mycols]

  genomesize <- mean(colSums(mat>0)) # mean genome size measured in gene families
  ng <- dim(mat)[2] #number of genomes

  # load phylogenetic tree table for this taxonomic group
  treename <- paste("treetable",taxaname,sep="-")
  treetable <- read.table(paste(datadir,treename,sep=""),sep="\t",row.names=1)
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

  # Set up plot with 2 rows, 1 column
  par(mfrow=c(2,1))

  #calculate gene family frequency spectrum
  Gk <- f.getspectrum(mat)

  #plot gene family frequency spectrum
  bp.x <- barplot(Gk, ylab="Frequency of Gene Families",xlab="Number of Genomes",names.arg=1:ng,main=paste(taxaname,"Gene Family Frequency Spectrum",sep="\n"),log='y')
  legend("top",c("1D+E on coalescent","2D+E on coalescent","2D+E on fixed tree"),lty=c(2,1,1),col=c('red','red','blue'),lwd=2)

  ## Predict the gene family frequency spectrum
  ## using model 1D+E on a coalescent

  mymaxit <- 10000
  myreltol <- 1e-6

  mymodel <- "coalescent.spec" #use the coalescent tree w/G(k)
  myfitting <- "chi2" #fit it using this error measurement
  constr <- 1 # G0 constrained to the mean genome size during fitting
  # !! if constr changes, need to change input parameters as well
  mymethod <- "Nelder" # alternative fitting routine: "BFGS"

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
  myfitting <- "chi2" #fit it using this error function
  constr <- 1 # G0 constrained to the mean genome size during fitting
  # !! if constr changes, need to change input parameters as well
  mymethod <- "Nelder" # alternative fitting routine: "BFGS"

  #set initial parameters and recursively optimize
  pinitial.spec.c2de <- c(params.spec.cde[1]/10,params.spec.cde[2]/100,params.spec.cde[3],params.spec.cde[1]*10000)
  opt.spec.c2de <- f.recurse(pinitial.spec.c2de,r.data=Gk)

  # get optimized parameters, calculate constrained parameter
  params.spec.c2de <- c(opt.spec.c2de$par,opt.spec.c2de$par[4]*(genomesize-opt.spec.c2de$par[2]/opt.spec.c2de$par[1]-opt.spec.c2de$par[3]))
  print(params.spec.c2de)
  spec.c2de <- f.coalescent.spec(params.spec.c2de,ng)
  lines(bp.x,spec.c2de,col='red',lwd=2)

  ## Predict the gene family frequency spectrum
  ## using model 2D+E on a fixed tree

  mymodel <- "fixed.spec" #use a fixed tree w/G(k)
  myfitting <- "chi2" #fit it using this error function
  constr <- 1 # G0 constrained to the mean genome size during fitting;
  # !! if constr changes, need to change input parameters as well
  mymethod <- "Nelder" # alternative fitting routines: "BFGS","Nelder"

  #set initial parameters and recursively optimize
#   pinitial.spec.f2de <- opt.spec.c2de$par/c(1,1,1,1)
#   mymaxit <- 200
#   myreltol <- 1e-4
#   opt.spec.f2de <- f.recurse(pinitial.spec.f2de,r.data=Gk)

  pinitial.spec.f2de <- opt.spec.c2de$par
  mymaxit <- 10000
  myreltol <- 1e-6
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

  ##----------------------------------------
  ## Table 4 -- Chi^2 values for Gk fits
  ##----------------------------------------

  # add results to Table 4
  table4 <- rbind(table4,cbind(
    taxaname,
    opt.spec.cde$value,
    opt.spec.c2de$value,
    opt.spec.f2de$value
    ))

  ##-----------------------------------------------------
  ## Fit to the Pangenome and core genome curves
  ##-----------------------------------------------------

  # Calculate 100 permutations each of the pangenome and core genome
  perm.pangenome <- f.pangenome(mat,500)
  perm.core <- f.core(mat,500)

  # Calculate the exact mean pan and core genome curves
  # from the gene frequency spectrum G(k)
  mean.pangenome <- f.meanpancore(Gk)$pan
  mean.core <- f.meanpancore(Gk)$core
  pancore <- c(mean.pangenome,mean.core)

  # Calculate the RMS value for the permutations
  rms.perm <- mean(f.rms(c(mean.pangenome,mean.core),rbind(perm.pangenome,perm.core)))

  # Prepare a new plot window
  plot(1:ng,xlim=c(1,ng),ylim=c(0.9*min(mean.core), 1.1*max(mean.pangenome)),log='y',xlab="Genomes added", ylab="Clusters of Gene Families",main=paste(taxaname,"Pangenome and Core genome",sep="\n"),pch='')
  legend("right",c("1D+E on coalescent","2D+E on coalescent","2D+E on star tree","2D+E on a fixed tree"),lty=c(2,1,3,1),col=c('red','red','grey','blue'),lwd=2)

  # Plot polygons outlining permutations
  polygon(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))), col="gray88",border=NA)
  polygon(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.pangenome,1,max))), col="gray88",border=NA)

  # Add the mean pan and core genome curves to the plot
  points(1:ng,mean.pangenome)
  points(1:ng,mean.core)

  ## Predict the core and pangenome curves
  ## using model 1D+E on a coalescent tree

  mymodel <- "coalescent" #use the coalescent tree w/pancore
  myfitting <- "rms" #fit it using the Chi^2
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
  lines(1:ng,pancore.cde$pan,col='red',lty=2,lwd=2)
  lines(1:ng,pancore.cde$core,col='red',lty=2,lwd=2)

  ## Predict the core and pangenome curves
  ## using model 2D+E on a coalescent tree

  mymodel <- "coalescent" #use the coalescent tree w/G(k)
  myfitting <- "rms" #fit it using the Chi^2
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
  myfitting <- "rms" #fit it using the Chi^2
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
  myfitting <- "rms" #fit it using the Chi^2
  constr <- 1 # G0 constrained to the mean genome size during fitting
  mymethod <- "Nelder" # alternative fitting routine: "BFGS","Nelder"

  # recursively optimize using initial parameters and dataset to fit
  pinitial.f2de <- opt.spec.c2de$par/c(1,1,1,1)
  mymaxit <- 10000
  myreltol <- 1e-6
  opt.f2de <- f.recurse(pinitial.f2de,r.data=pancore,r.maxit=1000,r.reltol=1e-8)

  # get optimized parameters and calculate constrained param
  params.f2de <- c(opt.f2de$par,opt.f2de$par[4]*(genomesize-opt.f2de$par[2]/opt.f2de$par[1]-opt.f2de$par[3]))
  print(params.f2de)
  pancore.f2de <- f.meanpancore(f.fixed.spec(params.f2de,treetable))
  rms.f2de <- f.rms(c(mean.pangenome,mean.core),c(pancore.f2de$pan,pancore.f2de$core))

  lines(1:ng,pancore.f2de$pan,col='blue',lwd=2,lty=1)
  lines(1:ng,pancore.f2de$core,col='blue',lwd=2,lty=1)

  # print to PDF if interactive
  if(.Internal(dev.displaylist())) {
    dev.copy2pdf(file=paste("pancore-",taxaname,".pdf",sep=""),width=6,height=6)
  }


  ##----------------------------------------
  ## Plot Figure 3
  ##----------------------------------------

  par(mfrow=c(2,1)) #2 rows, 1 column

  #pangenome
  plot(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.core,1,max))),pch='',xlab=paste(taxaname,"genomes, (n)",sep=" "), ylab="Clusters of Gene Families, Gcore(n)",main="Pangenome")
  polygon(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.pangenome,1,max))), col="gray88",border=NA)
  points(1:ng,mean.pangenome)
  lines(1:ng,pancore.cde$pan,col='red',lty=2,lwd=2)
  lines(1:ng,pancore.c2de$pan,col='red',lwd=2)

  #core genome
  plot(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))),pch='',xlab=paste(taxaname,"genomes, (n)",sep=" "), ylab="Clusters of Gene Families, Gpan(n)",main="Core genome")
  polygon(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))), col="gray88",border=NA)
  points(1:ng,mean.core)
  lines(1:ng,pancore.cde$core,col='red',lty=2,lwd=2)
  lines(1:ng,pancore.c2de$core,col='red',lwd=2)

  # print to EPS if interactive, for import to SVG editor
  if(.Internal(dev.displaylist())) {
    dev.copy2eps(file=paste("fig3-",taxaname,".eps",sep=""),width=6,height=12)
  }

  ##----------------------------------------
  ## Plot Figure 4
  ##----------------------------------------

  par(mfrow=c(2,1)) #2 rows, 1 column

  #pangenome
  plot(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.core,1,max))),pch='',xlab=paste(taxaname,"genomes, (n)",sep=" "), ylab="Clusters of Gene Families, Gcore(n)",main="Pangenome")
  polygon(c(1:ng, rev(1:ng)), c(apply(perm.pangenome,1,min), rev(apply(perm.pangenome,1,max))), col="gray88",border=NA)
  points(1:ng,mean.pangenome)
  lines(1:ng,pancore.c2de$pan,col='red',lwd=2)
  lines(1:ng,pancore.s2de$pan,col='black',lwd=2)
  lines(1:ng,pancore.f2de$pan,col='blue',lwd=2)

  #core genome 
  plot(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))),pch='',xlab=paste(taxaname,"genomes, (n)",sep=" "), ylab="Clusters of Gene Families, Gpan(n)",main="Core genome")
  polygon(c(1:ng, rev(1:ng)), c(apply(perm.core,1,min), rev(apply(perm.core,1,max))), col="gray88",border=NA)
  points(1:ng,mean.core)
  lines(1:ng,pancore.c2de$core,col='red',lwd=2)
  lines(1:ng,pancore.s2de$core,col='black',lwd=2)
  lines(1:ng,pancore.f2de$core,col='blue',lwd=2)

  # print to EPS if interactive, for import to SVG editor  
  if(.Internal(dev.displaylist())) {
    dev.copy2eps(file=paste("fig4-",taxaname,".eps",sep=""),width=6,height=12)
  }

  ##----------------------------------------
  ## Plot Figure 5
  ##----------------------------------------

  # Add the G(k) plot window and best fit lines
  bp.x <- barplot(Gk, ylim=c(10,10000),ylab="Gene family frequency, G(k)",xlab=paste(taxaname,"genomes, k",sep=" "),names.arg=1:ng,main="Gene Family Frequency Spectrum",log='y',col="gray88",border="gray88")
  lines(bp.x,spec.cde,col='red',lty=3,lwd=2)
  lines(bp.x,spec.c2de,col='red',lty=2,lwd=2)
  lines(bp.x,spec.f2de,col='blue',lty=2,lwd=2)
  lines(bp.x,f.coalescent.spec(params.c2de,ng),col='red',lty=1,lwd=2)
  lines(bp.x,f.fixed.spec(params.f2de,treetable),col='blue',lty=1,lwd=2)
  lines(bp.x,f.fixed.spec(params.spec.f2de,treetable),col='grey',lty=2,lwd=2)

  # print to EPS if interactive, for import to SVG editor  
  if(.Internal(dev.displaylist())) {
    dev.copy2eps(file=paste("fig5-",taxaname,".eps",sep=""),width=6,height=12)
  }

  ##----------------------------------------
  ## Table 2 -- RMS values for pancore fits
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
  Gess <- params.spec.c2de[3]

  theta1 <- params.spec.c2de[2]
  rho1 <- params.spec.c2de[1]
  theta2 <- params.spec.c2de[5]
  rho2 <- params.spec.c2de[4]

  fslow <- theta1/rho1/G0
  ffast <- theta2/rho2/G0
  fess <- params.spec.c2de[3]/G0

  Gnew100 <- f.coalescent(params.spec.c2de,100)$pan[100]-f.coalescent(params.spec.c2de,100)$pan[99]
  Gnew1000 <- f.coalescent(params.spec.c2de,1000)$pan[1000]-f.coalescent(params.spec.c2de,1000)$pan[999]
  Gcore100 <- f.coalescent(params.spec.c2de,100)$core[100]
  Gcore1000 <- f.coalescent(params.spec.c2de,1000)$core[1000]
  Gpan100 <- f.coalescent(params.spec.c2de,100)$pan[100]
  Gpan1000 <- f.coalescent(params.spec.c2de,1000)$pan[1000]
  
  # add results to Table 3
  table3 <- rbind(table3,cbind(
    taxaname,
    Gess,
    theta1/rho1,
    rho1,
    theta2/rho2,
    rho2,
    fess,
    fslow,
    ffast,
    Gnew100,
    Gnew1000,
    Gpan100,
    Gpan1000,
    Gcore100,
    Gcore1000
))


colnames(table1) <- c("taxaname","ng","Ngenes","G0","Ngenes/G0","Gcore","Gcore/G0","Gpan","Gpan/G0","dprot")
colnames(table2) <- c("taxaname","RMS(data)","Coalescent 1D+E","Coalescent 2D+E","Star 2D+E","Fixed 2D+E")
colnames(table3) <- c("taxaname","Gess","theta1/rho1","rho1","theta2/rho2","rho2","fess","fslow","ffast","Gnew(100)","Gnew(1000)","Gcore(100)","Gcore(1000)","Gpan(100)","Gpan(1000)")
colnames(table4) <- c("taxaname","Coalescent 1D+E","Coalescent 2D+E","Fixed 2D+E")
write.table(table1,f="table1_new3.csv",sep="\t")
write.table(table2,f="table2_new3.csv",sep="\t")
write.table(table3,f="table3_new3.csv",sep="\t")
write.table(table4,f="table4_new3.csv",sep="\t")

} #end allkeys


print("done")
