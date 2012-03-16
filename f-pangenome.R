##---------------------------------------------------
f.getspectrum <- function(mat) {
#----------------------------------------------------
# Function to compute gene family frequency spectrum
# given a matrix of gene clusters
#----------------------------------------------------
  ng <- dim(mat)[2]
  clusters <- apply(mat,1,function(x) sum(x>0))
  Gk <- hist(clusters,breaks=0:(ng+1),right=F,plot=F)$counts[-1] #observed G(k)
  return(Gk)
}

##---------------------------------------------------
f.meanpancore <- function(Gk) {
#----------------------------------------------------
# Function to compute mean pan and core genome curves
# given the gene frequency spectrum
#----------------------------------------------------
  ng <- length(Gk)
  gk.pan <- (1:ng)*0 #mean pangenome curve from G(k)
  gk.core <- (1:ng)*0 #mean core curve from G(k)

  #a gene is present in k genomes out of ng
  #the probability that the gene is absent in n genomes out of ng is...
  Pabs <- ((1:ng) %o% (1:ng))*0
  #the probability that the gene is present in all n genomes out of ng is...
  Pall <- ((1:ng) %o% (1:ng))*0

  #calculate Pabs(n,k)
  for (k in 1:ng) {
    for (n in 1:ng) {
      if (n<=(ng-k)) {
	Pabs[n,k] <- prod((ng+1-k-1:n)/(ng+1-1:n))
      }
    }
  }

  #calculate Pall(n,k)
  for (n in 1:ng) { 
    for (k in n:ng) {
      Pprod <- 1:n
      for (m in 1:n) {
	Pprod[m] <- (k-m+1)/(ng-m+1)
      }
      Pall[n,k] <- prod(Pprod)
    }
  }

  #calculate gk.pan and gk.core
  for (n in 1:ng) {
      gk.pan[n] <- sum(Gk * (1-Pabs[n,]))
      gk.core[n] <- sum(Gk[n:ng] * Pall[n,n:ng])
  }

  return(list("pan"=gk.pan,"core"=gk.core))
}

##---------------------------------------------------
f.pangenome <- function(mat,reps) {
#----------------------------------------------------
# function to compute pangenome curve permutations
# given 1) a matrix of gene clusters
#       2) the number of permutations 
#----------------------------------------------------

  ng <- dim(mat)[2] #number of genomes
  mat.ret <- NULL #matrix to return
  mat.pang <- (1:ng)*0 # initialize permutation pangenome
  mat.tmp <- mat[rowSums(mat)>0,]>0 #use only clusters that have data, presence/absence
  for (i in 1:reps) {
    cat(i," ") #progress
    samples <- sample(ng,ng,replace=F)
    seen <- mat.tmp[,samples[1]]
    mat.pang[1] <- sum(seen)
    for (j in 2:ng) {
      mat.pang[j] <- sum(mat.tmp[,samples[j]] - seen >0)
      seen <- mat.tmp[,samples[j]] | seen
    }
    mat.ret <- cbind(mat.ret,cumsum(mat.pang))
  }
  cat("\n")
  return(mat.ret)
}

##---------------------------------------------------
f.core <- function(mat,reps) {
#----------------------------------------------------
# function to compute core genome curve permutations
# given 1) a matrix of gene clusters
#       2) the number of permutations 
#----------------------------------------------------

  ng <- dim(mat)[2] #number of genomes
  mat.ret <- NULL #matrix to return
  mat.core <- (1:ng)*0 #initialize
  mat.save <- (1:ng)*0 #initialize
  mat.tmp <- 1*(mat[rowSums(mat)>0,]>0) #use only clusters that have data, presence/absence
  for (i in 1:reps) {
    cat(i, " ")
    samples <- sample(ng,ng,replace=F)
    seen <- mat.tmp[,samples]
    mat.core <- which(seen[,1]==1)
    mat.save[1] <- length(mat.core)
    for (j in 2:ng) {
	mat.core <- intersect(which(seen[,j]==1),mat.core)
	mat.save[j] <- length(mat.core)
    }
    mat.ret <- cbind(mat.ret,mat.save)
  }
  cat("\n")
  return(mat.ret)
}

##---------------------------------------------------
f.coalescent <- function(x,ng) {
#----------------------------------------------------
# function to calculate stationary distribution of 
# the IMG model on a coalescent
# given A1) deletion rate class 1 (rho1)
#       A2) insertion rate class 1 (theta1)
#       A3) number of essential genes (gess)
#       A4) deletion rate class 2 (rho2)
#       A5) insertion rate class 2 (theta2)
#	B) number of genomes
# 1D (rho1,theta1), 1D+E (rho1,theta1,gess), 2D (rho1,theta1,0,rho2,theta2), 2D+E (rho1,theta1,gess,rho2,theta2)
#----------------------------------------------------

  rho1 <- x[1] 
  theta1 <- x[2]

  if (is.na(x[3])) {gess<-0} else {gess <- x[3]}
  if (is.na(x[4])) {theta2<-0;rho2<-1} else {rho2 <- x[4];theta2 <- x[5]}

  core <- (1:ng)*0
  pangenome <- (1:ng)*0

  for (k in 1:ng) {
  #size of pan-genome after sampling K genomes from N
    pangenome[k] <- theta1 * sum(1/(rho1-1+1:k)) + theta2 * sum(1/(rho2-1+1:k)) + gess

    specprod1 <- (k-1:k+1)/(k-1:k+rho1)
    specprod2 <- (k-1:k+1)/(k-1:k+rho2)

    core[k] <- (theta1/k)*prod(specprod1[1:k]) + (theta2/k)*prod(specprod2[1:k]) + gess
  }

  return(list("pan"=pangenome,"core"=core))
}


##---------------------------------------------------
f.coalescent.spec <- function(x,ng) {
#----------------------------------------------------
# function to calculate gene family frequency spectrum
# for the IMG model on a coalescent
# given A1) deletion rate class 1 (rho1)
#       A2) insertion rate class 1 (theta1)
#       A3) number of essential genes (gess)
#       A4) deletion rate class 2 (rho2)
#       A5) insertion rate class 2 (theta2)
#	B) number of genomes
# 1D (rho1,theta1), 1D+E (rho1,theta1,gess), 2D (rho1,theta1,0,rho2,theta2), 2D+E (rho1,theta1,gess,rho2,theta2)
#----------------------------------------------------

  rho1 <- x[1]
  theta1 <- x[2]

  if (is.na(x[3])) {gess<-0} else {gess <- x[3]}
  if (is.na(x[4])) {theta2<-0;rho2<-1} else {rho2 <- x[4];theta2 <- x[5]}

  spec <- (1:ng)*0

  specprod1 <- (ng-1:ng+1)/(ng-1:ng+rho1)
  if (theta2 > 0) {
    specprod2 <- (ng-1:ng+1)/(ng-1:ng+rho2)
  } else {
    specprod2 <- 0*1:ng
  }

  for (k in 1:ng) {
  spec[k] <- (theta1/k)*prod(specprod1[1:k]) + (theta2/k)*prod(specprod2[1:k])
  }

  spec[ng] <- spec[ng] + gess

  return(spec)
}

##---------------------------------------------------
f.star <- function(x,ng) {
#----------------------------------------------------
# function to calculate stationary distribution of 
# the IMG model on a star tree
# given A1) deletion rate class 1 (vt1)
#       A2) insertion rate class 1 (ut1)
#       A3) number of essential genes (gess)
#       A4) deletion rate class 2 (vt2)
#       A5) insertion rate class 2 (ut2)
#	B) number of genomes
# 1D (vt1,ut1), 1D+E (vt1,ut1,gess), 2D (vt1,ut1,0,vt2,ut2), 2D+E (vt1,ut1,gess,vt2,ut2)
#----------------------------------------------------

  vt1 <- x[1] #v*t
  ut1 <- x[2]

  if (is.na(x[3])) {gess<-0} else {gess<- x[3]}
  if (is.na(x[4])) {ut2<-0;vt2<-1} else {vt2 <- x[4];ut2 <- x[5]}

  core <- (1:ng)*0
  pangenome <- (1:ng)*0
  core[1] <- ut1/vt1 + ut2/vt2 + gess
  pangenome[1] <- ut1/vt1 + ut2/vt2 + gess

  for (k in 2:ng) {
    core[k] <- (ut1/vt1)*exp(-k*vt1) + (ut2/vt2)*exp(-k*vt2) + gess
    pangenome[k] <- (ut1/vt1)*(1 + k*(1 - exp(-vt1)) - (1-exp(-vt1))^k) + (ut2/vt2)*(1 + k*(1 - exp(-vt2)) - (1-exp(-vt2))^k) + gess
  }

  return(list("pan"=pangenome,"core"=core))
}

##---------------------------------------------------
f.getdesc <- function(tree,node) {
#----------------------------------------------------
# Function to get the number of descendent nodes
# from a tree in treetable format
#----------------------------------------------------

  nodes <- node
  diffnodes <- 1
  while(diffnodes>0) {
    oldnodes <- length(nodes)
    nodes <- c(nodes, tree$desc_a[nodes],tree$desc_b[nodes])
    nodes <- sort(unique(nodes))
    diffnodes <- length(nodes) - oldnodes
  }
  nodes <- nodes[-which(nodes==0)]
  nodes <- nodes[-which(nodes==node)]
  return(length(nodes))
}

##---------------------------------------------------
f.dprot <- function(tree,ng) {
#----------------------------------------------------
# Function to get the average distance from the ancestor
# to each node in a tree in treetable format
#----------------------------------------------------

  dists <- 1:ng*0;
  for (leaf in 1:ng) {
      dists[leaf] <- dists[leaf] + tree[leaf,'dist']
      ancestor <- which(tree[,1:2]==leaf,arr.ind=TRUE)[1]
      while(!is.na(ancestor)) {
	dists[leaf] <- dists[leaf] + tree[ancestor,'dist']
	ancestor <- which(tree[,1:2]==ancestor,arr.ind=TRUE)[1]
      }
  }
  return(mean(dists))
}

##---------------------------------------------------
f.rms <- function(expected,observed) {
#----------------------------------------------------
# Function to calculate root-mean-square value
#----------------------------------------------------

  if(is.null(dim(observed))) {
      sqrt( sum( (expected - observed)^2 ) /length(expected) )
  } else {
      sqrt( colSums( (expected - observed)^2 ) /length(expected) )
  }
}

##---------------------------------------------------
f.fixed.spec <- function(x,treetable) {
#----------------------------------------------------
# Function to compute unconstrained gene frequency spectrum on a fixed tree
# given x[1] deletion rate class 1 (vt1)
#       x[2] insertion rate class 1 (ut1)
#       x[3] number of essential genes (gess)
#       x[4] deletion rate class 2 (vt2)
#       x[5] insertion rate class 2 (ut2)
#	treetable table of distances and daughters
# 1D (vt1,ut1), 1D+E (vt1,ut1,gess), 2D (vt1,ut1,0,vt2,ut2), 2D+E (vt1,ut1,gess,vt2,ut2)
#----------------------------------------------------

  v1 <- x[1]
  u1 <- x[2]
  if (is.na(x[3])) {gess<-0} else {gess<-x[3]}
  if (is.na(x[4])) {v2<-1;u2<-0} else {v2<-x[4];u2<-x[5]}

  losegain <- cbind(lose(v1,treetable$dist),lose(v2,treetable$dist),retain(v1,treetable$dist),retain(v2,treetable$dist))
  nn <- nrow(treetable) #number of nodes
  ng <- (nn+1)/2 #number of genomes if bifurcating tree
  nd <- NA*1:nn #number of descendent nodes not counting itself
  ndg <- NA*1:nn #number of descendent genomes
  gexp1 <- 0*1:nn #expected number of genes
  gexp2 <- 0*1:nn #expected number of genes
  gprob1 <- array(0,dim=c(nn,ng+1)) #probability of k genomes from k=0
  gprob2 <- array(0,dim=c(nn,ng+1)) #probability of k genomes from k=0

  for (a in 1:nn) {
    # Let a be a node in this fixed tree
    # let treetable$dist[a] be the length of the branch leading to node a
    # let ndg[a] be the number of genomes in the data that descend from node a
    nd[a] <- f.getdesc(treetable,a)
    ndg[a] <- (nd[a]+2)/2 #number of descendent genomes (tips)

    # For a single class of dispensable gene families in the IMG model,
    # gexp1[a] is the expected number of families present in 'a' that 
    # arose on the branch leading to 'a'
    # For a second class of dispensable gene families in the IMG model,
    # gexp2[a] is the expected number of families present in 'a' that 
    # arose on the branch leading to 'a'

    #if root else
    if (ndg[a] == ng) {
      gexp1[a] <- u1/v1
      gexp2[a] <- u2/v2
    } else {
      gexp1[a] <- (u1/v1) * losegain[a,1]
      gexp2[a] <- (u2/v2) * losegain[a,2]
    }

    # gprob[a,k+1] is the probability that a family present at 'a'
    # is present in k out of ndg[a] genomes that descend from 'a'

    # if tip else
    k=0
    if (ndg[a] == 1) {
      # If 'a' is a tip node, then gprob[a,2] = 1
      # and gprob[a,k+1] = 0 for k neq 1 (initial array of zeros)
      gprob1[a,k+1+1] <- 1
      gprob2[a,k+1+1] <- 1
    } else {
      # If 'a' is an internal node, gprob[a,k+1] = 0 for k > ndg[a]
      # calculate gprob[a,k+1] for 0 <= k <= ndg[a]
	desc_a <- treetable$desc_a[a]
	desc_b <- treetable$desc_b[a]
# 	t_a <- treetable$dist[desc_a]
# 	t_b <- treetable$dist[desc_b]

	#for k = 0
	gprob1[a,k+1] <- (losegain[desc_a,1] + losegain[desc_a,3] * gprob1[desc_a,k+1]) *
			 (losegain[desc_b,1] + losegain[desc_b,3] * gprob1[desc_b,k+1])
	gprob2[a,k+1] <- (losegain[desc_a,2] + losegain[desc_a,4] * gprob2[desc_a,k+1]) *
			 (losegain[desc_b,2] + losegain[desc_b,4] * gprob2[desc_b,k+1])

	#for 1 <= k <= ndg[a]
	for (k in 1:ndg[a]) {
	  gprob1[a,k+1] <- 
	    losegain[desc_a,3]*losegain[desc_b,1]*gprob1[desc_a,k+1] +
	    losegain[desc_b,3]*losegain[desc_a,1]*gprob1[desc_b,k+1] +
	    losegain[desc_a,3]*losegain[desc_b,3]*
	    sum(unlist(sapply(0:k,function(j) gprob1[desc_a,j+1]*gprob1[desc_b,k-j+1])))

	  gprob2[a,k+1] <- 
	    losegain[desc_a,4]*losegain[desc_b,2]*gprob2[desc_a,k+1] +
	    losegain[desc_b,4]*losegain[desc_a,2]*gprob2[desc_b,k+1] +
	    losegain[desc_a,4]*losegain[desc_b,4]*
	    sum(unlist(sapply(0:k,function(j) gprob2[desc_a,j+1]*gprob2[desc_b,k-j+1])))
	}
    } #end p_a(k)
  } #end for a

  Gk1 <- unlist(sapply(1:ng,function(k) {
    sum(unlist(sapply(1:nn,function(a) {
      gexp1[a]*gprob1[a,k+1]
      })))
  }))

  Gk2 <- unlist(sapply(1:ng,function(k) {
    sum(unlist(sapply(1:nn,function(a) {
      gexp2[a]*gprob2[a,k+1]
      })))
  }))

  Gk <- Gk1 + Gk2
  Gk[ng] <- Gk[ng]+gess

  return(Gk)
}

retain <- function(v,t) {
# Let the probability that a family is retained for a time $t$ be $r(t) = e^{-v t}$
  exp(-v*t)
}
lose <- function(v,t) {
# let the probability that it is lost during time $t$ be $l(t) = 1 - e^{-v t}$
  1 - exp(-v*t)
}


##---------------------------------------------------
f.fit <- function(x,data,constr,modeltype=mymodel,fitting=myfitting,genomesize=NA,treetable=NA,ng=NA) {
#--------------------------------------
# Fitting function
# using Chi^2 or Sum of Squares
# using G(k) or core and pangenome curve functions
# on models 1D, 1D+E, 2D, or 2D+E
# param fitting=chi2 or sumsq
# param tree=coalescent, star, or fixed
# param constr=1 (constrained) or 0 (unconstrained) to G0
#--------------------------------------

  if (any(x<0)) {return(NA)}
  if (constr==1) {
    if (length(x)==1) { # 1D constrained
      v1 <- x[1]
      u1 <- v1*genomesize
      v2 <- 255
      gess <- 0
      u2 <- 0
    } else if (length(x)==2) { # 1D+E constrained
      v1 <- x[1]
      gess <- x[2]
      u1 <- v1*(genomesize-gess)
      v2 <- 255
      u2 <- 0
    } else if (length(x)==3) { # 2D constrained
      v1 <- x[1]
      u1 <- x[2]
      v2 <- x[3]
      gess <- 0
      u2 <- v2*(genomesize-u1/v1)
    } else if (length(x)==4) { # 2D+E constrained
      v1 <- x[1]
      u1 <- x[2]
      gess <- x[3]
      v2 <- x[4]
      u2 <- v2*(genomesize-u1/v1-gess)
    } else {return(NA)}
  } else if (constr==0) {
    if (length(x)==2) { # 1D unconstrained
      v1 <- x[1]
      u1 <- x[2]
      v2 <- 255
      gess <- 0
      u2 <- 0
    } else if (length(x)==3) { # 1D+E unconstrained
      v1 <- x[1]
      u1 <- x[2]
      gess <- x[3]
      v2 <- 255
      u2 <- 0
    } else if (length(x)==4) { # 2D unconstrained
      v1 <- x[1]
      u1 <- x[2]
      gess <- 0
      v2 <- x[3]
      u2 <- x[4]
    } else if (length(x)==5) { # 2D+E unconstrained
      v1 <- x[1]
      u1 <- x[2]
      gess <- x[3]
      v2 <- x[4]
      u2 <- x[5]
    } else {return(NA)}
  }

  pars <- c(v1,u1,gess,v2,u2) #rho1,theta1,gess,rho2,theta2
#   print(pars)
  if (v1>v2 || any(pars<0) ) {return(NA)}

  if (modeltype=="fixed.spec") {
    theory <- f.fixed.spec(pars,treetable)
  } else if (modeltype=="fixed") {
    result <- f.meanpancore(f.fixed.spec(pars,treetable))
    theory <- c(result$pan,result$core)
  } else if (modeltype=="coalescent.spec") {
    theory <- f.coalescent.spec(pars,ng)
  } else if (modeltype=="coalescent") {
    result <- f.coalescent(pars,ng)
    theory <- c(result$pan,result$core)
  } else if (modeltype=="star") {
    result <- f.star(pars,ng)
    theory <- c(result$pan,result$core)
  } else {return(NA)}

  if (fitting=="sumsq") { value <- sum((theory-data)^2)
  } else if (fitting=="rms") { value <- sqrt( sum( (theory - data)^2 ) /length(theory) )

  } else if (fitting=="chi2") { value <- sum((theory-data)^2/theory)
  } else {return(NA)}
  # lines(theory,col='grey')
  return(value)
}

##-------------------------------------------------
f.recurse <- function(pinitial,r.data,r.constr=constr,r.modeltype=mymodel,r.fitting=myfitting,r.genomesize=genomesize,r.treetable=treetable,r.ng=ng,r.method=mymethod,r.maxit=10000,r.reltol=1e-6) {
#--------------------------------------
# Recursive Fitting metafunction
# takes initial parameters and data to fit
# uses best fit from previous iteration
# to start new iteration up to 10x
#--------------------------------------
  saved <- NULL
  while(length(saved$value)<10) {
    if(is.null(saved)) { params <- pinitial
    } else { params <- opt$par #recursive
    }
    opt <- optim(
      params,
      function(x) f.fit(x,data=r.data,constr=r.constr,modeltype=r.modeltype,fitting=r.fitting,treetable=r.treetable,genomesize=r.genomesize,ng=r.ng),
      method=r.method,
      control=list(trace=1,parscale=params+1e-6,reltol=r.reltol,maxit=r.maxit)
    )
    if(opt$convergence>0) {warning("did not converge!")}

    saved <- rbind(saved,as.data.frame(t(unlist(opt))))
    print(saved$value)

    if (dim(saved)[1] > 1) {
      if (identical(all.equal(sort(saved$value)[2],sort(saved$value)[1]),TRUE)) {
	break;
      }
    }
  }

return(opt)
}
