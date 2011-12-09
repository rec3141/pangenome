#----------------------------------------------------
# Function to compute gene frequency spectrum
# given a matrix of gene clusters
#----------------------------------------------------
f.getspectrum <- function(mat) {
ng <- dim(mat)[2]
clusters <- apply(mat,1,function(x) sum(x>0))
Gk <- hist(clusters,breaks=0:(ng+1),right=F,plot=F)$counts[-1] #observed G(k)
return(Gk)
}

#----------------------------------------------------
# Function to compute mean pan and core genome curves
# given the gene frequency spectrum
#----------------------------------------------------
f.meanpancore <- function(Gk) {
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

#----------------------------------------------------
# function to compute pangenome curve permutations
# given 1) a matrix of gene clusters
#       2) the number of permutations 
#----------------------------------------------------
f.pangenome <- function(mat,reps) {
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

#----------------------------------------------------
# function to compute core genome curve permutations
# given 1) a matrix of gene clusters
#       2) the number of permutations 
#----------------------------------------------------
f.core <- function(mat,reps) {
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

#----------------------------------------------------
# function to calculate stationary distribution of 
# the IMG model on a coalescent
# given A1) deletion rate class 1 (rho1)
#       A2) insertion rate class 1 (theta1)
#       A3) number of essential genes
#       A4) deletion rate class 2 (rho2)
#       A5) insertion rate class 2 (theta2)
#	B) number of genomes
# 1D (rho1,theta1), 1D+E (rho1,theta1,ness), 2D (rho1,theta1,0,rho2,theta2), 2D+E (rho1,theta1,ness,rho2,theta2)
#----------------------------------------------------
f.coalescent <- function(x,ng) {

rho1 <- x[1] 
theta1 <- x[2]

if (is.na(x[3])) {ness<-0} else {ness <- x[3]}
if (is.na(x[4])) {theta2<-0;rho2<-1} else {rho2 <- x[4];theta2 <- x[5]}

core <- (1:ng)*0
pangenome <- (1:ng)*0

for (k in 1:ng) {
#size of pan-genome after sampling K genomes from N
  pangenome[k] <- theta1 * sum(1/(rho1-1+1:k)) + theta2*sum(1/(rho2-1+1:k)) + ness

  specprod1 <- (k-1:k+1)/(k-1:k+rho1)
  specprod2 <- (k-1:k+1)/(k-1:k+rho2)

  core[k] <- (theta1/k)*prod(specprod1[1:k]) + (theta2/k)*prod(specprod2[1:k]) + ness
}

return(list("pan"=pangenome,"core"=core))
}


#----------------------------------------------------
# function to calculate gene family frequency spectrum
# for the IMG model on a coalescent
# given A1) deletion rate class 1 (rho1)
#       A2) insertion rate class 1 (theta1)
#       A3) number of essential genes
#       A4) deletion rate class 2 (rho2)
#       A5) insertion rate class 2 (theta2)
#	B) number of genomes
# 1D (rho1,theta1), 1D+E (rho1,theta1,ness), 2D (rho1,theta1,0,rho2,theta2), 2D+E (rho1,theta1,ness,rho2,theta2)
#----------------------------------------------------
f.coalescent.spec <- function(x,ng) {
  rho1 <- x[1]
  theta1 <- x[2]

if (is.na(x[3])) {ness<-0} else {ness <- x[3]}
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

spec[ng] <- spec[ng] + ness

return(spec)
}

#----------------------------------------------------
# function to calculate stationary distribution of 
# the IMG model on a star tree
# given A1) deletion rate class 1 (vt1)
#       A2) insertion rate class 1 (ut1)
#       A3) number of essential genes
#       A4) deletion rate class 2 (vt2)
#       A5) insertion rate class 2 (ut2)
#	B) number of genomes
# 1D (vt1,ut1), 1D+E (vt1,ut1,ness), 2D (vt1,ut1,0,vt2,ut2), 2D+E (vt1,ut1,ness,vt2,ut2)
#----------------------------------------------------
f.star <- function(x,ng) {

vt1 <- x[1] #v*t
ut1 <- x[2]

if (is.na(x[3])) {ness<-0} else {ness<- x[3]}
if (is.na(x[4])) {ut2<-0;vt2<-1} else {vt2 <- x[4];ut2 <- x[5]}

core <- (1:ng)*0
pangenome <- (1:ng)*0
core[1] <- ut1/vt1 + ut2/vt2
pangenome[1] <- ut1/vt1 + ut2/vt2

for (k in 2:ng) {
  core[k] <- (ut1/vt1)*exp(-k*vt1) + (ut2/vt2)*exp(-k*vt2) + ness
  pangenome[k] <- (ut1/vt1)*(1 + k*(1 - exp(-vt1)) - (1-exp(-vt1))^k) + (ut2/vt2)*(1 + k*(1 - exp(-vt2)) - (1-exp(-vt2))^k) + ness
}

return(list("pan"=pangenome,"core"=core))
}

#----------------------------------------------------
# Function to get the number of descendent nodes
# from a tree in treetable format
#----------------------------------------------------
f.getdesc <- function(tree,node) {
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

#----------------------------------------------------
# Function to compute unconstrained gene frequency spectrum on a fixed tree
# given x[1] deletion rate class 1 (vt1)
#       x[2] insertion rate class 1 (ut1)
#       x[3] number of essential genes
#       x[4] deletion rate class 2 (vt2)
#       x[5] insertion rate class 2 (ut2)
#	treetable table of distances and daughters
# 1D (vt1,ut1), 1D+E (vt1,ut1,ness), 2D (vt1,ut1,0,vt2,ut2), 2D+E (vt1,ut1,ness,vt2,ut2)
#----------------------------------------------------
f.fixed.spec <- function(x,treetable) {
  v1 <- x[1]
  u1 <- x[2]
  if (is.na(x[3])) {ness<-0} else {ness<-x[3]}
  if (is.na(x[4])) {v2<-1;u2<-0} else {v2<-x[4];u2<-x[5]}

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
      gexp1[a] <- (u1/v1) * (1-exp(-v1*treetable$dist[a]))
      gexp2[a] <- (u2/v2) * (1-exp(-v2*treetable$dist[a]))
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
	t_a <- treetable$dist[desc_a]
	t_b <- treetable$dist[desc_b]

	#for k = 0
	gprob1[a,k+1] <- (lose(v1,t_a) + retain(v1,t_a) * gprob1[desc_a,k+1]) *
			 (lose(v1,t_b) + retain(v1,t_b) * gprob1[desc_b,k+1])
	gprob2[a,k+1] <- (lose(v2,t_a) + retain(v2,t_a) * gprob2[desc_a,k+1]) *
			 (lose(v2,t_b) + retain(v2,t_b) * gprob2[desc_b,k+1])

	#for 1 <= k <= ndg[a]
	for (k in 1:ndg[a]) {
	  gprob1[a,k+1] <- 
	    retain(v1,t_a)*lose(v1,t_b)*gprob1[desc_a,k+1] +
	    retain(v1,t_b)*lose(v1,t_a)*gprob1[desc_b,k+1] +
	    retain(v1,t_a)*retain(v1,t_b)*
	    sum(unlist(sapply(0:k,function(j) gprob1[desc_a,j+1]*gprob1[desc_b,k-j+1])))

	  gprob2[a,k+1] <- 
	    retain(v2,t_a)*lose(v2,t_b)*gprob2[desc_a,k+1] +
	    retain(v2,t_b)*lose(v2,t_a)*gprob2[desc_b,k+1] +
	    retain(v2,t_a)*retain(v2,t_b)*
	    sum(unlist(sapply(0:k,function(j) gprob2[desc_a,j+1]*gprob2[desc_b,k-j+1])))
	}
    } #end p_a(k)
  } #end for a

#can also get k=0 out of this?? number of extinct genes?
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
Gk[ng] <- Gk[ng]+ness

return(Gk)
}
# Let the probability that a family is retained for a time $t$ be $r(t) = e^{-v t}$
retain <- function(v,t) {exp(-v*t)}
# let the probability that it is lost during time $t$ be $l(t) = 1 - e^{-v t}$
lose <- function(v,t) {1 - exp(-v*t)}


#--------------------------------------
# Fitting function
# using Chi^2 or Sum of Squares
# using G(k) or core and pangenome curve functions
# on models 1D, 1D+E, 2D, or 2D+E
# param fitting=chi2 or sumsq
# param tree=coalescent, star, or fixed
# param constr=1 (constrained) or 0 (unconstrained) to G0
#--------------------------------------
f.fit <- function(x,constr,data,treetype,fitting,genomesize=NA,treetable=NA,ng=NA) {

if (constr==1) {
  if (length(x)==1) { # 1D constrained
    v1 <- x[1]
    u1 <- v1*genomesize
    v2 <- 255
    ness <- 0
    u2 <- 0
  } else if (length(x)==2) { # 1D+E constrained
    v1 <- x[1]
    u1 <- v1*genomesize
    v2 <- 255
    ness <- x[2]
    u2 <- 0
  } else if (length(x)==3) { # 2D constrained
    v1 <- x[1]
    u1 <- x[2]
    v2 <- x[3]
    ness <- 0
    u2 <- v2*(genomesize-u1/v1)
  } else if (length(x)==4) { # 2D+E constrained
    v1 <- x[1]
    u1 <- x[2]
    ness <- x[3]
    v2 <- x[4]
    u2 <- v2*(genomesize-u1/v1-ness)
  } else {return(NA)}
} else if (constr==0) {
  if (length(x)==2) { # 1D unconstrained
    v1 <- x[1]
    u1 <- x[2]
    v2 <- 255
    ness <- 0
    u2 <- 0
  } else if (length(x)==3) { # 1D+E unconstrained
    v1 <- x[1]
    u1 <- x[2]
    ness <- x[3]
    v2 <- 255
    u2 <- 0
  } else if (length(x)==4) { # 2D unconstrained
    v1 <- x[1]
    u1 <- x[2]
    ness <- 0
    v2 <- x[3]
    u2 <- x[4]
  } else if (length(x)==5) { # 2D+E unconstrained
    v1 <- x[1]
    u1 <- x[2]
    ness <- x[3]
    v2 <- x[4]
    u2 <- x[5]
  } else {return(NA)}
}

pars <- c(v1,u1,ness,v2,u2)

if (v1>v2 || any(pars<0) ) {return(6.022e23)}

if (treetype=="fixed.spec") {
  theory <- f.fixed.spec(pars,treetable)
} else if (treetype=="fixed") {
  result <- f.meanpancore(f.fixed.spec(pars,treetable))
  theory <- c(result$pan,result$core)
} else if (treetype=="coalescent.spec") {
  theory <- f.coalescent.spec(pars,ng)
} else if (treetype=="coalescent") {
  result <- f.coalescent(pars,ng)
  theory <- c(result$pan,result$core)
} else if (treetype=="star") {
  result <- f.star(pars,ng)
  theory <- c(result$pan,result$core)
} else {return(NA)}

if (fitting=="chi2") {value <- sum((theory-data)^2/data)}
else if (fitting=="sumsq") {value <- sum((theory-data)^2)}
else {return(NA)}
# lines(theory,col='grey')
return(value)

}