cluster <- function(Y_data,X_data=NULL,logit.model=NULL,ordering=NULL,method="pvalue",verbose=FALSE) {
  # validate input ----------------------------

  if (is.null(X_data)) X_data <- as.data.frame(matrix(0,nrow(Y_data),0))
  if (is.null(logit.model)) logit.model <- .~1
  if (is.null(ordering)) ordering <- rep(1,nrow(Y_data))

  if (!is.data.frame(Y_data)) stop("Y_data argument (responses) must be a data frame")
  if (!is.data.frame(X_data)) stop("X_data argument (predictors) must be a data frame")
  if (nrow(Y_data)!=nrow(X_data)) stop("Y_data and X_data must have the same number of rows")
  if (class(logit.model)!="formula") stop("logit.model argument must be a formula")
  if (!is.numeric(ordering)) stop("order.model argument must be a formula")
  for (ii in 1:ncol(Y_data)) {if (!is.factor(Y_data[,ii])) stop("All columns in Y_data must be factors")}

  # extract category names and column range
  # also make columns for encoding of missingness
  category.names <- colnames(Y_data)
  category.last <- cumsum(c(apply(Y_data,2,function(x){length(unique(x[!is.na(x)]))}),
                            ncol(Y_data)))
  category.first <- c(1,1+category.last[-length(category.last)])
  names(category.first) <- names(category.last)

  # make test data frame
  test.data <- matrix(0,nrow(Y_data),category.last[length(category.last)])
  # encoding of non-missing events
  for (ii in 1:ncol(Y_data)) {
    tmp <- levels(Y_data[,ii])
    for (jj in 1:length(tmp)) test.data[Y_data[,ii]==tmp[jj],category.first[ii]+jj-1] <- 1
  }
  # encoding of missing events
  for (ii in 1:ncol(Y_data)) {
    test.data[is.na(Y_data[,ii]),category.last[length(category.last)-1]+ii] <- 1
  }
  # as data frame
  test.data <- as.data.frame(test.data)
  for (ii in 1:ncol(Y_data)) {
    colnames(test.data)[category.first[ii]:category.last[ii]] <- levels(Y_data[,ii])
  }
  colnames(test.data)[category.first[length(category.first)]:category.last[length(category.last)]] <- category.names


  # make likelihood
  logL <- function(mydata) {
    prob <- apply(mydata[,category.first[1]:category.last[length(category.last)]],2,mean)
    Nparm <- sum(prob[category.first[1]:category.last[length(category.last)-1]]>0)
    res <- sum(prob[prob>0]*log(prob[prob>0]))*nrow(mydata)
    return(list(logL=res,Nparm=Nparm))
  }


  # Initialize likelihoods and clusters ---------------------
  Ntemp       <- nrow(test.data)
  nlogL       <- rep(0,Ntemp)
  nlogL.pair  <- array(Inf,dim=c(Ntemp,Ntemp))
  Nparm       <- rep(0,Ntemp)
  Nparm.pair  <- array(0,dim=c(Ntemp,Ntemp))

  i.list        <- NULL
  j.list        <- NULL
  LR.list       <- NULL
  df.list       <- NULL
  p.list        <- NULL
  nlogL.list    <- NULL            # -logL(pair) +  logL(pair separation)
  nlogL.initial <- rep(0,Ntemp)    # Initial -logL
  nlogL.final   <- rep(0,Ntemp)    # Final   -logL

  clusters.initial <- as.list(1:Ntemp)  # Initial clusters
  clusters <-         clusters.initial  # Initialize clusters

  clusters.sequence          <- as.list(1:Ntemp) # Used to collect sequence of clusters: added December 21, 2012
  clusters.sequence[[Ntemp]] <- clusters


  # Initialize fits with one response ----------------------------
  if (verbose) cat("Separate fits of observations\n")
  for (i in 1:Ntemp) {
    tmp <- logL(test.data[i,])
    nlogL[i] <- -tmp$logL
    Nparm[i] <- tmp$Nparm
  }

  # Fit all selections of two responses --------------------------
  if (verbose) cat("Initialize pairwise fits observations\n")
  for (i in 1:(Ntemp-1)) {
    if (verbose) cat("fit i=",i,"out of",Ntemp,"\n")
    for (j in (i+1):Ntemp) {
      ii  <- c(clusters[[i]],clusters[[j]])
      tmp <- logL(test.data[ii,])
      nlogL.pair[i,j] <- -tmp$logL
      Nparm.pair[i,j] <- tmp$Nparm
      # Add Nparm and logL for separation of pair
      Y <- c(rep(0,length(clusters[[i]])),rep(1,length(clusters[[j]])))
      tmp <- glm(reformulate(strsplit(deparse(logit.model[[3]]), " \\+ ")[[1]],response=Y),binomial,X_data[ii,])
      Nparm.pair[i,j] <- Nparm.pair[i,j] + sum(names(coef(tmp))=="(Intercept)") - sum(!is.na(coef(tmp)))
      tmp <- fitted.values(tmp)
      tmp[1:length(i)] <- 1-tmp[1:length(i)]
      nlogL.pair[i,j]  <- nlogL.pair[i,j]+sum(log(tmp))
    }
  }

  # Agglomerate successively -------------------------------------
  while (Ntemp > 1) {
    # Select minimum and collapse nlogL array and clusters
    if (method=="pvalue") {
      ii <- which.min(pchisq(2*nlogL.pair-2*outer(nlogL,nlogL,'+'),pmax(0.01,outer(Nparm,Nparm,'+')-Nparm.pair),log.p=TRUE))
    } else {
      ii <- which.min(nlogL.pair-outer(nlogL,nlogL,'+'))
    }
    j.min <- ceiling(ii/Ntemp)
    i.min <- ii-j.min*Ntemp+Ntemp
    tmp   <- nlogL.pair[ii]-nlogL[i.min]-nlogL[j.min]
    tmp2  <- max(0,Nparm[i.min]+Nparm[j.min]-Nparm.pair[ii])
    if (i.min >= j.min) stop("Error i=",i.min," >= ",j.min,"=j")

    # update results
    i.list     <- c(i.list,i.min)
    j.list     <- c(j.list,j.min)
    LR.list    <- c(LR.list,tmp)
    df.list    <- c(df.list,tmp2)
    p.list     <- c(p.list,1-pchisq(2*tmp,tmp2))
    nlogL.list <- c(nlogL.list,nlogL.pair[ii])

    # Add negative -logL's
    nlogL.initial[c(clusters[[i.min]],clusters[[j.min]])] <-
      nlogL.initial[c(clusters[[i.min]],clusters[[j.min]])] + min(0,tmp)
    nlogL.final[c(clusters[[i.min]],clusters[[j.min]])] <-
      nlogL.final[c(clusters[[i.min]],clusters[[j.min]])] + max(0,tmp)

    if (verbose) {
      cat("clusters left=",Ntemp,": collapsing",i.min,"and",j.min,": logL=",nlogL.pair[ii],
          ", log(LR)=",tmp,", df=",tmp2,", p=",1-pchisq(2*tmp,tmp2),"\n")
    }

    # Collapse clusters, nlogL and nlogL.pairwise
    if (mean(ordering[clusters[[i.min]]]) > mean(ordering[clusters[[j.min]]])) {
      clusters[[i.min]] <- c(clusters[[j.min]],clusters[[i.min]])
    } else {
      clusters[[i.min]] <- c(clusters[[i.min]],clusters[[j.min]])
    }
    clusters          <- clusters[-j.min]
    nlogL             <- nlogL[-j.min]
    tmp               <- logL(test.data[clusters[[i.min]],])
    nlogL[i.min]      <- -tmp$logL
    nlogL.pair        <- nlogL.pair[-j.min,-j.min]
    Nparm             <- Nparm[-j.min]
    Nparm[i.min]      <- tmp$Nparm
    Nparm.pair        <- Nparm.pair[-j.min,-j.min]
    Ntemp             <- Ntemp-1

    # Update list of consequtive clusters, added December 21, 2012
    clusters.sequence[[Ntemp]] <- clusters


    # Compute new pairs
    for (j in (1:Ntemp)[-i.min]) {
      ii  <- c(clusters[[i.min]],clusters[[j]])
      tmp <- logL(test.data[ii,])
      nlogL.pair[min(i.min,j),max(i.min,j)] <- -tmp$logL
      Nparm.pair[min(i.min,j),max(i.min,j)] <- tmp$Nparm
      # Add Nparm and logL for separation of pair
      Y <- c(rep(0,length(clusters[[i.min]])),rep(1,length(clusters[[j]])))
      tmp <- glm(reformulate(strsplit(deparse(logit.model[[3]]), " \\+ ")[[1]],response=Y),binomial,X_data[ii,])
      Nparm.pair[min(i.min,j),max(i.min,j)] <- Nparm.pair[min(i.min,j),max(i.min,j)] + sum(names(coef(tmp))=="(Intercept)") - sum(!is.na(coef(tmp)))
      tmp <- fitted.values(tmp)
      tmp[1:length(clusters[[i.min]])] <- 1-tmp[1:length(clusters[[i.min]])]
      nlogL.pair[min(i.min,j),max(i.min,j)] <- nlogL.pair[min(i.min,j),max(i.min,j)]+sum(log(tmp))
    }

    # End loop
  }


  # Reorder cluster sequence
  for (i in 1:(length(i.list)+1)) {
    tmp <- clusters.sequence[[i]]
    ii  <- order(match(sapply(clusters.sequence[[i]],function(x){x[[1]]}),clusters[[1]]))
    for (j in 1:i) clusters.sequence[[i]][[j]] <- tmp[[ii[j]]]
  }

  # number of clusters
  tmp <- p.adjust(rev(p.list),"holm") < 0.05
  if (!tmp[1])            number.clusters <- 1
  if (all(tmp))           number.clusters <- length(tmp)
  if (tmp[1] & any(!tmp)) number.clusters <- which.min(tmp)

  # return result -------------------
  res <- list(Y_data=Y_data,X_data=X_data,logit.model=logit.model,ordering=ordering,test.data=test.data,
              category.first=category.first,category.last=category.last,category.names=category.names,
              nlogL.initial=nlogL.initial,nlogL.final=nlogL.final,
              clusters.initial=clusters.initial,clusters=clusters,clusters.sequence=clusters.sequence,
              LR.list=LR.list,df.list=df.list,p.list=p.list,i.list=i.list,j.list=j.list,
              number.clusters=number.clusters)
  class(res) <- "multinomClust"
  return(res)
}
