predict.multinomClust <- function(object,number.clusters=NULL,new_Y_data=NULL,new_X_data=NULL,...) {
  # number of clusters
  if (is.null(number.clusters)) number.clusters <- object$number.clusters

  # make likelihood
  logL <- function(mydata) {
    prob <- apply(mydata[,object$category.first[1]:object$category.last[length(object$category.last)]],2,mean)
    Nparm <- sum(prob[object$category.first[1]:object$category.last[length(object$category.last)-1]]>0)
    res <- sum(prob[prob>0]*log(prob[prob>0]))*nrow(mydata)
    return(list(logL=res,Nparm=Nparm))
  }

  # initialize cluster estimates
  nlogL       <- rep(0,number.clusters)
  Nparm       <- rep(0,number.clusters)
  for (i in 1:number.clusters) {
    tmp <- logL(object$test.data[object$clusters.sequence[[number.clusters]][[i]],])
    nlogL[i] <- -tmp$logL
    Nparm[i] <- tmp$Nparm
  }

  # prediction within dataset
  if (is.null(new_Y_data)) {
    # initialize result matrix
    res <- matrix(0,nrow(object$Y_data),2+number.clusters)
    colnames(res) <- c("id","from.cluster",paste("prob",1:number.clusters,sep="."))

    # insert probabilities
    ii <- 1
    for (i in 1:number.clusters) for (j in object$clusters.sequence[[number.clusters]][i][[1]]) {
      res[ii,1] <- j
      res[ii,2] <- i

      # make log(likelihood ratio)
      for (k in 1:number.clusters) {
        kk <- object$clusters.sequence[[number.clusters]][k][[1]]
        if (k==i) kk <- setdiff(kk,j)
        # multinomial part
        if (k==i) res[ii,2+k] <- nlogL[k] + logL(object$test.data[kk,])$logL
        if (k!=i) res[ii,2+k] <- -logL(object$test.data[c(j,kk),])$logL - nlogL[k]
        # logistic regression part
        Y <- c(0,rep(1,length(kk)))
        tmp <- fitted.values(glm(reformulate(strsplit(deparse(object$logit.model[[3]]), " \\+ ")[[1]],response=Y),
                                 binomial,object$X_data[c(j,kk),]))
        tmp[1] <- 1-tmp[1]
        res[ii,2+k] <- res[ii,2+k]+sum(log(tmp))
      }

      # convert into probabilities
      res[ii,-(1:2)] <- res[ii,-(1:2)] - min(res[ii,-(1:2)])
      res[ii,-(1:2)] <- exp(-res[ii,-(1:2)])
      res[ii,-(1:2)] <- res[ii,-(1:2)]/sum(res[ii,-(1:2)])

      # step forward
      ii <- ii+1
    }
  }

  # prediction on new data
  if (!is.null(new_Y_data)) {
    # validate data
    if (ncol(new_Y_data)!=ncol(object$Y_data)) stop("new_Y_data must have the same number of variables as Y_data")
    if (any(colnames(new_Y_data)!=colnames(object$Y_data))) stop("new_Y_data must have the same variables as Y_data")
    for (i in 1:ncol(new_Y_data)) {
      if (any(!is.factor(new_Y_data[,i]))) stop("variables in new_Y_data must be factors")
      if (any(length(levels(new_Y_data[,i]))!=length(levels(object$Y_data[,i])))) stop("new_Y_data variables must have the same number of levels as Y_data")
      if (any(levels(new_Y_data[,i])!=levels(object$Y_data[,i]))) stop("new_Y_data variables must have the same levels as Y_data")
    }
    if (is.null(new_X_data)) new_X_data <- data.frame(intercept=rep(1,nrow(new_Y_data)))
    if (ncol(new_X_data)!=ncol(object$X_data)) stop("new_X_data must have the same number of variables as X_data")
    if (any(colnames(new_X_data)!=colnames(object$X_data))) stop("new_Y_data must have the same variables as Y_data")

    # make new.test.data frame
    new.test.data <- matrix(0,nrow(new_Y_data),object$category.last[length(object$category.last)])
    # encoding of non-missing events
    for (ii in 1:ncol(new_Y_data)) {
      tmp <- levels(new_Y_data[,ii])
      for (jj in 1:length(tmp)) new.test.data[new_Y_data[,ii]==tmp[jj],object$category.first[ii]+jj-1] <- 1
    }
    # encoding of missing events
    for (ii in 1:ncol(new_Y_data)) {
      new.test.data[is.na(new_Y_data[,ii]),object$category.last[length(object$category.last)-1]+ii] <- 1
    }
    # as data frame
    new.test.data <- as.data.frame(new.test.data)
    for (ii in 1:ncol(new_Y_data)) {
      colnames(new.test.data)[object$category.first[ii]:object$category.last[ii]] <- levels(new_Y_data[,ii])
    }
    colnames(new.test.data)[object$category.first[length(object$category.first)]:object$category.last[length(object$category.last)]] <- object$category.names


    # initialize result matrix
    Nparm.pair <- rep(0,number.clusters)
    res <- matrix(0,nrow(new_Y_data),2+number.clusters)
    colnames(res) <- c("id","pvalue",paste("prob",1:number.clusters,sep="."))

    # insert probabilities
    for (ii in 1:nrow(new_Y_data)) {
      res[ii,1] <- ii

      # number of parameters for ii
      Nparm.ii <- logL(new.test.data[ii,])$Nparm

      # make log(likelihood ratio)
      for (k in 1:number.clusters) {
        kk <- object$clusters.sequence[[number.clusters]][k][[1]]
        # multinomial part
        tmp <- logL(rbind(new.test.data[ii,],object$test.data[kk,]))
        res[ii,2+k] <- -tmp$logL - nlogL[k]
        Nparm.pair[k] <- tmp$Nparm
        # logistic regression part
        Y <- c(0,rep(1,length(kk)))
        tmp <- glm(reformulate(strsplit(deparse(object$logit.model[[3]]), " \\+ ")[[1]],response=Y),
                   binomial,rbind(new_X_data[ii,],object$X_data[kk,]))
        Nparm.pair[k] <- Nparm.pair[k] + sum(names(coef(tmp))=="(Intercept)") - sum(!is.na(coef(tmp)))
        tmp <- fitted.values(tmp)
        tmp[1] <- 1-tmp[1]
        res[ii,2+k] <- res[ii,2+k]+sum(log(tmp))
      }

      # find p-value
      res[ii,2] <- max(1-pchisq(2*res[ii,-(1:2)],pmax(0.01,Nparm.ii+Nparm-Nparm.pair)))
      #res[ii,2] <-min(1,exp(-min(res[ii,-(1:2)])))

      # convert into probabilities
      res[ii,-(1:2)] <- res[ii,-(1:2)] - min(res[ii,-(1:2)])
      res[ii,-(1:2)] <- exp(-res[ii,-(1:2)])
      res[ii,-(1:2)] <- res[ii,-(1:2)]/sum(res[ii,-(1:2)])
    }
  }

  # return result
  return(res)
}
