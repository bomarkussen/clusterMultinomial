plot.multinomClust <- function(x,what="dendrogram",add.names=TRUE,number.clusters=NULL,variables=NULL,pval.correction="holm",significance.level=0.05,...) {
  # what to plot?
  what <- charmatch(what,c("dendrogram","pvalues","histogram"),0)

  # number of clusters
  if (is.null(number.clusters)) number.clusters <- x$number.clusters

  # variables in histogram
  if (is.numeric(variables)) variables <- names(x$Y_data)[variables]
  if (is.null(variables)) variables <- names(x$Y_data)


  # dendrogram ----------------------------------
  if (is.element(1,what)) {
    # Define list of colors via significance level: red=1, green=2, blue=3
    col.list <- 1+(p.adjust(x$p.list,pval.correction)>significance.level)
    col.list[x$LR.list<0] <- 3

    # Initialize likelihood positions
    x.nlogL  <- x$nlogL.final-x$nlogL.initial
    y.pos    <- order(unlist(x$clusters))
    cluster2 <- x$clusters.initial

    # Make x scale
    tmp <- pretty(c(0,0.8*max(x.nlogL)))
    axis.x <- max(x$nlogL.final)-rev(tmp)
    axis.labels <- paste(rev(tmp))

    # Add names?
    mynames <- NULL
    if (is.logical(add.names)) {
      if (add.names) mynames <- data.frame(x=x.nlogL[unlist(x$clusters)],y=1:length(y.pos),label=paste(unlist(x$clusters)),hjust=0)
    } else {
      mynames <- data.frame(x=x.nlogL[unlist(x$clusters)],y=1:length(y.pos),label=add.names[unlist(x$clusters)],hjust=0)
    }

    # Make dendrogram
    mylines <- data.frame(x=NULL,xend=NULL,y=NULL,yend=NULL,col=NULL)
    for (i in 1:length(x$LR.list)) {
      group1 <- cluster2[[x$i.list[i]]]
      group2 <- cluster2[[x$j.list[i]]]
      pos1   <- unique(y.pos[group1])
      pos2   <- unique(y.pos[group2])
      nlogL1 <- unique(x.nlogL[group1])
      nlogL2 <- unique(x.nlogL[group2])
      new.nlogL <- nlogL1-abs(x$LR.list[i])

      # Add lines
      mylines <- rbind(mylines,
                       data.frame(x=c(nlogL1,nlogL2,new.nlogL),xend=rep(new.nlogL,3),y=c(pos1,pos2,pos1),yend=c(pos1,pos2,pos2),col=rep(col.list[i],3)))

      # Collapse clusters
      if (mean(x$ordering[cluster2[[x$i.list[i]]]]) > mean(x$ordering[cluster2[[x$j.list[i]]]])) {
        cluster2[[x$i.list[i]]] <- c(cluster2[[x$j.list[i]]],cluster2[[x$i.list[i]]])
      } else {
        cluster2[[x$i.list[i]]] <- c(cluster2[[x$i.list[i]]],cluster2[[x$j.list[i]]])
      }
      cluster2                       <- cluster2[-x$j.list[i]]
      x.nlogL[cluster2[[x$i.list[i]]]] <- new.nlogL
      y.pos[cluster2[[x$i.list[i]]]]   <- (pos1+pos2)/2
    }

    # return result
    m <- ggplot() +
      geom_segment(aes(x=x,y=y,xend=xend,yend=yend,col=factor(col)),mylines) +
      scale_color_manual(values=c("red","green","blue")) +
      scale_x_continuous(breaks=axis.x,labels=axis.labels) +
      scale_y_continuous(breaks=NULL) +
      xlab("log(likelihood ratio)") +
      ylab("observation index") +
      guides(col=FALSE) +
      theme_classic()
    if (!is.null(mynames)) m <- m + geom_text(aes(x=x,y=y,label=label,hjust=hjust),mynames)
  }


  # p-values --------------------------
  if (is.element(2,what)) {
    # Define list of colors via significance level
    col.list <- c("FALSE","TRUE")[1+(p.adjust(x$p.list,pval.correction)>significance.level)]
    col.list[x$LR.list<0] <- "TRUE (LR<0)"


    # compute p-values
    pval <- data.frame(theoretical=(length(x$p.list):1)/(1+length(x$p.list)),
                       observed=x$p.list,
                       accepted=col.list)

    # make plot of p-values
    #theme_set(theme_bw())
    m <- ggplot(pval,aes(x=1/theoretical,y=1/observed,color=accepted)) + geom_point() +
      scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
      scale_color_manual(values=c("red","green","blue")) +
      geom_abline(intercept=0,slope=1,col="green") +
      xlab("1/uniform distribution") + ylab("1/observed p-values") +
      theme_bw()
  }


  # histograms -------------------------
  if (is.element(3,what)) {
    # cluster, category, response, proportion
    mydf <- data.frame()
    for (i in 1:number.clusters) {
      for (j in 1:ncol(x$Y_data)) {
        tmp <- x$test.data[x$clusters.sequence[[number.clusters]][[i]],x$category.first[j]:x$category.last[j],drop=FALSE]
        tmp <- apply(tmp,2,mean)
        mydf <- rbind(mydf,data.frame(cluster=i,
                                      category=x$category.names[j],response=names(tmp),
                                      x=paste(x$category.names[j],names(tmp),sep=":"),prop=tmp))
      }
    }

    mydf2 <- subset(mydf,is.element(category,variables))
    mydf2$category <- droplevels(mydf2$category)
    mydf2$x <- factor(mydf2$x,levels=unique(mydf2$x))

    m <- ggplot(mydf2,aes(x=x,y=prop),ylim=c(0,1)) + scale_y_continuous(breaks=seq(0,1,1),labels=scales::percent) +
      facet_grid(cluster~category,scales="free_x",space="free_x") + geom_col() +
      ylab("Proportion within cluster") + xlab("") +
      scale_x_discrete(breaks=mydf2$x,labels=mydf2$response) +
      theme_bw() + theme(axis.text.x=element_text(angle=90,vjust=0.5))
  }

  # return ggplot-object
  return(m)
}
