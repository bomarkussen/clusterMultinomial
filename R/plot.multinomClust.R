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
    # Define list of colors via significance level
    col.list <- c("red","green")[1+(p.adjust(x$p.list,pval.correction)>significance.level)]
    col.list[x$LR.list<0] <- "blue"

    # Initialize likelihood positions
    x.nlogL  <- x$nlogL.final-x$nlogL.initial
    y.pos    <- order(unlist(x$clusters))
    cluster2 <- x$clusters.initial

    # Make window
    x.length <- max(x.nlogL)
    plot(c(-0.1*x.length,1.1*x.length),
         c(0,length(y.pos)),type="n",axes=FALSE,
         xlab="log(likelihood ratio)",ylab="observation index")
#         main="Joint likelihood: red=(*), green=NS, blue=increase")
    tmp <- axTicks(1)
    i <- which.min(abs(max(x$nlogL.final)-tmp))
    if (max(x$nlogL.final)<tmp[i]) {tmp <- tmp[-1]; i <- i-1} #else {tmp <- tmp[-length(tmp)]}
    axis(1,max(x$nlogL.final)+tmp-tmp[i],labels=paste(abs(tmp-tmp[i])))

    # Add names?
    if (is.logical(add.names)) {
      if (add.names) text(x.nlogL[unlist(x$clusters)],1:length(y.pos),paste(unlist(x$clusters)),pos=4)
    } else {
      text(x.nlogL[unlist(x$clusters)],1:length(y.pos),add.names[unlist(x$clusters)],pos=4)
    }

    # Add lines
    for (i in 1:length(x$LR.list)) {
      group1 <- cluster2[[x$i.list[i]]]
      group2 <- cluster2[[x$j.list[i]]]
      pos1   <- unique(y.pos[group1])
      pos2   <- unique(y.pos[group2])
      nlogL1 <- unique(x.nlogL[group1])
      nlogL2 <- unique(x.nlogL[group2])
      new.nlogL <- nlogL1-abs(x$LR.list[i])

      # Plot lines
      lines(c(nlogL1,new.nlogL),rep(pos1,2),col=col.list[i])
      lines(c(nlogL2,new.nlogL),rep(pos2,2),col=col.list[i])
      lines(rep(new.nlogL,2),c(pos1,pos2),col=col.list[i])

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

    # empty result
    m <- NULL
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
    theme_set(theme_bw())
    m <- ggplot(pval,aes(x=1/theoretical,y=1/observed,color=accepted)) + geom_point() +
      scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
      scale_color_manual(values=c("red","green","blue")) +
      geom_abline(intercept=0,slope=1,col="green") +
      xlab("1/uniform distribution") + ylab("1/observed p-values")
  }


  # histograms -------------------------
  if (is.element(3,what)) {
    # cluster, category, response, proportion
    mydf <- data.frame(cluster=NULL,category=NULL,response=NULL,prop=NULL)
    for (i in 1:number.clusters) {
      for (j in 1:ncol(x$Y_data)) {
        tmp <- x$test.data[x$clusters.sequence[[number.clusters]][[i]],x$category.first[j]:x$category.last[j],drop=FALSE]
        tmp <- apply(tmp,2,mean)
        mydf <- rbind(mydf,data.frame(cluster=i,category=x$category.names[j],response=names(tmp),prop=tmp))
      }
    }

    theme_set(theme_bw())

    mydf2 <- subset(mydf,(response!=0)&(is.element(category,variables)))
    mydf2$category <- droplevels(mydf2$category)
    mydf2$response <- factor(mydf2$response,levels=unique(mydf2$response))
    m <- ggplot(mydf2,aes(x=response,y=prop),ylim=c(0,1)) + scale_y_continuous(breaks=seq(0,1,1),labels=scales::percent) +
      facet_grid(cluster~category,scales="free_x",space="free_x") + geom_col() +
      ylab("Proportion within cluster") + xlab("") + theme(axis.text.x=element_text(angle=90,vjust=0.5))
  }

  # return NULL
  return(m)
}
