
kohparse<-function(som,data,train.expr,data.missing=NA,data.threshold=c(-Inf,Inf),quiet=FALSE,n.cores=1,max.na.frac) { #{{{
  #Starting Prompt
  if (!quiet) { 
    timer<-proc.time()[3]
    cat("Running kohparse\n")
  }
  data.len<-nrow(data)
  if (data.len==0) { 
    stop("There is no data to parse! nrow(data)=0!") 
  }
  if (missing(train.expr)) { 
    if (!quiet) { 
      cat("Loading factor names from SOM\n")
    }
    train.expr<-colnames(som$codes[[1]])
  } else if (any(train.expr!=colnames(som$codes[[1]]))) { 
    cat(paste("Assuming provided train.expr entries map directly onto the training labels!\n",
            'i.e.\n',paste(collapse=' ',train.expr[1:min(length(train.expr),3)]),
            '\n         are the equivalent of \n',
            paste(collapse=' ',colnames(som$codes[[1]])[1:min(length(train.expr),3)]),
            '\n'))
    warning(paste("Assuming provided train.expr entries map directly onto the training labels!\n"))
    colnames(som$whiten.param)<-train.expr
  }
  #Scale the data for use in the SOM
  kohwhiten.output<-kohwhiten(data=data,train.expr=train.expr,whiten.param=som$whiten.param,
                              data.missing=data.missing,data.threshold=data.threshold)
  data.white<-kohwhiten.output$data.white
  whiten.param<-kohwhiten.output$whiten.param
  #Check for any data rows with fully NA data
  bad.data<-rowAlls(is.na(data.white))
  if (any(bad.data)) { 
    warning("There are fully NA rows in the whitened data!")
  }
  if (any(is.na(bad.data))) { 
    stop("there are NA's in the bad data check")
  }
  #check the max NA fraction parameter 
  if (!missing(max.na.frac)) { 
    if (length(som$maxNA.fraction)==0 & !missing(max.na.frac)) { 
      warning(paste("SOM has no maxNA.fraction parameter! Using requested frac:",max.na.frac))
      som$maxNA.fraction<-max.na.frac
    } else if (length(som$maxNA.fraction)==0 & missing(max.na.frac)) { 
      stop(paste("SOM has no maxNA.fraction parameter, and there isn't one specified! Please specify a value for max.na.frac!"))
    } else if (som$maxNA.fraction != max.na.frac) { 
      warning(paste("Overwriting SOM maxNA.fraction with requested frac:",som$maxNA.fraction,"->",max.na.frac))
      som$maxNA.fraction<-max.na.frac
    }
  }
  bad.data<-rowSums(ifelse(is.na(data.white),1,0))/ncol(data.white) > som$maxNA.fraction
  if (any(bad.data) & som$maxNA.fraction==1) { 
    warning("There are fully NA rows in the whitened data!")
  }
  if (!quiet) { 
    #close the progress bar and prompt
    cat("Running data prediction\n")
  }
  #Setup the number of processes
  if (n.cores<1) {
    num_splits<-try(as.numeric(system("grep -c ^processor /proc/cpuinfo",intern=T)))
    if (class(num_splits)!='try-error' || length(num_splits)==0 || !is.finite(num_splits)) { 
      cat("core number lookup failed. Running in serial!\n")
      warning("core number lookup failed. Running in serial!")
      num_splits<-1
    }
  } else { 
    num_splits<-n.cores
  }
  if (num_splits > 1) { 
    #Run the data prediction in parallel
    registerDoParallel(cores=num_splits)
    #Save the training classification
    if (length(som$training.classif)==0) { 
      som$training.classif<-som$unit.classif
      som$training.distances<-som$distances
    }
    #Remove data entry (can be large and is defunct after parse)
    if (length(som$data)!=0) { 
      for (i in 1:length(som$data)) { 
        som$data[[i]]<-som$data[[i]][1:2,]
      }
    }
    som$unit.classif<-rep(NA,nrow(data.white))
    som$distances<-rep(NA,nrow(data.white))
    good.unit.classif <-
     foreach(d=isplitRows(data.white[!bad.data,], chunks=num_splits),
     .combine=rbind, .packages=c("stats")) %dopar% {
        pred<-predict(som, newdata=d)
        return=cbind(unit.classif=pred$unit.classif,distances=pred$distances)
    }
    #Check that the parse worked 
    if (nrow(good.unit.classif)!=length(which(!bad.data))) {
      stop("Parse failed to predict entries for all good data")
    }
    som$unit.classif[which(!bad.data)]<-good.unit.classif[,1]
    som$distances[which(!bad.data)]<-good.unit.classif[,2]
  } else { 
    #Run the data prediction in serial 
    pred.som<-predict(som,newdata=data.white[!bad.data,])
    if (is.null(som$training.classif)) { 
      som$training.classif<-som$unit.classif
      som$training.distances<-som$distances
    }
    som$distances<-som$unit.classif<-rep(NA,nrow(data.white))
    som$unit.classif[which(!bad.data)]<-pred.som$unit.classif
    som$distances[which(!bad.data)]<-pred.som$distances
  } 
  #Check that the parse worked 
  if (length(som$unit.classif)!=nrow(data.white)) {
    stop("Parse failed to predict entries for all data")
  }
  #Save the data 
  som$data[[1]]<-data.white
  #Add a vector for training rows all good data
  som$good.phot<-rowAlls(!is.na(data.white))
  if (!is.null(som$n.cluster.bins)) { 
    cat("Generating kohgroups for previously defined number of clusters\n")
    som<-generate.kohgroups(som,n.cluster.bins=som$n.cluster.bins,n.cores=n.cores,quiet=quiet)
  }
  if (!quiet) { 
    #Prompt
    cat("Ending\n")
  }
  if (!quiet) { 
    #Closing prompt
    cat(paste0("Done in ",as.time(proc.time()[3]-timer),'\n'))
  }
  #Return the SOM
  return=som
}#}}}

generate.kohgroups<-function(som,n.cluster.bins=Inf,n.cores=1,new.data,subset,quiet=FALSE,
                             method='complete',hclust_in,...) { #{{{
  # The function generates SOM groupings based on the unit classifications
  # returned from the SOM training. If new.data is specified, this data is 
  # used to generate new unit classifications (additional options to kohparse 
  # can be supplied via "...") and generate groupings. 

  #Check if we're analysing a different dataset
  if (!missing(new.data)) { 
    som<-kohparse(som=som,data=new.data,quiet=quiet,n.cores=n.cores,...)
  }

  #Do we want a subset of the data? 
  if (missing(subset)) { 
    subset<-which(is.finite(som$unit.classif))
    somcells<-somclust<-som$unit.classif[subset]
  } else { 
    #If needed, convert subset from logical
    if (class(subset)=="logical") {
      subset<-which(subset)
    }
    #Is the subset the NULL set?
    if (length(subset)!=0) { 
      somclust<-somcells<-som$unit.classif[subset]
    } else { 
      warning("Subset provided to generate_kohgroups is the null set!")
      somclust<-NULL
    }
  }
  if (length(somclust)!=0) { 
    som.dim<-c(som$grid$xdim,som$grid$ydim)
    #Convert the IDs from data-to-SOMcell into data-to-cluster {{{
    if (n.cluster.bins<prod(som.dim)) { 
      if (!quiet) { 
        cat(paste0("Constructing ",n.cluster.bins," clusters\n"))
      }
      #There are fewer cluster bins than SOM cells
      if (missing(hclust_in) & length(som$hclust)==0) { 
        if (!quiet) { 
          cat("Generating hierarchical clustering\n")
        }
        #The cell clustering has not been done yet
        som$hclust<-hclust(dist(x=som$codes[[1]]),method=method)
      } else { 
        if (!quiet) { 
          if (!missing(hclust_in)) { 
            cat(paste0("Using input hierarchical clustering\n"))
            som$hclust<-hclust_in
          } else { 
            cat(paste0("Using existing hierarchical clustering in SOM structure\n"))
          }
        }
      }
      #Cut the hclust dendrogram at the desired number of groups
      som.hc = cutree(tree=som$hclust, k=n.cluster.bins)
      #Do we have data to group?
      if (length(somclust)!=0) {
        #Group the data
        somclust<-unsplit(som.hc,factor(somcells,levels=seq(prod(som.dim))),drop=FALSE)
      }
    } else { 
      #All cells are a cluster
      n.cluster.bins<-prod(som.dim)
      som.hc<-seq(n.cluster.bins)
    }
    #}}}
  } else { 
    #There is nothing to do
    som.hc<-NULL
  }
  #If using a subset, reconstruct the full unit.classif
  somclust.full<-rep(NA,length(som$unit.classif))
  somclust.full[subset]<-somclust
  somclust<-somclust.full
  #Update the SOM structure
  som$clust.classif<-somclust
  som$n.cluster.bins<-n.cluster.bins
  som$cell.clust<-som.hc
  #Return
  return=som
}#}}}

generate.kohgroup.property<-function(som,data,expression,expr.label=NULL,n.cores=1,n.cluster.bins,subset,quiet=FALSE,
                                     returnMatrix=FALSE,...) { #{{{
  
  if (missing(n.cluster.bins)) { 
    #If missing, read the n.cluster.bins
    n.cluster.bins<-som$n.cluster.bins
  }
  if (is.null(n.cluster.bins)) { 
    n.cluster.bins<-prod(som$grid$xdim,som$grid$ydim)
  }
  #Check that the group number is finite
  if (!is.finite(n.cluster.bins)) { 
    stop("The SOM has non-finite n.cluster.bins?! (If you wanted max bins, set to NULL)")
  }
  #convert the expression(s) to a single command
  if (length(expression)>1) { 
    expression<-paste0('c(',paste(expression,collapse=','),')')
  }
  #Setup the vector of groups 
  factors<-seq(n.cluster.bins)
  #Check that the SOM has the correct classifications
  if (nrow(data)!=length(som$unit.classif)) { 
    if (!quiet) { 
      cat("Data length is not equal to SOM unit classifications. Regenerating groups!\n") 
    }
    #Rerun the parse and grouping
    som<-generate.kohgroups(som=som,new.data=data,n.cluster.bins=n.cluster.bins,subset=subset,
                            quiet=quiet,n.cores=n.cores,...)
  } else if (nrow(data)!=length(som$clust.classif)) { 
    if (length(som$clust.classif)==0) { 
      if (!quiet) { 
        cat("There are no SOM cluster classifications. Generating groups for the first time!\n") 
      }
    } else { 
      if (!quiet) { 
        cat("Data length is not equal to SOM cluster classifications. Regenerating groups!\n") 
      }
    }
    #Rerun the grouping
    som<-generate.kohgroups(som=som,n.cluster.bins=n.cluster.bins,subset=subset,
                            quiet=quiet,n.cores=n.cores,...)
  } else if (n.cluster.bins!=som$n.cluster.bins) { 
    if (!quiet) { 
      cat("Requested n.cluster.bins is different to SOM n.cluster.bins. Regenerating groups!\n") 
    }
    #Rerun the grouping
    som<-generate.kohgroups(som=som,n.cluster.bins=n.cluster.bins,subset=subset,
                            quiet=quiet,n.cores=n.cores,...)
  }
  #Prepare the SOM groupings
  som.group<-som$clust.classif
  if (!missing(subset)) { 
    if (is.logical(subset)) { 
      som.group[!subset]<-NA
    } else { 
      som.group[!(1:length(som.group))%in%subset]<-NA
    }
  }
  
  #Prepare the expression per-group 
  expression<-gsub("data","data.tmp",expression)
  expression<-gsub("full.data.tmp","data",expression)
  ##Prepare the parallelisation
  #registerDoParallel(cores=n.cores)
  ##Run the expression per group
  #property<-foreach(i=factors,.combine=rbind,
  #               .export=c('som.group','expression','data'),
  #               .inorder=TRUE)%dopar%{
  #  group.index<-which(som.group==i)
  #  data.tmp<-data[group.index,]
  #  #evaluate the expression
  #  val<-eval(parse(text=expression)) 
  #  #Return the expression result(s)
  #  frame.tmp<-data.frame(group.id=i,rbind(val))
  #  #Prepare the expression label(s)
  #  if (length(expr.label)==0) { 
  #    colnames(frame.tmp)<-c("group.id",paste0("value.",seq(length(val))))
  #  } else { 
  #    if (length(val)!=length(expr.label)) { 
  #      stop(paste0("The expression returned more/less values than there are expression labels!\n",
  #                  paste(expr.label,collapse=' : '),"\n",length(val)))
  #    }
  #    colnames(frame.tmp)<-c("group.id",expr.label)
  #  }
  #  return=frame.tmp
  #}
  evalprop<-function(data.tmp) eval(parse(text=expression))
  som.group.fact<-factor(som.group,levels=seq(n.cluster.bins))
  property<-t(sapply(split(data,som.group.fact),evalprop))
  if (is.null(dim(property))) { 
    property<-cbind(property)
  }
  if (length(expr.label)==0) { 
    if (nrow(property)==n.cluster.bins) { 
      colnames(property)<-paste0('value.',1:ncol(property))
    } else if (ncol(property)==n.cluster.bins) { 
      property<-t(property)
      colnames(property)<-paste0('value.',1:ncol(property))
    } else { 
      stop(paste0("The output property(s) have no dimension (",dim(property)[1],
                  ",",dim(property)[2],") that matches the input number of clusters (",
                  n.cluster.bins,")?!"))
    }
  } else { 
    if (ncol(property)==length(expr.label)) { 
      colnames(property)<-expr.label
    } else if (nrow(property)==length(expr.label)) { 
      property<-t(property)
      colnames(property)<-expr.label
    } else { 
      warning("The expression labels are not the same length as the expression outputs!!")
    }
  }
  #if (bycell & (nrow(som$grid$pts)!=n.cluster.bins)) { 
  #  cellproperty<-unsplit(property,factor(som$cell.clust,levels=seq(nrow(som$grid$pts))))
  #}
  #Check for parallel errors {{{ 
  if (nrow(property)!=n.cluster.bins) { 
    stop("Error in parallelisation: Try Rerunning in serial!")
  } 
  if (returnMatrix) { 
    return=list(property=property,som=som)
  } else { 
    property<-as.data.frame(property)
    property$group.id<-1:n.cluster.bins
    return=list(property=as.data.frame(property),som=som)
  } 
  #}}}
}#}}}

kohgroup.loop<-function(som,data,expression,expr.label=NULL,n.cores=1,n.cluster.bins,quiet=FALSE,...) { #{{{
  
  if (missing(n.cluster.bins)) { 
    #If missing, read the n.cluster.bins
    n.cluster.bins<-som$n.cluster.bins
    if (is.null(n.cluster.bins)) { 
      n.cluster.bins<-prod(som$grid$xdim,som$grid$ydim)
    }
  }
  #Check that the group number is finite
  if (!is.finite(n.cluster.bins)) { 
    stop("The SOM has non-finite n.cluster.bins?!")
  }
  #convert the expression(s) to a single command
  if (length(expression)>1) { 
    expression<-paste0('cbind(',paste(expression,collapse=','),')')
  }
  #Setup the vector of groups 
  factors<-seq(n.cluster.bins)
  #Check that the SOM has the correct classifications
  if (nrow(data)!=length(som$unit.classif)) { 
    if (!quiet) { 
      cat("Data length is not equal to SOM unit classifications. Regenerating groups!\n") 
    }
    #Rerun the parse and grouping
    som<-generate.kohgroups(som=som,new.data=data,n.cluster.bins=n.cluster.bins,quiet=quiet,n.cores=n.cores,...)
  } else if (nrow(data)!=length(som$clust.classif)) { 
    if (!quiet) { 
      cat("Data length is not equal to SOM cluster classifications. Regenerating groups!\n") 
    }
    #Rerun the grouping
    som<-generate.kohgroups(som=som,n.cluster.bins=n.cluster.bins,quiet=quiet,n.cores=n.cores,...)
  } else if (n.cluster.bins!=som$n.cluster.bins) { 
    if (!quiet) { 
      cat("Requested n.cluster.bins is different to SOM n.cluster.bins. Regenerating groups!\n") 
    }
    #Rerun the grouping
    som<-generate.kohgroups(som=som,n.cluster.bins=n.cluster.bins,quiet=quiet,n.cores=n.cores,...)
  }
  #Prepare the SOM groupings
  som.group<-som$clust.classif
  
  #Prepare the expression per-group 
  expression<-gsub("data","data.tmp",expression)
  expression<-gsub("full.data.tmp","data",expression)
  #Prepare the parallelisation
  registerDoParallel(cores=n.cores)
  ##Run the expression per group
  #values<-foreach(i=factors,.combine=rbind,
  #               .export=c('som.group','expression','data'),
  #               .inorder=TRUE)%dopar%{
  #  group.index<-which(som.group==i)
  #  data.tmp<-data[group.index,]
  #  #evaluate the expression
  #  vals<-eval(parse(text=expression)) 
  #  #Return the expression result(s)
  #  frame.tmp<-data.frame(index=group.index,vals)
  #  #Prepare the expression label(s)
  #  if (length(expr.label)==0) { 
  #    colnames(frame.tmp)<-c("source.id",paste0("value.",ncol(frame.tmp)-1))
  #  } else { 
  #    if (is.null(dim(vals)) && length(vals)!=length(group.index)) { 
  #      stop(paste0("The expression returned more/less values than there are expression labels!\n",
  #                  paste(expr.label,collapse=' : '),"\n",length(vals)))
  #    } else if (nrow(vals)!=length(group.index)) { 
  #      stop(paste0("The expression returned more/less values than there are expression labels!\n",
  #                  paste(expr.label,collapse=' : '),"\n",length(vals)))
  #    }
  #    colnames(frame.tmp)<-c("group.id",expr.label)
  #  }
  #  return=frame.tmp
  #}
  evalprop<-function(data.tmp) eval(parse(text=expression))
  som.group.fact<-factor(som.group,levels=seq(n.cluster.bins))
  property<-rep(NA,length(som.group))
  if (length(evalprop(data[1:10,]))!=10) { stop("kohgroup.loop is only designed to return 1 expression at a time") }
  split(property,som.group.fact)<-sapply(split(data,som.group.fact),evalprop)
  #Check for parallel errors {{{ 
  #if (nrow(values)!=length(which(is.finite(som.group)))) { 
  #  stop("Error in parallelisation: Try Rerunning in serial!")
  #} 
  return=list(values=property,som=som)
  #}}}
}#}}}


