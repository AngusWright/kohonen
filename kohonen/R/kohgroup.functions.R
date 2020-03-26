
kohparse<-function(som,data,train.var.expr,data.missing=NA,data.threshold=c(-Inf,Inf),quiet=F,n.cores=1,max.na.frac=1) { #{{{
  #Starting Prompt
  if (!quiet) { 
    timer<-proc.time()[3]
    cat("Running kohparse\n")
  }
  data.len<-nrow(data)
  if (missing(train.var.expr)) { 
    if (!quiet) { 
      cat("Loading factor names from SOM\n")
    }
    train.var.expr<-colnames(som$codes[[1]])
  } else if (any(train.var.expr!=colnames(som$codes[[1]]))) { 
    cat(paste("Assuming provided train.var.expr entries map directly onto the training labels!\n",
            'i.e.\n',paste(collapse=' ',train.var.expr[1:min(length(train.var.expr),3)]),
            '\n         are the equivalent of \n',
            paste(collapse=' ',colnames(som$codes[[1]])[1:min(length(train.var.expr),3)]),
            '\n'))
    warning(paste("Assuming provided train.var.expr entries map directly onto the training labels!\n"))
  }
  #Scale the data for use in the SOM
  kohwhiten.output<-kohwhiten(data,train.expr,data.missing=data.missing,data.threshold=data.threshold)
  data.white<-kohwhiten.output$data.white
  whiten.param<-kohwhiten.output$whiten.param
  #check the max NA fraction parameter 
  if (length(som$maxNA.fraction)==0) { 
    warning(paste("SOM has no maxNA.fraction parameter! Using requested frac:",max.na.frac))
  } else if (som$maxNAfrac != max.na.frac) { 
    warning(paste("Overwriting SOM maxNA.fraction with requested frac:",som$max.na.frac,"->",max.na.frac))
  }
  if (!quiet) { 
    #close the progress bar and prompt
    cat("Running data prediction\n")
  }
  #Setup the number of processes
  if (n.cores<1) {
    num_splits<-try(as.numeric(system("grep -c ^processor /proc/cpuinfo",intern=T)))
    if (length(num_splits)==0||!is.finite(num_splits)) { 
      cat("core number lookup failed. Running in serial!\n")
      warning("core number lookup failed. Running in serial!")
      num_splits<-1
    }
  } else { 
    num_splits<-n.cores
  }
  if (class(num_splits)!='try-error') { 
    #Run the data prediction in parallel
    registerDoParallel(cores=num_splits)
    #Save the training classification
    if (length(som$training.classif)==0) { 
      som$training.classif<-som$unit.classif
    }
    som$unit.classif <-
     foreach(d=isplitRows(data.white, chunks=num_splits),
     .combine=c, .packages=c("stats")) %dopar% {
        return=predict(som, newdata=d)$unit.classif
    }
  } else { 
    #Run the data prediction in serial 
    pred.som<-predict(som,newdata=data.white)
    som$training.classif<-som$unit.classif
    som$unit.classif<-pred.som$unit.classif
  } 
  #Check that the parse worked 
  if (length(som$unit.classif)!=nrow(data.white)) {
    stop("Parse failed to predict entries for all data")
  }
  if (!quiet) { 
    #Prompt
    cat("Ending\n")
  }
  #Add a vector for training rows all good data
  som$good.phot<-rowAlls(!is.na(data.white))
  if (!quiet) { 
    #Closing prompt
    cat(paste0("Done in ",as.time(proc.time()[3]-timer),'\n'))
  }
  #Return the SOM
  return=som
}#}}}

generate.kohgroups<-function(som,n.cluster.bins=Inf,n.cores=1,new.data,subset,...) { #{{{
  # The function generates SOM groupings based on the unit classifications
  # returned from the SOM training. If new.data is specified, this data is 
  # used to generate new unit classifications (additional options to kohparse 
  # can be supplied via "...") and generate groupings. 

  #Check if we're analysing a different dataset
  if (!missing(new.data)) { 
    som<-kohparse(som=som,data=new.data,...)
  }

  #Do we want a subset of the data? 
  if (missing(subset)) { 
    somcells<-somclust<-som$unit.classif
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
      #There are fewer cluster bins than SOM cells
      if (length(som$hclust)==0) { 
        #The cell clustering has not been done yet
        som$hclust<-hclust(dist(x=som$codes[[1]]))
      }
      #Cut the hclust dendrogram at the desired number of groups
      som.hc = cutree(tree=som$hclust, k=n.cluster.bins)
      #Do we have data to group?
      if (length(somclust)!=0) {
        #Group the data in parallel
        registerDoParallel(cores=n.cores)
        somind<-foreach(i=seq(n.cluster.bins),.combine=rbind,
                        .export=c("somcells","som.hc"),.inorder=FALSE)%dopar%{ 
          #Assign the data to the cluster if it's SOMcell belongs to the cluster
          ind<-which(somcells%in%which(som.hc==i))
          return=cbind(ind,rep(i,length(ind)))
        } 
        if (length(somind)==0 || nrow(somind)!=length(somcells)) { 
          cat("Error in parallelisation! Must rerun in serial.")
          somind<-matrix(NA,ncol=2,nrow=length(somcells)) 
          for(i in seq(n.cluster.bins)){ 
            #Assign the data to the cluster if it's SOMcell belongs to the cluster
            ind<-which(somcells%in%which(som.hc==i))
            somind[ind,]<-cbind(ind,rep(i,length(ind)))
          } 
        }
        #Update the SOMclust vector
        somclust[somind[,1]]<-somind[,2]
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
  if (!missing(subset)) { 
    somclust.full<-rep(NA,length(som$unit.classif))
    somclust.full[subset]<-somclust
    somclust<-somclust.full
  } 
  #Update the SOM structure
  som$clust.classif<-somclust
  som$n.cluster.bins<-n.cluster.bins
  som$cell.clust<-som.hc
  #Return
  return=som
}#}}}

generate.kohgroup.property<-function(som,data,expression,expr.label=NULL,n.cores=1,n.cluster.bins,quiet=F,...) { #{{{
  
  if (missing(n.cluster.bins)) { 
    #If missing, read the n.cluster.bins
    group.max<-som$n.cluster.bins
    if (is.null(group.max)) { 
      group.max<-prod(som$grid$xdim,som$grid$ydim)
    }
  }
  #Check that the group number is finite
  if (!is.finite(group.max)) { 
    stop("The SOM has non-finite n.cluster.bins?!")
  }
  #convert the expression(s) to a single command
  if (length(expression)>1) { 
    expression<-paste0('c(',paste(expression,collapse=','),')')
  }
  #Setup the vector of groups 
  factors<-seq(group.max)
  #Check that the SOM has the correct classifications
  if (nrow(data)!=length(som$unit.classif)) { 
    if (!quiet) { 
      cat("Data length is not equal to SOM unit classifications. Regenerating groups!\n") 
    }
    #Rerun the parse and grouping
    som<-generate.kohgroups(som=som,new.data=data,n.cluster.bins=group.max,quiet=quiet,...)
  } else if (group.max!=som$n.cluster.bins) { 
    if (!quiet) { 
      cat("Requested n.cluster.bins is different to SOM n.cluster.bins. Regenerating groups!\n") 
    }
    #Rerun the grouping
    som<-generate.kohgroups(som=som,n.cluster.bins=group.max,quiet=quiet,...)
  }
  #Prepare the SOM groupings
  som.group<-som$unit.classif
  
  #Prepare the expression per-group 
  expression<-gsub("data","data.tmp",expression)
  expression<-gsub("full.data.tmp","data",expression)
  #Prepare the parallelisation
  registerDoParallel(cores=n.cores)
  #Run the expression per group
  property<-foreach(i=factors,.combine=rbind,
                 .export=c('som.group','expression','data'),
                 .inorder=TRUE)%dopar%{
    group.index<-which(som.group==i)
    data.tmp<-data[group.index,]
    #evaluate the expression
    val<-eval(parse(text=expression)) 
    #Return the expression result(s)
    frame.tmp<-data.frame(group.id=i,rbind(val))
    #Prepare the expression label(s)
    if (length(expr.label)==0) { 
      colnames(frame.tmp)<-c("group.id",paste0("value.",seq(length(val))))
    } else { 
      if (length(val)!=length(expr.label)) { 
        stop(paste0("The expression returned more/less values than there are expression labels!\n",
                    paste(expr.label,collapse=' : '),"\n",length(val)))
      }
      colnames(frame.tmp)<-c("group.id",expr.label)
    }
    return=frame.tmp
  }
  #Check for parallel errors {{{ 
  if (nrow(property)!=group.max) { 
    stop("Error in parallelisation: Try Rerunning in serial!")
  } 
  return=list(property=property,som=som)
  #}}}
}#}}}


