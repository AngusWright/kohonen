#
#
# Function for training a SOM given input data
#
#

kohtrain<-function(data,train.expr,
                   som.dim=c(10,10),som.topo='hexagonal',som.toroidal=TRUE,som.iter=100,
                   som.rate=c(0.05,0.01),n.cores=1,som.method='pbatch',max.na.frac=1,
                   train.sparse=FALSE,sparse.frac=0.1,sparse.min.density=3,sparse.var=NULL,
                   data.missing=NA,data.threshold=c(-Inf,Inf),
                   quiet=FALSE,seed,keep.data=TRUE,...) {
  #Function trains a SOM from input data. Performs a number of 
  #preparatory steps beforehand. 

  #Check input parameters {{{
  #Data {{{
  if (nrow(data)==0) { 
    stop("input data has no rows")
  }
  #}}}
  #Training expressions {{{
  if (length(train.expr)<=1) { 
    stop("there are 1 or fewer training expressions")
  }
  #}}}
  #SOM dim {{{
  if (length(som.dim)!=2){ 
    stop("som.dim must be of length 2")
  }
  if (any(!is.finite(som.dim)|som.dim<=0)){ 
    stop("som.dim values must be finite & > 0")
  }
  #}}}
  #Seed {{{ 
  if (!missing(seed)) { 
    if (is.finite(seed)) { set.seed(seed) }
  }
  #}}}
  #}}}
  #Notify /*fold*/ {{{
  if (!quiet) { 
    cat("    -> whitening the input data")
    short.timer<-proc.time()[3]
  }#/*fend*/}}}
  #Whiten the input data {{{
  kohwhiten.output<-kohwhiten(data=data,train.expr=train.expr,data.missing=data.missing,data.threshold=data.threshold)
  data.white<-kohwhiten.output$data.white
  whiten.param<-kohwhiten.output$whiten.param
  #}}}
  #Check for bad rows in the input data {{{
  row.warning<-''
  if (any(rowAlls(is.na(data.white)))) { 
    #There are rows which are all NA {{{
    if (!quiet) { 
      warning("There are fully NA rows in the input data after whitening! Removing them!")
      row.warning<-" [ Warning: There are fully NA rows in the training catalogue! Removing them! ]\n"
    }
    #Remove them {{{
    bad.index<-which(rowAlls(is.na(data.white)))
    data<-data[-bad.index,]
    data.white<-data.white[-bad.index,]
    #}}}
    #}}}
  } else { 
    bad.index<-NULL
  }
  #Check that we didn't lose all the data {{{
  if (nrow(data)==0) { 
    #The input data has zero rows?!
    stop("Training data is completely filled with NA rows after whitening!")
  }
  #}}}
  #}}} 
  #Notify /*fold*/ {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
    if (row.warning!='') {
      cat(paste0('\n',row.warning)) 
      cat(paste0("    -> constructing som grid (",som.dim[1],'x',som.dim[2],')'))
    } else { 
      cat(paste0("\n    -> constructing som grid (",som.dim[1],'x',som.dim[2],')'))
    }
    short.timer<-proc.time()[3]
  }#/*fend*/}}}
  #Create the SOM grid /*fold*/ {{{
  data.grid<-somgrid(xdim = som.dim[1], ydim=som.dim[2], topo=som.topo, toroidal=som.toroidal)
  #/*fend*/}}}
  #Generate the SOM /*fold*/ {{{
  if (train.sparse) { 
    #Generate a SOM from sparse data sampling /*fold*/ {{{
    #Notify /*fold*/ {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
    }#/*fend*/}}}
    #Calculate the SOM using a sparse sampling of the data /*fold*/ {{{
    #Determine the sparse sampling fraction {{{
    if (sparse.frac*nrow(data.white)<prod(som.dim,sparse.min.density)) { 
      if (!quiet) { 
        cat(paste0("- WARNING: Sparse sampling creates fewer than ",sparse.min.density," sources per SOM cell!\n"))
        cat(paste0("       |-> forcing at least ",sparse.min.density," sources per cell.\n"))
      }
      sparse.frac<-som.dim[1]*som.dim[2]*sparse.min.density/nrow(data.white)
      if (sparse.frac > 1) { 
        if (!quiet) { 
          cat("           |-> ERROR: the SOM grid is too fine to sparse sample. Using full data vector\n")
        }
        sparse.frac<-1
      } else { 
        if (!quiet) { 
          cat("    -> continuing constructing SOM from sparse sampling of data vector")
        }
      }
    }
    #}}}
    #Sparse Sample {{{
    if (!is.null(sparse.var)) {
      if (any(colnames(data)==sparse.var)) { 
        #Sparse sample using weights from [[sparse.var]] {{{
        if (!quiet) { 
          cat(paste("\n    -> constructing SOM from sampling of",sparse.var,"in data vector"))
        }
        #Sparse indicies are weighted random draws from rows of data {{{
        sparse.index<-sample(nrow(data.white),size=ceiling(sparse.frac*nrow(data.white)),prob=data[[sparse.var]])
        #}}}
        #}}}
      } else { 
        #Error: sparse.var is not found {{{
        stop(paste("Cannot sparse sample data vector using",sparse.var,"because it is not in the data vector!"))
        #}}}
      }
    } else {
      #Sparse sample randomly from data {{{
      if (!quiet) { 
        cat("\n    -> constructing SOM from sparse sampling of data vector")
      }
      sparse.index<-sample(nrow(data.white),size=ceiling(sparse.frac*nrow(data.white)))
      #}}}
    }
    #}}}
    #Train the SOM {{{
    train.som<-som(data.white[sparse.index,], grid=data.grid, rlen=som.iter, alpha=som.rate, cores=n.cores,
    mode=som.method,maxNA=max.na.frac,...)
    if (!keep.data) train.som$data<-NULL
    #}}}
    #}}}
    #/*fend*/}}}
  } else { 
    #Generate a SOM from the full data vector /*fold*/ {{{
    #Notify /*fold*/ {{{
    if (!quiet) { 
      cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0)," ")) 
      short.timer<-proc.time()[3]
      cat("\n    -> constructing SOM from full data vector")
    }#/*fend*/}}}
    #Calculate the SOM using the full data vector /*fold*/ {{{
    train.som<-try(som(data.white, grid=data.grid, rlen=som.iter, alpha=som.rate, cores=n.cores,
                mode=som.method,maxNA=max.na.frac,...))
    if (!keep.data) train.som$data<-NULL
    if (class(train.som)=='try-error') { 
      cat("Error in SOM training\n")
      cat("Input variables were:\n") 
      cat("data.white:\n") 
      print(str(data.white))
      cat("data.grid:\n") 
      print(str(data.grid))
      cat("som.iter:\n") 
      print(str(som.iter))
      cat("som.rate:\n") 
      print(str(som.rate))
      cat("n.cores:\n") 
      print(str(n.cores))
      cat("som.method:\n") 
      print(str(som.method))
      cat("max.na.frac\n") 
      print(str(max.na.frac))
      stop("Error in SOM training")
    }
    #}}}
    #}}}
  }
  #Initialise the n.cluster.bins variable (used elsewhere){{{
  train.som$n.cluster.bins<-prod(som.dim)
  #}}}
  #Save the bad indices {{{
  train.som$bad.index<-bad.index
  #}}}
  #Save the whitening parameters {{{
  train.som$whiten.param<-whiten.param
  #}}}
  #Save the training classifications {{{
  train.som$training.classif<-train.som$unit.classif
  #}}}
  #Save the max.na.fraction {{{ 
  train.som$maxNA.fraction.training<-train.som$maxNA.fraction
  #}}}
  #/*fend*/}}}
  #Notify /*fold*/ {{{
  if (!quiet) { 
    cat(paste(" - Done",as.time(proc.time()[3]-short.timer,digits=0),"\n")) 
  }#/*fend*/}}}
  #Return the trained som {{{
  return=train.som
  #}}}
}

kohwhiten<-function(data,train.expr,whiten.param,data.missing,data.threshold,factor.weight) {
  #Check for character columns /*fold*/ {{{
  seperated.labels<-unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",train.expr),' ')))
  seperated.labels<-seperated.labels[which(seperated.labels!="")]
  if (any(sapply(seperated.labels,function(C) class(try(silent=T,eval(parse(text=C)))))=='function')) { 
    func.ind<-which(sapply(seperated.labels,function(C) class(try(silent=T,eval(parse(text=C)))))=='function')
    seperated.labels<-seperated.labels[-func.ind]
    if (!quiet) { cat(" (removed functions from checks)") }
  }
  if (!any(class(data)=='data-table')) { 
    data<-as.data.table(data)
  }
  if (any(sapply(data[,seperated.labels,with=F],class)=='character')) { 
    if (!quiet) { cat(" (converting character cols to numeric!)") }
    for (i in sperated.labels[which(sapply(data[1,seperated.labels,with=F],class)=='character')]) { 
      data[,i,with=F]<-as.numeric(data[,i,with=F])
    }
  }
  #/*fend*/}}}
  #Prepare the whitened data matrix /*fold*/ {{{
  data.white<-matrix(NA,nrow=nrow(data),ncol=length(train.expr))
  colnames(data.white)<-train.expr
  #/*fend*/}}}
  #Check for provided whiten.param and factor.weight {{{
  if (!missing(whiten.param) & !missing(factor.weight)) { 
    warning("Both whiten.param and factor.weight are provided. For consistent usage, factor.weight will be ignored!")
  } 
  #}}}
  if (missing(whiten.param)) { 
    #Prepare the whitening parameters {{{
    whiten.param<-matrix(NA,nrow=2,ncol=length(train.expr))
    colnames(whiten.param)<-train.expr
    #Check for provided dimensional weights {{{
    if (!missing(factor.weight)) { 
      if (length(factor.weight)!=length(train.expr)) { 
        stop("provided dimensional weights are not the same length as the training expressions") 
      }
      if (any(!is.finite(factor.weight))) { 
        stop("Some/all provided dimensional weights are not finite?!") 
      } 
      if (any(factor.weight==0)) { 
        stop("Some/all provided dimensional weights are zero?!") 
      }
    } else { 
      #If none provided, set dimensional weights to 1 {{{
      factor.weight<-rep(1,length(train.expr)) 
      #}}}
    }
    #}}}
    #Set labels for dimensional weights 
    names(factor.weight)<-train.expr
    #}}}
  } else if (any(colnames(whiten.param)!=train.expr)) { 
    #Check that the whiten parameters match the training parameters {{{
    if (ncol(whiten.param)==length(train.expr)) { 
      warning(paste0("Whitening parameter names do not match training expressions!\n",
                     "They are the same length, so we are assuming that they match 1:1!\n"))
      colnames(whiten.param)<-train.expr
    } else { 
      stop(paste0("Whitening parameter names do not match training expressions!\n",
                  "They are not the same length, so we cannot continue!\n"))
    }
    #}}}
  }
  #Loop through the factor expressions (they could be columns or expressions!)/*fold*/ {{{
  for (factor.expr in train.expr) { 
    #Calculate the expression result
    white.value<-data[,eval(parse(text=factor.expr))]
    #seperate out the components
    seperated.labels<-unique((vecsplit(gsub('[-+*\\/\\)\\(]'," ",factor.expr),' ')))
    seperated.labels<-seperated.labels[which(seperated.labels!="")]
    if (any(sapply(seperated.labels,function(C) class(try(silent=T,eval(parse(text=C)))))=='function')) { 
      func.ind<-which(sapply(seperated.labels,function(C) class(try(silent=T,eval(parse(text=C)))))=='function')
      seperated.labels<-seperated.labels[-func.ind]
    }
    #Start with everything being OK, and then...
    data.beyond.thresh<-data.is.missing<-rep(FALSE,nrow(data))
    #...Loop through the components 
    for (label in seperated.labels) { 
      #Get the detected and non detected sources
      if (is.na(data.missing)) { 
        data.is.missing<-data.is.missing | is.na(data[[label]])
      } else if (is.infinite(data.missing)) { 
        data.is.missing<-data.is.missing | (is.infinite(data[[label]]))
      } else if (is.nan(data.missing)) { 
        data.is.missing<-data.is.missing | (is.nan(data[[label]]))
      } else { 
        data.is.missing<-data.is.missing | (data[[label]]==data.missing)
      } 
      data.beyond.thresh<-data.beyond.thresh | 
        ((data[[label]]>max(data.threshold) | data[[label]]<min(data.threshold)) & 
          !data.is.missing)
    }
    #Set the values with missing parts to dummy values 
    if (any(is.na(data.is.missing))) { 
      data.is.missing[which(is.na(data.is.missing))]<-TRUE
    }
    if (any(is.na(data.beyond.thresh))) { 
      data.beyond.thresh[which(is.na(data.beyond.thresh))]<-TRUE
    }
    white.value[data.is.missing]<-NA
    white.value[data.beyond.thresh]<-NA

    #Calculate the detected median and mad
    if (!is.na(whiten.param[1,factor.expr])) {
      med.tmp<-whiten.param[1,factor.expr]
      mad.tmp<-whiten.param[2,factor.expr]
    } else { 
      med.tmp<-median(white.value,na.rm=T)
      mad.tmp<-mad(white.value,na.rm=T)*factor.weight[factor.expr]
    }
    #Whiten the data
    white.value<-(white.value-med.tmp)/mad.tmp
    #Save the whitening parameters 
    whiten.param[,factor.expr]<-c(med.tmp,mad.tmp)
    #Save the whitened data 
    data.white[,factor.expr]<-white.value
  }
  #/*fend*/}}}
  #Output the whitened data and the whitening parameters {{{
  return=list(data.white=data.white,whiten.param=whiten.param) 
  #/*fend*/}}}
}
