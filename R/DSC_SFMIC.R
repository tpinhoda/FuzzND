#' DSC_SFMIC
#' @export
#' @import stream  methods foreach magrittr e1071 cluster stringr

## SFMiC - Management
DSC_SFMiC <- function(maxMiC = 30, m = 2, theta = 0.8,theta_adapt = 0.1 ,Theta = 0.8, thetaClass = 0.7,minNumSFMiC = 4, timeWindow = 10, P = 9,  ts =50, thresholdTemp = 100, k_temp=25, minWeight = 10,silhouetteThreshold = 0,fr_threshold = 0.6,description = NULL){



  if(!is.null(description)) desc <- description
  else desc <-"Fuzzy Micro-Cluster - Management"
  DSC_SFMiC <- fMicro_refClass$new(maxMiC = maxMiC, m = m, theta = theta, theta_adapt = theta_adapt, Theta = Theta, thetaClass=thetaClass, minNumSFMiC = minNumSFMiC, timeWindow = timeWindow , P = P, ts= ts, thresholdTemp = thresholdTemp, k_temp=k_temp , minWeight = minWeight, silhouetteThreshold = silhouetteThreshold, fr_threshold=fr_threshold )
  structure(list(description = desc, RObj = DSC_SFMiC), class = c("DSC_SFMiC","DSC_Micro","DSC_R","DSC"))
}

fMicro_refClass <-
  setRefClass("fuzzyMicroCluster",
              fields = list(
                ## args
                maxMiC = "numeric",
                minNumSFMiC = "numeric",
                m = "numeric",
                theta = "numeric", # add to existing FMiC
                theta_adapt = 'numeric',
                Theta = "numeric", # merge two FMiC
                thetaClass = "numeric",
                ## store microcluster
                allSFMiC = "list",
                microCenters = "matrix", # needed for Macro-Clustering
                microWeights = "numeric", # needed for Macro-Clustering
                temporaryMemory = "data.frame",
                numMiC = "numeric",
                prediction = "vector",
                currentTime = "numeric",
                sleepMemory = "list",
                timeWindow = "numeric",
                P = "numeric",
                ts = "numeric",
                thresholdTemp = "numeric",
                k_temp = "numeric",
                silhouetteThreshold = "numeric",
                fr_threshold = "numeric",
                PN = "numeric",
                minWeight = "numeric",
                evaluation_hist = "list"
              ),

              methods = list(
                initialize = function(
                  maxMiC = 10, m = 2, theta = 0.8, Theta = 0.8, theta_adapt = 1, thetaClass = 0.7,minNumSFMiC = 4, timeWindow = 10, P = 9, ts=50, thresholdTemp = 100, k_temp=25, minWeight = 10, silhouetteThreshold = 0, fr_threshold = 0.6
                ) {
                  maxMiC <<- maxMiC
                  minNumSFMiC <<- minNumSFMiC
                  m <<- m
                  theta <<- theta
                  Theta <<- Theta
                  theta_adapt <<- theta_adapt
                  thetaClass <<- thetaClass
                  timeWindow <<- timeWindow
                  P <<- P
                  ts <<- ts
                  minWeight <<- minWeight

                  silhouetteThreshold <<- silhouetteThreshold
                  allSFMiC <<- list()
                  microCenters <<- matrix()
                  microWeights <<- numeric()
                  temporaryMemory <<- data.frame()
                  thresholdTemp <<- thresholdTemp
                  fr_threshold <<- fr_threshold
                  k_temp <<- k_temp
                  numMiC <<- numeric()
                  prediction <<- vector()
                  currentTime <<- 0
                  sleepMemory <<- list()
                  PN <<- 0
                  evaluation_hist <<- list()
                  .self
                }
              )
  )

fMicro_refClass$methods(

  train = function(STREAM,training_number,k){
    training_points <- get_points(STREAM,n=training_number,class = TRUE)
    splitted_points <- split_class(training_points)                              #Separa os pontos por classe
    #Utilizar o cmeans em cada grupo de classe e retornar k microgrupos para cada classe
    splitted_microclusters <- lapply(splitted_points, function(class_set){
      if(nrow(class_set) > 0){
        list(cluster = cmeans ( class_set[,-ncol(class_set)], centers = k, iter.max=100, verbose=FALSE, dist="euclidean", method="cmeans", m=2, rate.par = NULL), class=as.character(class_set[1,ncol(class_set)]))
      }else NA
    })
    complete_cases <- !is.na(splitted_microclusters)
    splitted_microclusters <- splitted_microclusters[complete_cases]
    #Sumariza os grupos em micro-grupos
    if(length(splitted_microclusters) > 1){
      allSFMiC <<- unlist(lapply(splitted_microclusters, function(class_set_mc){
        apply(class_set_mc$cluster$centers,1,function(center){
          SFMiC(x = center, membership = 1,class_set_mc$class, time = 1)
        })
      }))
    }else{
      allSFMiC <<- apply(splitted_microclusters[[1]]$cluster$centers,1,function(center){
        SFMiC(x = center, membership = 1,splitted_microclusters[[1]]$class,time = 1)})

    }
    microCenters <<- t(sapply(allSFMiC,function(sfmic){sfmic$get_center()}))
    numMiC <<- nrow(microCenters)
  },


  predict = function(STREAM){
    evaluation <- 1
    eval_data <- c()
    memb_greater_theta <- c(theta)
    while(STREAM$state$counter <= nrow(STREAM$strm)){
      #Point arrival#################################################
      currentTime <<- currentTime +1
      data_point <- as.numeric(get_points(STREAM, n=1, class = TRUE))
      eval_data <- rbind(eval_data,data_point)
      class_teste <- data_point[length(data_point)]
      data_point <- data_point[-length(data_point)]
      #Calculate membership of point to all micro clusters###########
      all_memberships <- get_memberships(microCenters,data_point)
      #Find max  membership sum of microclusters from the same class
      maxMicroclusters <- find_maxClass(allSFMiC,all_memberships)
      if(maxMicroclusters["membership"] == 1){
        center <- colMeans(microCenters)
        mics_memberships <- get_memberships(microCenters,center)
        dispersion <- calculate_dpf(mics_memberships,center,microCenters)
        distance <- dist(rbind(center,data_point))
        class_mics <- c(1:nrow(microCenters))
        if(distance <= dispersion["dpf"]){
          prediction <<- c(prediction,maxMicroclusters["class"])
          class_memberships <- get_memberships(microCenters,data_point)
          add_SFMiC(data_point,class_memberships,class_mics)
        }else if(distance <= dispersion["dpf"]*2){
          prediction <<- c(prediction,maxMicroclusters["class"])
          #Concept drift - Create a new microcluster for the class
          create_SFMiC(data_point,maxMicroclusters["class"],currentTime)
          #Check maximun of class microclusters
          if(length(class_mics) == maxMiC){
            remove_oldest(class_mics)
          }
        }else{
          #Updating prediction Unknows data discovered
          prediction <<- c(prediction,"unknown")
          unknown_data <- c(data_point, timestamp = currentTime)
          ##Store in temporary memory outlier
          temporaryMemory <<- rbind(temporaryMemory,unknown_data)
          #Novelty Detection

          if(nrow(temporaryMemory) >= thresholdTemp){
            cluster <- cmeans (temporaryMemory[-ncol(temporaryMemory)], centers = k_temp, iter.max=100, verbose=FALSE, dist="euclidean", method="cmeans", m=2, rate.par = NULL)
            mean_sil <- calculate_sil(cluster,temporaryMemory[,-ncol(temporaryMemory)])
            histo <- table(factor(cluster$cluster, levels = c(1:k_temp) ))
            for(clus_index in c(1:k_temp)){
              clus_silhouette <- mean_sil[clus_index]
              if(clus_silhouette > silhouetteThreshold && histo[clus_index] >= minWeight){
                clus_center <- cluster$centers[clus_index,]
                clus_memberships <- cluster$membership[,clus_index]
                instances_cluster <- which(cluster$cluster == clus_index)
                #verify novelty pattern####################################################
                check_PN(clus_memberships,clus_center,instances_cluster)
              }
            }
          }
        }

        ## if maximum membership >= theta
      }else if(maxMicroclusters["membership"] >= theta){
        #Updating Prediction
        memb_greater_theta <- c(memb_greater_theta,maxMicroclusters["membership"])
        theta <<- mean(as.numeric(memb_greater_theta))-theta_adapt
        prediction <<- c(prediction,maxMicroclusters["class"])
        #Calculate membership to same class microclusters
        class_mics <- sameClass_sfmics(allSFMiC,maxMicroclusters["class"])
        class_memberships <- get_memberships(microCenters[class_mics,],data_point)
        if(max(class_memberships) >= thetaClass){
          #Add data_point by membership
          add_SFMiC(data_point,class_memberships,class_mics)
        } else {
          #Concept drift - Create a new microcluster for the class
          create_SFMiC(data_point,maxMicroclusters["class"],currentTime)
          #Check maximun of class microclusters
          if(length(class_mics) == maxMiC){
            remove_oldest(class_mics)
          }
        }
      }else {
        #Updating prediction Unknows data discovered
        prediction <<- c(prediction,"unknown")
        unknown_data <- c(data_point, timestamp = currentTime)
        ##Store in temporary memory outlier
        temporaryMemory <<- rbind(temporaryMemory,unknown_data)
        #Novelty Detection
        if(nrow(temporaryMemory) >= thresholdTemp){
          cluster <- cmeans (temporaryMemory[-ncol(temporaryMemory)], centers = k_temp, iter.max=100, verbose=FALSE, dist="euclidean", method="cmeans", m=2, rate.par = NULL)
          mean_sil <- calculate_sil(cluster,temporaryMemory[,-ncol(temporaryMemory)])
          histo <- table(factor(cluster$cluster, levels = c(1:k_temp) ))
          length(temporaryMemory)
          for(clus_index in c(1:k_temp)){
            clus_silhouette <- mean_sil[clus_index]
            if(clus_silhouette >= silhouetteThreshold && histo[clus_index] >= minWeight){
              clus_center <- cluster$centers[clus_index,]
              clus_memberships <- cluster$membership[,clus_index]
              instances_cluster <- which(cluster$cluster == clus_index)
              #verify novelty pattern####################################################
              check_PN(clus_memberships,clus_center,instances_cluster)
             # temporaryMemory <<- temporaryMemory[-instances_cluster,]
            }
          }
        }
      }


      #save in sleep memory sfmics that do not receive any example for a P period of time
      toSleep <- get_forgottenSFMiC(P)
      if(length(toSleep) > 0){ ##DEVO LIMITAR PARA O MINUCO DE MICRO GRUPOS????????????????????
        allSFMiC <<- allSFMiC[-toSleep]
        microCenters <<- microCenters[-toSleep,]
        numMiC <<- nrow(microCenters)
      }
      remove_fromTemporary()

      if(mod(currentTime,timeWindow)==0){
        #cat("\n Evaluation Moment: ",evaluation," size: ",length(allSFMiC))
        evaluation <- evaluation + 1
        evaluation_now <- list(centers = microCenters, predictions = prediction, data = eval_data, PN = PN)
        evaluation_hist <<- append(evaluation_hist, list(evaluation_now))
        eval_data <- c()
      }
    }

    ## Merging Step
    #if(mod(currentTime,100)==0)
    #  merge_SFMiC()
  },

#######################################################################################################################
####################################################Predict Utils######################################################
#######################################################################################################################

oldest_mic = function(indexes){

  allmean <- sapply(allSFMiC, function(sfmic){sfmic$meant})
  classmean <- allmean[indexes]
  minindex <- which.min(classmean)
  return(minindex[1])
},

split_class = function(dataset){
  #Cria uma lista de pontos separados pela classe
  class_list <-as.character(dataset[,ncol(dataset)])
  splitted <- split(dataset,class_list)
  return(splitted)
},
get_memberships = function(centers,p) {
  if(!is.null(nrow(centers))){
    distmatrix <- as.matrix(dist(rbind(centers, p))) ### distance matrix for p and all FMiC centers
    only_pdist <- distmatrix[nrow(centers)+1,-(nrow(centers)+1)]
  }else{
    distmatrix <- as.matrix(dist(rbind(centers, p))) ### distance matrix for p and all FMiC centers
    only_pdist <- distmatrix[2,-2]
  }
  return(calculate_memberships(only_pdist)) ## Calculating membership values for p in all availableFMiC
},
calculate_memberships = function(dxC){
  dxC <- dxC[which(dxC != 0)]
  membership <- sapply(dxC,function(dx){sum((dx/dxC)^(2/(m-1)))})
  return(1/membership)
},

add_SFMiC =  function(data_point,memberships,mics_class){
  ## Update all FMiC where membership is > 0.1
  ## tem que classificar essa bagaca!!!!!!!!!
  membindex <- 0
  sapply(mics_class,function(j){
    membindex <- membindex + 1
    if(memberships[membindex] > 0.1){
      microCenters[j,] <<- allSFMiC[[j]]$c
      allSFMiC[[j]]$addPoints(x = data_point, membership = memberships[membindex], time = currentTime, m)
    }
  })
},

merge_SFMiC = function(...){
  can_mergerSFMiCs <- rep(TRUE, times=numMiC)
  merged <- numeric()
  same_class <- lapply(1:numMiC,function(index){sameClass_sfmics(allSFMiC,allSFMiC[[index]]$class_id)})
  same_class <- lapply(1:numMiC, function(index){same_class[[index]][-which(same_class[[index]]==index)]})
  FRs <- lapply(1:numMiC, function(index){calculate_fuzzyRadio(index,same_class[[index]])})
  maxFRs <- sapply(FRs,function(fr){max(fr)})
  mergeable <- which(maxFRs > 1)
  #len_mergeable <- sapply(mergeable,function(index)length(same_class[[index]]))
  if(length(mergeable) > 0){
    sapply(mergeable,function(index){
      if(can_mergerSFMiCs[index] && numMiC > minNumSFMiC){
        numMiC <<- numMiC - 1
        merged <- c(merged,index)
        sapply(same_class[[index]],function(index_same){
          if(numMiC > minNumSFMiC){
            can_mergerSFMiCs[index_same] <- FALSE
            allSFMiC[[index_same]]$merge(allSFMiC[[index]])
          }
        })
      }

    })
    microCenters <<- t(sapply(allSFMiC,function(sfmic){sfmic$get_center()}))
  }

  # foreach(merge = can_mergerSFMiCs, j = 1:numMiC) %do%{
  #   if(merge){
  #     same_class <- sameClass_sfmics(allSFMiC,allSFMiC[[j]]$get_class())
  #     same_class <<- same_class[same_class!=j]
  #     if(length(same_class) > 4){
  #       FR <- calculate_fuzzyRadio(j,same_class)
  #       maxdpf <- max(FR)
  #       if(maxdpf >= 2 && numMiC > minNumSFMiC){
  #         ## merge os sfmics
  #         sapply(same_class, function(index){
  #           if(can_mergerSFMiCs[index]){
  #             allSFMiC[[index]]$merge(allSFMiC[[j]])
  #           }
  #         })
  #         microCenters <<- t(sapply(allSFMiC,function(sfmic){sfmic$get_center()}))
  #         can_mergerSFMiCs[same_class] <- FALSE
  #         can_mergerSFMiCs[j] <- FALSE
  #         merged <- c(merged,j)
  #         numMiC <<- numMiC - 1
  #       }
  #     }
  #   }
  # }
  #delete merged sfmics
  if(length(merged) > 0 ){
    #print(merged)
    allSFMiC <<- allSFMiC[-merged]
    microCenters <<- microCenters[-merged,]
    #numMiC <<- nrow(microCenters)
  }
},

calculate_fuzzyRadio = function(center_index,same_class){
  mic_dist <- sapply(same_class, function(index_class){ dist(rbind(microCenters[index_class],microCenters[center_index,]))})
  class_dpf <- sapply(same_class,function(index){allSFMiC[[index]]$get_dpf()})
  FR <- (allSFMiC[[center_index]]$get_dpf()+class_dpf)/mic_dist
  return(FR)
},

sameClass_sfmics = function(MICS,classe){
  #retorna a lista de ids de micro-grupos da mesma classe
  allClasses <- get_allClasses(MICS)
  ids_class <- which(allClasses == classe)
  return(ids_class)
},
get_forgottenSFMiC = function(P){
  alltimestamp <- sapply(allSFMiC, function(sfmic){sfmic$get_lastTimestamp()})
  alltimestamp <- currentTime - alltimestamp
  forgotten <- which(alltimestamp > P)
  return(forgotten)
},

get_microclusters = function(...) { microCenters },
get_microweights = function(...) { microWeights },

membership_byClass = function(MICS, memberships) {
  classes <- get_classes(MICS)
  allClasses <- get_allClasses(MICS)
  sameclass <- lapply(classes, function(classe){which(allClasses==classe)})
  class_memberships <- unlist(lapply(sameclass,function(row){sum(memberships[row])}))
  class_memberships <- cbind(membership = class_memberships, class = classes)
  return(class_memberships)
},

get_PNs = function(){
  allclasses <- get_allClasses(allSFMiC)
  pnbool <- str_detect(allclasses,"PN")
  PNs <- which(pnbool==TRUE)
  return(PNs)

},

get_classes = function(MICS){

  return(unique(sapply(MICS,function(sfmic){sfmic$get_class()})))
},


get_allClasses = function(MICS){

  return(sapply(MICS,function(sfmic){sfmic$get_class()}))
},

find_maxClass = function(MiCs, memberships){
  class_memberships <- membership_byClass(MiCs,memberships)
  maxInd <- which.max(class_memberships[,"membership"])
  maxMem <- class_memberships[maxInd,"membership"]
  maxClass <- class_memberships[maxInd,"class"]
  return(c(maxMem, maxClass))
},

create_SFMiC = function(x,class,time){
  newSFMiC <- SFMiC(x = x, membership = 1, class = class,time = time)
  #calculating initial SSD
  newdist <- as.matrix(dist(rbind(microCenters, x)))[nrow(microCenters)+1,-(nrow(microCenters)+1)] ### distance matrix for p and all FMiC centers
  newmemb <- calculate_memberships(newdist)
  newdpf <- min(newdist)*newmemb[which.min(newdist)]
  newSFMiC$dpf <- newdpf

  #Updating SFMiCs
  allSFMiC <<- c(allSFMiC, newSFMiC)
  microCenters <<- rbind(microCenters, x)
  numMiC <<- nrow(microCenters)

},

remove_oldest = function(indexes){
  #Find oldest
  remove_index <- oldest_mic(indexes)
  #Updating SFMiCs
  allSFMiC <<- allSFMiC[-indexes[remove_index]]
  microCenters <<- microCenters[-indexes[remove_index],]
  numMiC <<- nrow(microCenters)
},

get_sfmiccenters = function(MiCs){
  return(t(sapply(MiCs, function(sfmic){sfmic$get_center()})))
},

is.empty = function(list){
  if(length(list) > 0)
    return(FALSE)
  else
    return(TRUE)
},

calculate_dpf = function(clus_memberships, clus_center, centers){
  clus_memberships <- clus_memberships^m
  distance <- apply(centers,1, function(row){sum((row-clus_center)^2)})

    if(length(clus_memberships) != length(distance))
      stop()
    ssde <- sum(clus_memberships*distance)
    dpf <- sqrt(ssde/nrow(centers))
    return(c(dpf = dpf, ssde = ssde))

},

calculate_pnFR = function(clus_center, clus_dispersion, all_PNs){
  FR <- c()
  if(!is.empty(all_PNs)){
    all_PNsCenters <- microCenters[all_PNs,]
    all_ssde <- sapply(allSFMiC, function(sfmic){sfmic$get_dpf()})
    all_PNsSsde <- all_ssde[all_PNs]
    all_distance <- as.matrix(dist(rbind(all_PNsCenters, clus_center)))[nrow(all_PNsCenters)+1,-(nrow(all_PNsCenters)+1)]
    FR <- (clus_dispersion["dpf"]+all_PNsSsde)/all_distance
  }
  return(FR)
},

create_pnSFMiC = function(clus_center,maxClass,ssde){

  allSFMiC <<- c(allSFMiC, SFMiC(x = clus_center, membership = 1, class = maxClass,time = currentTime, dpf=ssde))
  microCenters <<- rbind(microCenters, clus_center)
  numMiC <<- nrow(microCenters)
},

remove_fromTemporary = function(){
  if(nrow(temporaryMemory) > 0){
    toDelete <- numeric()
    eTimestamp <- currentTime - temporaryMemory[,ncol(temporaryMemory)]
    toDelete <- which(eTimestamp > ts)
    if(!is.empty(toDelete) > 0)
      temporaryMemory <<- temporaryMemory[-toDelete,]
  }
},

check_recurrentClass = function(clus_center,data_point){
  recurrent <- FALSE
  if(length(sleepMemory)>=maxMiC){
    centerSleep <- get_sfmiccenters(sleepMemory)
    sleep_membership <- get_memberships(centerSleep,clus_center)
    sleep_maxMicroclusters <- find_maxClass(sleepMemory,sleep_membership)

    allclasses <- unique(get_allClasses(allSFMiC))
    existingClass <- which(allclasses == sleep_maxMicroclusters["class"])

    if(sleep_maxMicroclusters["membership"] > theta){
      if(!is.empty(existingClass)){
        create_SFMiC(data_point,sleep_maxMicroclusters["class"],currentTime)
        recurrent <- TRUE
      }
    }
  }
  return(recurrent)
},

check_PN = function(clus_memberships,clus_center,instance_cluster){
  criar <- TRUE
  #calculate dpf
  clus_dispersion <- calculate_dpf(clus_memberships,clus_center,temporaryMemory[,-ncol(temporaryMemory)])
  #calculate fr
  all_PNs <- get_PNs()
  FR <- calculate_pnFR(clus_center,clus_dispersion, all_PNs)
  if(!is.empty(FR)){
    maxFR <- max(FR)
    if(maxFR >= fr_threshold){
      criar <- FALSE
      index_max <- which.max(FR)
      maxClass <- allSFMiC[[all_PNs[index_max]]]$get_class()
      create_pnSFMiC(clus_center,maxClass,clus_dispersion["ssde"])
      PN_mics <- sameClass_sfmics(allSFMiC,maxClass)
      if(length(PN_mics) > maxMiC){
        remove_oldest(PN_mics)
      }
    }

  }
  if(criar){
    create_pnSFMiC(clus_center,paste0("PN",PN),clus_dispersion["ssde"])
    PN <<- PN + 1
  }
},

calculate_sil = function(cl,Xcal){
 dist_temporary <- as.matrix(dist(Xcal))
 sil <- silhouette(cl$cluster, dmatrix =  dist_temporary)
 classes <- c(1:k_temp)
 point_cluster <- sapply(classes, function(class){which(sil[,"cluster"]==class)})
 mean_sil <- sapply(point_cluster, function(points){mean(sil[points,"sil_width"])})
 mean_sil[is.nan(mean_sil)] <- 0
 return(mean_sil)
}


)


