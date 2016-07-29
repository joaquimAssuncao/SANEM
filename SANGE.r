#Sample data
x <- c(76,84,82,1,94,5,59,73,87,88,24,40,61,26,32,41,52,19,0,16,18,8,22,99, 42,85,3,11,27,30,31,47,95,13,68,79,75,35,53,43,29,10,57,97,92,38,65,81,2,49,20,72,17,64,63,71,100,50,60,36,67,44,51,80,28,70,58,46,12,83,86,90,6,66,62,91,9,14,55,25,78,96,54,39,48,34,56,93,21,37,45,33,23,4,98,77,89,7,74,15)
set.seed(9)
x <- floor(runif(100, min = 0, max = 100))
set.seed(55)
x2 <- floor(runif(100, min = 0, max = 80))
set.seed(90)
x3 <- floor(runif(100, min = 0, max = 90))
myDF <- data.frame(x,x2,x3)

#################### SANEM ##########################
# 2015 by Joaquim Assuncao (www.joaquim.pro.br)
# jassuncao@outlook.com
#####################################################

#################### SANEM parameters ##########################
# observations: The set of observations (N variables in a Data Frame) to be modeled
# SANname: The name of your file and network
# alphz: The number of characteres to PAA, also the number of states per automaton
# chron: The number of states to represent time
# intFunc: A vector of state to the integration functions.
# sync: The minimum value to create synchronizing events
# timeChron: Vector of taxes to state change in time (Chron automaton), if NULL random numbers are applied
# lseg: The level of segmentation divide the total number of observations (more=less dimensions & less accuracy)
#OUTPUT: .san code
#################### #################### #################### 

#SANEM(myDF,"mySAN", alphz=4,chron=3)
SANEM <<- function(observations, SANname="mySAN", alphz=4, chron=4, intFunc=NULL,
                  sync=0.7, timeChron=NULL,lseg=2) { 
     require(memisc)
     require(ggplot2)
     if ((intFunc=="A,A,A...") || (intFunc==NULL)){
          intFunc <- rep("A",alphz) 
     }  
     if (is.null(timeChron)){
          timeChron<-abs(rnorm(chron,mean=1/chron, sd=1/(chron+chron)))
     }
     #Create a basic .san file
     fileConn<-file(paste(SANname,".san",sep="")) #create a san file
     sink(fileConn)
     cat("identifiers",sep="\n")
     cat("tx_0 = 0.00001;")

     
     #SYMBOLIC DATA
     obs <- observations
     symb <- data.frame(matrix(ncol = ncol(obs), nrow = nrow(obs)/lseg))
     for (i in 1:ncol(obs)){
          symb[i] <- MakeSymbolicData(obs[i],alphasize=alphz,lseg=2)
     }
          # Makes .san partial reachability string
     pReach <- "( "
     for (c in 1:(ncol(symb))){
          pReach <- paste(pReach,"(st A",c," == ",symb[1,c],")", " && ", sep = "")
     }
     pReach <- paste(pReach,"(st Chron"," == ","A",")", sep = "")   
     pReach <-  paste(pReach,");")
     
     
     # CODER
          #Getting the Probability Matrix for each serie (TPM)

          #--create a list
          # //Gamma is a list of matrices//
     blocSize <- nrow(symb)/chron #to the correlation
     gamma <- list()
     for (c in 1:ncol(symb)){
          gamma[[c]] <- list(matrix(rep(0,alphz), nrow=alphz,ncol=alphz))
          gamma.new <- gamma
     }
          #--assign the probabilities
     for (c in 1:ncol(symb)){
          for (k in 1:nrow(symb)){
               tmp <- as.numeric(symb[,c])
               gam <- matrix(unlist(gamma[[c]]),ncol=alphz)
               gam[tmp[k],tmp[k+1]] <- gam[tmp[k],tmp[k+1]] +1#i?k
               gamma[[c]] <- gam
          }
     }
     
          #--normalize
     gam.new <- gam
     for (c in 1:ncol(symb)){
          for (k in 1:ncol(gam)){
               gam <- matrix(unlist(gamma[[c]]),ncol=alphz)
               s_gam <- sum(gam[k,])
               gam.new[k,] <- gam[k,]/s_gam
               if (is.na(gam.new[k,c])) { 
                    gam.new[k,] <- 1/nrow(gam)
               }
               gamma.new[[c]] <- gam.new
          }
     }
          
          #IDENTIFIERS
          #Taxes
     tax_id <- ""
     ev_id <- ""
     for (ts in 1:ncol(symb)){ #Number of Time Series
          tmp <- gamma.new[[ts]]
          for (l in 1:nrow(tmp)){ #Line to...
               for (c in 1:ncol(tmp)){ #Column
                    curValue <- tmp[l,c]
                    if (curValue > 0){
                         tax_id <- paste(tax_id,"\n","tx_ts",ts,"l",l,"c",c," = ", curValue,";",sep="")
                         ev_id <- paste(ev_id,"\n","loc ","e_ts",ts,"l",l,"c",c," (", 
                                        "tx_ts",ts,"l",l,"c",c,")",";", sep="")
                    }
               } 
          } 
     }
     #chron tx
     for (iChron in 1:chron){
          tax_id <- paste(tax_id,"\n","tx_timeChr",iChron,"_",iChron+1," = ", timeChron[iChron],";", sep="")
     }
     

     #EVENTS
     sync_id <- ""
     for (iChron in 1: chron){ #chron div.
          ini <- blocSize*(iChron-1)+1
          fin <- blocSize*(iChron)
          for (ts in 1:(ncol(symb)-1)){ #Foreach TS-1
               for (tsComp in (ts+1):ncol(symb)) {#against...
                    curCor <- cor(as.numeric(symb[ini:fin,ts]),as.numeric(symb[ini:fin,tsComp]))
                     if (curCor > sync){
                          #sync
                          sync_id <- paste(sync_id,"\n","tx_chr",iChron,"_ts",ts,"l",l,"c",c," = ", curCor,";",sep="")
                          ev_id <- paste(ev_id,"\n","syn ","e_chr",iChron,"_ts",ts,"l",l,"c",c," (", 
                                         "tx_chr",iChron,"_ts",ts,"l",l,"c",c,")",";", sep="")
                     }   
               }
          }
          ev_id <- paste(ev_id,"\n","syn ","e_timeChr",iChron,"_",iChron+1," (", 
                         "tx_timeChr",iChron,"_",iChron+1,")",";", sep="")    
     }
     cat(tax_id)
     cat(sync_id, sep="\n")
     cat("events",sep="\n")  
     cat(ev_id, sep="\n")
          
     #PARTIAL REACHABILITY
     cat(paste("partial reachability = ",pReach,sep = ""), sep="\n")
     
     #NETWORK
     network_id <-""
     for (ts in 1: ncol(symb)) { #ts or Automaton
          aPrinted <- FALSE;
          gam <- matrix(unlist(gamma.new[[ts]]),ncol=alphz)
          for (l in 1:alphz){ #
               stPrinted <- FALSE;
               for (c in 1:alphz){#
                    if (aPrinted==FALSE){
                         network_id <- paste(network_id,"\n","aut A",ts,sep="")
                         aPrinted<-TRUE
                    }
                    if (gam[l,c]>0){
                         if (stPrinted == FALSE){
                              network_id <- paste(network_id,"\n"," stt ",LETTERS[l]," to (",LETTERS[c],") ", "e_ts",ts,"l",l,"c",c, sep="")
                              stPrinted<-TRUE  
                         } else{
                              network_id <- paste(network_id,"\n","     to (",LETTERS[c],") ", "e_ts",ts,"l",l,"c",c, sep="")
                         }                               
                    }
                    
               }
          }
          
     }
     #Chronos internal network
     network_id <- paste(network_id,"\n","aut Chron ", sep="")
     for (i in 1:chron){
          network_id <- paste(network_id,"\n"," stt ",LETTERS[i]," to (", ifelse(i==chron,"A",LETTERS[i]),") ",
                              "e_timeChr",iChron,"_",iChron+1, sep="")
     }
     
     cat(paste("network ",SANname, "(continuous)"), sep="\n")
     cat(network_id,sep="\n")
     cat("results", sep="\n")
     
     #INTEGRATION FUNCTION -- Results
      function_id <- "My_integration_function = ("
      for (i in 1:(ncol(symb)-1)){
           function_id <-  paste(function_id,"\n","(st A",i,"==",intFunc[i],") && ", sep="")
      }
      function_id <-  paste(function_id,"\n","(st A",ncol(symb),"==",intFunc[ncol(symb)],")", sep="")
      function_id <-  paste(function_id,");")
      cat(function_id, sep="\n")
   
     sink()
     close(fileConn)
     
     return(scan(paste(SANname,".san",sep=""), what="character", sep="\n"))
}
#END SANEM






MakeSymbolicData <<- function(filename, serie='ALL', alphasize=10, lseg=2, haveheader = TRUE) {
     require(memisc)
     if (class(filename)=="data.frame"){
          myDF <- filename;
     }else {
          myDF <- read.csv(filename, haveheader);
     }
     if (serie=='ALL'){
          serie <- 1:ncol(myDF);
     }
     
     #---Must create TS and do every steps for each one
     for (i in (serie)){
          myVec <- myDF[,i];
          # lseg must be divisible by data length
          if ((length(myVec) %% lseg==0)){
               dim(myVec) <- c(lseg,length(myVec)/lseg) 
          }
          #Now, I shall take the mean of each column, let's change to my TS
          myTS <- colMeans(myVec);
          #cutSize determines where to put a symbol
          cutSize <- ((max(myVec)) - (min(myVec)))/alphasize;
          
          #Now, I map the TS according the cutSize
          symbolicData <- cases(
               "A" = myTS <= cutSize+min(myVec),
               "B" = myTS <= cutSize*2+min(myVec),
               "C" = myTS <= cutSize*3+min(myVec),
               "D" = myTS <= cutSize*4+min(myVec),
               "E" = myTS <= cutSize*5+min(myVec),
               "F" = myTS <= cutSize*6+min(myVec),
               "G" = myTS <= cutSize*7+min(myVec),
               "H" = myTS <= cutSize*8+min(myVec),
               "I" = myTS <= cutSize*9+min(myVec),
               "J" = myTS <= cutSize*10+min(myVec),
               "K" = myTS <= cutSize*11+min(myVec),
               "L" = myTS <= cutSize*12+min(myVec),
               "M" = myTS <= cutSize*13+min(myVec),
               "N" = myTS <= cutSize*14+min(myVec),
               "O" = myTS <= cutSize*15+min(myVec),
               "P" = myTS <= cutSize*16+min(myVec),
               "Q" = myTS <= cutSize*17+min(myVec),
               "R" = myTS <= cutSize*18+min(myVec),
               "S" = myTS <= cutSize*19+min(myVec),
               "T" = myTS <= cutSize*20+min(myVec)
          )
          #just to change the type
          symbolicData <- data.frame(symbolicData);
          return(symbolicData)
          #Then, print the output
          #print(symbolicData);
          #Must do something if there are multiple entries...  print("--- Time Serie  ----");
     }
}


#################### GetFittingsMeasurements ##########################
# 2015 by Joaquim Assuncao (www.joaquim.pro.br)
# jassuncao@outlook.com
#####################################################

#################### GetFittingsMeasurements parameters ##########################
# observations: The set of observations (N variables in a Data Frame) to be modeled
# ite: Number of iterations

#OUTPUT   list with 
     # mllk : The average Log likelihood
     # avgMatch : Vector with Averange number of matches in the exatly position
     # alldtwDist : Vector with Average DTW distance for the generated TS and the real data
#################### #################### #################### 
GetFittingsMeasurements <- function(observations, ite=10){
     #SYMBOLIC DATA
     obs <- observations
     symb <- data.frame(matrix(ncol = ncol(obs), nrow = nrow(obs)/lseg))
     for (i in 1:ncol(obs)){
          symb[i] <- MakeSymbolicData(obs[i],alphasize=alphz,lseg=2)
     }
     require(dtw)
     require(markovchain)
     llk <- NULL
     nTS <- ncol(symb)
     fp <- list()
     for (c in 1:nTS){
          mcFit <- markovchainFit(symb[,c])
          fp[[c]] <- mcFit$estimate[]
          llk[c] <- mcFit$logLikelihood
     }
     mllk <- mean(llk)
       
     alldtwDist <- rep(0,nTS)
     for (TS in 1:nTS){
          distance <- 0
          match <- 0
          dtwDist <- 0
          for (z in 1:ite){
               result <- NULL
               tProbs <- fp[[TS]]
               X <- as.numeric(symb[1:(nrow(obs)/3),TS])
               cumList <- unique(X)
               for( i in 1:(length(X))){ 
                    probVec <- tProbs[X[1],]
                    cProbVec <- cumsum(probVec)
                    rand <- runif(1)-0.01
                    whichVec <- which(rand < cProbVec)[1] 
                    result <- c(result,sample(cumList[[whichVec]],1))
                    tProbs <- fp[[TS]]%*%tProbs
               }
               
               #dist Measurement
               match <- match+sum(X==result)
               alignment <- dtw(result,X,step=asymmetricP1,keep=TRUE)
               dtwDist <- dtwDist+alignment$distance/length(X)
          }
          avgMatch[TS] <- match/ite/length(X)
          alldtwDist[TS] <- dtwDist/ite
     }

     return(list(mllk=-llk, avgMatch=avgMatch, dtwDist=alldtwDist)) 
}

     
