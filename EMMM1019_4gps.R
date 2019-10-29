
library(readr)
library(xlsx)
library(rmutil)
library(ggplot2)
library(geepack)

BoricAcidMousedata_processed <- read_csv("C:/Users/xinqi/Desktop/BoricAcidMousedata_processed.csv")
#BoricAcidRatdata_processed <- read_csv("C:/Users/xinqi/Desktop/BoricAcidRatdata_processed.csv")
#BoricAcidRatdata_processed <- read_csv("z:/EOGeorge/Data/Binary/BoricAcidRatdata_processed.csv")

#BBPdata_processed <- read_csv("C:/Users/xinqi/Desktop/BBPdata_processed.csv")
#DEHPmousedata_processed <- read_csv("C:/Users/xinqi/Desktop/DEHPmousedata_processed.csv")
#DEHPratdata_processed <- read_csv("C:/Users/xinqi/Desktop/DEHPratdata_processed.csv")
#DESdata_processed <- read_csv("C:/Users/xinqi/Desktop/DESdata_processed.csv")
#MEHPdata_processed <- read_csv("C:/Users/xinqi/Desktop/MEHPdata_processed.csv")

#BoricAcidRatdata_processed <- read_csv("C:/Users/Administrator/Desktop/BoricAcidRatdata_processed.csv")
#BBPdata_processed <- read_csv("C:/Users/Administrator/Desktop/BBPdata_processed.csv")
#BoricAcidMousedata_processed <- read_csv("C:/Users/Administrator/Desktop/BoricAcidMousedata_processed.csv")
#DEHPmousedata_processed <- read_csv("C:/Users/Administrator/Desktop/DEHPmousedata_processed.csv")
#DEHPratdata_processed <- read_csv("C:/Users/Administrator/Desktop/DEHPratdata_processed.csv")
#DESdata_processed <- read_csv("C:/Users/Administrator/Desktop/DESdata_processed.csv")
#MEHPdata_processed <- read_csv("C:/Users/Administrator/Desktop/MEHPdata_processed.csv")


# generate dataset
set.seed(123)

MaxClusterSize = 10
qbaseline = dbetabinom(0:MaxClusterSize,size=MaxClusterSize, m=0.5, s=1)
thetaPopulat = c(0.75,0.5,0.25)


DataGenerat = function(MarginalProb,THETAS,MaxClusize) {
  
  # retrieve lambda's from marginal probabilities 
  Lam0Base = rep(0, times=MaxClusize+1)
  for(kk in 0:MaxClusize) {
    jvec = seq(from=0,to=MaxClusize-kk,by=1)
    Lam0Base[kk+1] = sum(choose(MaxClusize-jvec,kk) * MarginalProb[MaxClusize-jvec+1]) / choose(MaxClusize,kk)
  }
  Lam1Resp = rep(0, times=MaxClusize+1)
  Lam2Resp = rep(0, times=MaxClusize+1)
  Lam3Resp = rep(0, times=MaxClusize+1)
  Lam1Resp[1] = 1
  Lam2Resp[1] = 1
  Lam3Resp[1] = 1
  for (jj in 1:MaxClusize) {
    Lam1Resp[jj+1] = Lam0Base[jj+1] * (THETAS[1]^jj)
    Lam2Resp[jj+1] = Lam0Base[jj+1] * (THETAS[2]^jj)
    Lam3Resp[jj+1] = Lam0Base[jj+1] * (THETAS[3]^jj)
  }
  
  # retrieve qt's using model estimate of lambda
  QfromLambda = function(clustersize,Lambdas) {
    QS = numeric(clustersize+1)
    for(s in 0:clustersize) {
      js = seq(from=0,to=(clustersize-s),by=1)
      QS[s+1] = choose(clustersize,s)*
        sum(ifelse(Lambdas[s+js+1]!=0,((-1)^js)*choose((clustersize-s),js)*Lambdas[s+js+1],0))
    }
    QS = ifelse(QS<0,0,QS)
    QS
  }

  # distribution of clustersize: uniform
  PERMUTATION = function(DOSELVL,Lambdas) {
    SUBSAMPLE = matrix(0,nrow=250,ncol=2)
    # clustersizes are uniformly distributed 
    ClusterSizes = seq(from=1,to=MaxClusize,by=1)
    SUBSAMPLE[,1] = sample(ClusterSizes,size=(dim(SUBSAMPLE)[1]),replace=TRUE,prob=rep(1,times=MaxClusize))
    for(ii in 1:(dim(SUBSAMPLE)[1])) {
      qts = QfromLambda(SUBSAMPLE[ii,1],Lambdas[1:(SUBSAMPLE[ii,1]+1)])
      multivar =  as.numeric(rmultinom(n=1,size=1,prob=qts))
      SUBSAMPLE[ii,2] = which(multivar>0)-1
    }
    SUBSAMPLE = as.data.frame(SUBSAMPLE)
    SUBSAMPLE = cbind.data.frame(DOSELVL,SUBSAMPLE,1)
  }
  
  sampleA = PERMUTATION(0,Lam0Base)
  sampleB = PERMUTATION(1,Lam1Resp)
  sampleC = PERMUTATION(2,Lam2Resp)
  sampleD = PERMUTATION(3,Lam3Resp)
  DATAS = rbind.data.frame(sampleA,sampleB,sampleC,sampleD)
  names(DATAS) = c("dose","nTotal","nResp","Freq")
  DATAS
}

GENERATE = DataGenerat(qbaseline,thetaPopulat,MaxClusterSize)


### Reverse coding of responses ###
ReverseCoding = function(DATASET) {
  DATASET$nResp =  DATASET$nTotal - DATASET$nResp
  DATASET
}
#BoricAcidRatdata_Reversed = ReverseCoding(BoricAcidRatdata_processed)
#BBPdata_Reversed = ReverseCoding(BBPdata_processed)
#BoricAcidMousedata_Reversed = ReverseCoding(BoricAcidMousedata_processed)
#DEHPmousedata_Reversed = ReverseCoding(DEHPmousedata_processed)
#DEHPratdata_Reversed = ReverseCoding(DEHPratdata_processed)
#DESdata_Reversed = ReverseCoding(DESdata_processed)
#MEHPdata_Reversed = ReverseCoding(MEHPdata_processed)



### Data dictionary ###
# dose:   boric acid concentration (%)
# nTotal: number of rats/mice/rabbits per group w.r.t a certain dose; 
#         exchangeable within cluster;
#         essentially the cluster size from real data;
# nResp:  number of resorptions, dead, or live fetuses examined at termination;
#         observation with one or multiple responses is recorded as "1";
# Freq:   frequency of a certain number of response (nResp) with offset (nTotal) w.r.t a certain dose
#         number of "success"

### Data manipulation ###
# Creating data set for analysis
DATA.Reversion = function(DataSets) {
  attach(DataSets)
  datareversion = function(columnname,reptimes) {
    repeatterm = rep(columnname,times=reptimes)
  }
  dataprocessed = cbind.data.frame(datareversion(dose,Freq),
                                   datareversion(nTotal,Freq),
                                   datareversion(nResp,Freq))
  names(dataprocessed) = c("dose","nTotal","nResp")
  detach(DataSets)
  dataprocessed
}



# doselvls are values of solution concentrations
DATA.Subset = function(DataProcessed,doselvls) {
  attach(DataProcessed)
  dosesubset = NULL
  # Clusters partitioning with respect to dose levels
  for (lvls in 1:length(doselvls)) {
    dosesubset = rbind.data.frame(dosesubset,
                                  subset(DataProcessed,dose==doselvls[lvls]))
  }
  names(dosesubset) = c("dose","nTotal","nResp")
  detach(DataProcessed)
  dosesubset
}



# returns a vector of solution concentrations
DoseLevels = function(Data.Analysis) {
  sort(unique(Data.Analysis$dose),decreasing=FALSE)
}



# Function to estimate qk
# where qk represents: for the complete data setting, given all cluster sizes are Rglobal,
# the probability of achieving k successes/responses
# doselevel are values of solution concentrations
qestimate = function(DATA.analysis,doselevel,RMAX,itermax=1000) {
  BARat = DATA.Subset(DATA.analysis,doselevel)
  ### Maximum Likelihood Estimation ###
  # Number of clusters within a doselevel
  K = dim(BARat)[1]
  # The kth cluster having rk observations (essentially cluster size)
  rk = BARat$nTotal
  meanclustersize = mean(rk)
  # Maximum cluster size within certain doselevel
  # Rglobal = max(rk)
  # Maximum cluster size within specified doselevels
  Rglobal = RMAX
  
  # Number of clusters with cluster size r (Kr)
  #r = sort(unique(rk)) # Domain of rk
  #Kr = as.numeric(xtabs(~nTotal, data=BARat))
  
  # Number of clusters of size r with exactly s reponses
  #rowname = sort(unique(rk))
  #colname = sort(unique(BARat$nResp)) # Domain of reponse s
  #nrs = xtabs(~nTotal+nResp, data=BARat)
  
  # For each cluster within doselevel, number of cluster size r & reponses s
  # Array "nrs.array" with dimension: {clustersize(r): 1,2,...,Rglobal}*{#response(s):0,1,...,Rglobal}*K
  # Matrix "nrs.matrix" with dimension: {clustersize(r): 1,2,...,Rglobal}*{#response(s):0,1,...,Rglobal}
  nrs.array = array(0,dim = c(Rglobal,(Rglobal+1),K))
  nrs.matrix = matrix(0,nrow=Rglobal,ncol=(Rglobal+1))
  for(i in 1:Rglobal) {
    for(ii in 0:Rglobal) {
      for(iii in 1:K) {
        nrs.array[i,ii+1,iii] = (BARat$nTotal[iii]==seq(1:Rglobal)[i] & 
                                   BARat$nResp[iii]==(seq(from=0,to=Rglobal,by=1)[ii+1]))
      }
      nrs.matrix[i,ii+1] = sum(nrs.array[i,ii+1,])
    }
  }
  #rownames(nrs.matrix) = seq(1:Rglobal)
  #colnames(nrs.matrix) = seq(from=0,to=Rglobal,by=1)
  
  ### With exchangeability assumption ###
  ## EM algorithm ##
  # Initial estimate of qs
  q0 = rep(1/(Rglobal+1), times = Rglobal+1)
  # Number of missing completely at random (MCAT) for each cluster size
  NMCAT = Rglobal - BARat$nTotal
  BARat = cbind.data.frame(BARat, NMCAT)
  
  # prst=Pr(T=tt|S=ss) for incompletedata clustersize=rr success "S=ss" & completedata clustersize=Rglobal success "T=tt"
  # Domain of tt: {0, 1, ..., Rglobal}; {s, s+1, s+2, ..., Rglobal-r+s}
  # Domain of rr: {1, 2, ..., Rglobal}; {1, 2, 3, ..., Rglobal}
  # Domain of ss: {0, 1, ..., Rglobal}; {max(0,t+r-Rglobal), ..., min(t,r)}
  # Lower & upper limits: values for "complete"; "imcomplete" situations
  # Scalar representation of "prst"
  prst = function(rr,ss,tt,qvec) {
    # prst is a scalar
    # qvec is a vector, rr/ss/tt are scalars
    if (rr==Rglobal) {
      prsTvec =ifelse(ss==tt,1,0)
    } else {
      numerator = choose(tt,ss) * choose((Rglobal-tt),(rr-ss)) * qvec[tt+1]
      tprime = seq(ss,Rglobal-rr+ss)
      denominator = sum(choose(tprime,ss) * choose(Rglobal-tprime,rr-ss) * qvec[tprime+1])
      prsTvec = numerator/denominator
    }
    prsTvec
  }
  # incorrect version
  #prst = function(rr,ss,tt,qvec) {
  #  # prst is a scalar
  #  # qvec is a vector, rr/ss/tt are scalars
  #  numerator = choose(tt,ss) * choose((Rglobal-tt),(rr-ss)) * qvec[tt+1]
  #  tprime = seq(ss,Rglobal-rr+ss)
  #  denominator = sum(choose(tprime,ss) * choose(Rglobal-tprime,rr-ss) * qvec[tprime+1])
  #  numerator/denominator
  #}
  
  # Vector representation of "prst" w.r.t "complete data" # of responses
  # Calculate prst for all values of T (completedata success) for fixed clustersize=rr incompletedata success "S=ss"
  # Domain of ss: {0, 1, ..., Rglobal}; {max(0,t+r-Rglobal), ..., min(t,r)}
  prst.vec = function(rr,ss,qvec) {
    # prst.vec is a vector
    # qvec is a vector, rr/ss/tt are scalars
    res = rep(0, length = Rglobal+1)
    # pRst=1 if ss=tt
    if (rr==Rglobal) { 
      res[ss+1] = 1 
    } else {
      tt = seq(ss, Rglobal-rr+ss)
      numerator = choose(tt,ss) * choose((Rglobal-tt),(rr-ss)) * qvec[tt+1]
      res[tt+1] = numerator / sum(numerator)
    }
    res
  }
  # incorrect version
  #prst.vec = function(rr,ss,qvec) {
  #  # prst.vec is a vector
  #  # qvec is a vector, rr/ss/tt are scalars
  #  res = rep(0, length = Rglobal+1)
  #  tt = seq(ss, Rglobal-rr+ss)
  #  numerator = choose(tt,ss) * choose((Rglobal-tt),(rr-ss)) * qvec[tt+1]
  #  res[tt+1] = numerator / sum(numerator)
  #  res
  #}
  
  ## Expectation Step
  # Scalar representation of nt(prime)=E[Nt|{Sk}]
  # Nt: #clusters with exactly t successes if the "missing data" were present,
  #                                    so that all clusters have common size Rglobal
  # Sk: observed data; Sk=Yk1+Yk2+Yk3+...+Ykr as "incomplete data"
  # Tk: "complete data"; Tk=Sk+ Yk(r+1) + Yk(r+2) + ... + YkR
  expectation = function(tt, qvec) {
    summ = 0
    for(rr in 1:Rglobal) {
      for(ss in max(0,tt+rr-Rglobal):min(tt,rr)) {
        if (nrs.matrix[rr,ss+1] != 0) {
          summ = summ + nrs.matrix[rr,ss+1] * prst(rr,ss,tt,qvec)
        }
      }
    }
    return(summ)
  }
  
  # Vector representation of nt(prime)=E[Nt|{Sk}] w.r.t "complete data" # of responses
  # Equivalent to the scalar representation: already proved via data
  expectation.vec = function(qvec) {
    p = array(0, dim=c(Rglobal,Rglobal+1,Rglobal+1))
    memory.matrix = array(0, dim=c(Rglobal,Rglobal+1,Rglobal+1))
    for (rr in 1:Rglobal){
      for (ss in 0:rr){
        p[rr,ss+1, ] = prst.vec(rr,ss,qvec)
      }
    }
    res = numeric(Rglobal+1)
    for (tt in 0:Rglobal) {
      for(rr in 1:Rglobal) {
        for(ss in max(0,tt+rr-Rglobal):min(tt,rr)) {
          if(nrs.matrix[rr,ss+1] != 0) {
            memory.matrix[rr,ss+1,tt+1] = nrs.matrix[rr,ss+1] * p[rr,ss+1,tt+1]
          }
        }
      }
      res[tt+1] = sum(memory.matrix[,,tt+1])
    }
    return(res)
  }
  
  
  ## Maximization step
  maximization = function(qvec) {
    ex = expectation(seq(from=0,to=Rglobal,by=1), qvec)
    new.qvec = ex/K
    new.qvec
  }
  
  # Iteration continues until convergence satisfied or itermax reached
  iter = 0
  difference = 10
  q0 = rep(1/(Rglobal+1),times=Rglobal+1)
  
  # Absolute convergence: difference <= 1e-8
  while((difference > (1e-8)) && (iter < itermax)) {
    iter = iter + 1
    qnext = maximization(q0)
    difference = sum(abs(qnext-q0))
    q0 = qnext
  }
  # cat("The estimated q0 is",q0, "\n")
  
  # Calculation of p(r)
  sequence = seq(0:Rglobal)
  probability = ifelse(q0!=0,
                       q0/choose(Rglobal,sequence),0) 
  # cat("The estimated pr is",probability, "\n")
  
  # Calculation of lambda(r)
  lambdaest = rep(0, times=Rglobal+1)
  for(kk in 0:Rglobal) {
    jvec = seq(from=0,to=Rglobal-kk,by=1)
    lambdaest[kk+1] = sum(choose(Rglobal-jvec,kk)*q0[Rglobal-jvec+1]) / choose(Rglobal,kk)
  }
  # cat("The estimated lambda is",lambdaest, "\n")
  
  # Calculation of ntprime
  ntprime = expectation.vec(q0)
  # Handling NaN's in the ntprime
  ntprime[is.na(ntprime)] = 0
  # ntprime = c(ntprime[!is.na(ntprime)],
  #            rep(0,times=(length(ntprime)-sum(ntprime>0,na.rm=TRUE))))
  
  list(q=q0,p=probability,lambda=lambdaest,ntprime=ntprime,nrs=nrs.matrix,
       meanclustersize = meanclustersize)
}
   
QestGenerate1 = qestimate(GENERATE,1,MaxClusterSize)
QestGenerate2 = qestimate(GENERATE,2,MaxClusterSize)
QestGenerate3 = qestimate(GENERATE,3,MaxClusterSize)
QestGenerate0 = qestimate(GENERATE,0,MaxClusterSize)

plot(QestGenerate1$lambda,col="black",ylim=c(0,1))
lines(QestGenerate2$lambda,col="red")
lines(QestGenerate3$lambda,col="green")
lines(QestGenerate0$lambda,col="blue")

plot(QestGenerate1$q,col="black",ylim=c(0,0.35))
lines(QestGenerate2$q,col="red")
lines(QestGenerate3$q,col="green")
lines(QestGenerate0$q,col="blue")


# inputs: Data.Analysis & doselvlsinput
# notice the order of input doselevels
# doselvlsinput: integers
# 1st is baseline A, 2nd is response B, 3rd is response C and so on 
modelestimate = function(Data.Analysis,doselvlsinput,SIMULATORNOT,RMAX) {
  # extracting the dataset corresponding to particular doselvls
  modeldat = DATA.Subset(Data.Analysis,doselvlsinput)
  
  # Maximum cluster size within the extracted doselevels
  # Rglobal = max(modeldat$nTotal)
  # if dataset is from simulation dataset, then Rglobal=original Rglobal
  # if dataset is from real observations, then Rglobal=max(modeldat$nTotal)
  Rglobal = ifelse(SIMULATORNOT==TRUE,RMAX,max(modeldat$nTotal))
  
  # assuming uniform effectiveness within two specifc doselevels
  # initial value of lambda for the baseline
  # baseline, higher values in this situation 
  #          resultpooled = qestimate(Data.Analysis,doselvlsinput,Rglobal)
  #          lambdaApooled = resultpooled$lambda
  # reason that qtA also pooled over all doselevels: lambdaA never used in EM algorithm
  #          qtApooled = resultpooled$q          # baseline
  
  # notice the order of input doselevels
  resultA = qestimate(Data.Analysis,doselvlsinput[1],Rglobal)  # baseline A
  resultB = qestimate(Data.Analysis,doselvlsinput[2],Rglobal)  # response B
  resultC = qestimate(Data.Analysis,doselvlsinput[3],Rglobal)  # response C
  resultD = qestimate(Data.Analysis,doselvlsinput[4],Rglobal)  # response D
  
  meanclustersizeA = resultA$meanclustersize
  meanclustersizeB = resultB$meanclustersize
  meanclustersizeC = resultC$meanclustersize
  meanclustersizeD = resultD$meanclustersize
  
  lambdaAnotpooled = resultA$lambda  # baseline A, nonparametric estimate for comparison
  lambdaB = resultB$lambda           # reponse B, nonparametric estimate for comparison
  lambdaC = resultC$lambda           # reponse C, nonparametric estimate for comparison
  lambdaD = resultD$lambda           # reponse D, nonparametric estimate for comparison
  
  # Initial values for qt's from group A & B
  qtA.notpooled = resultA$q          # baseline A, nonparametric estimate for comparison
  qtA.vec = qtA.notpooled            # baseline qtA.notpooled  # qtApooled before
  qtB.vec = resultB$q                # response B
  qtC.vec = resultC$q                # response C
  qtD.vec = resultD$q                # response D
  
  nrsA = resultA$nrs                 # baseline A
  nrsB = resultB$nrs                 # response B
  nrsC = resultC$nrs                 # response C
  nrsD = resultD$nrs                 # response D
  
  ###########################################################
  ##################### EM MM ALGORITHM #####################
  ######## Expectation Maximization Minorize-Maximize #######
  ###########################################################
  
  # Initial value for theta's
  theta1 = lambdaB[2] / lambdaAnotpooled[2]
  theta2 = lambdaC[2] / lambdaAnotpooled[2]
  theta3 = lambdaD[2] / lambdaAnotpooled[2]
  
  # Relation between qtB.vec and qtA.vec to obtain qtB
  qBqATransform <- function(thetas,qqttAA) {
    qqttBB = numeric(Rglobal+1)
    for(t in 0:Rglobal) {
      ssuumm = 0
      for(alpha in t:Rglobal) {
        ssuumm = ssuumm + 
          (thetas^t)*((1-thetas)^(alpha-t))*choose(alpha,t)*qqttAA[alpha+1]
      }
      qqttBB[t+1] = ssuumm
    }
    qqttBB
  }
  # response B/C/D initial value
  qtB.transformed = qBqATransform(theta1,qtA.vec) 
  qtC.transformed = qBqATransform(theta2,qtA.vec)
  qtD.transformed = qBqATransform(theta3,qtA.vec)
  
  ######################
  ## EXPECTATION STEP ##
  ######################
  # values and parameters (initial values) required in E Step
  # nrs for group A & B & C & D
  # qtA for group A
  
  # functions needed in the E step
  # qBqATransform: obtain qtB, qtC according to SEMIPARAMETRIC using qtA & theta's
  # prst.vec: obtain prst, which will be used in expectation
  # ntprime.est: expectation
  # Qfunction: expectation of log-likelihood
  
  # Vector representation of "prst" w.r.t "complete data" # of responses
  # Calculate prst for all values of T (completedata success) for fixed clustersize=rr incompletedata success "S=ss"
  # Domain of ss: {0, 1, ..., Rglobal}; {max(0,t+r-Rglobal), ..., min(t,r)}
  prst.vec = function(rr,ss,qvec) {
    # prst.vec is a vector
    # qvec is a vector, rr/ss/tt are scalars
    res = rep(0, length = Rglobal+1)
    # pRst=1 if ss=tt
    if (rr==Rglobal) { 
      res[ss+1] = 1 
    } else {
      tt = seq(ss, Rglobal-rr+ss)
      numerator = choose(tt,ss) * choose((Rglobal-tt),(rr-ss)) * qvec[tt+1]
      res[tt+1] = numerator / sum(numerator)
    }
    res
  }
  # incorrect version
  #prst.vec = function(rr,ss,qvec) {
  #  # prst.vec is a vector
  #  # qvec is a vector, rr/ss/tt are scalars
  #  res = rep(0, length = Rglobal+1)
  #  tt = seq(ss, Rglobal-rr+ss)
  #  numerator = choose(tt,ss) * choose((Rglobal-tt),(rr-ss)) * qvec[tt+1]
  #  res[tt+1] = numerator / sum(numerator)
  #  res
  #}
  
  # Calculate expectation of Nt for different groups (different qt's & nrs's)
  # Nt: # of clusters with exactly t successes if the "missing data" were present,
  #     so that all clusters have common size Rglobal
  # Vector representation of nt(prime)=E[Nt|{Sk}] w.r.t "complete data" # of responses
  ntprime.est = function(qvec,nnrrss) {
    p = array(0, dim=c(Rglobal,Rglobal+1,Rglobal+1))
    memory.matrix = array(0, dim=c(Rglobal,Rglobal+1,Rglobal+1))
    for (rr in 1:Rglobal){
      for (ss in 0:rr){
        p[rr,ss+1, ] = prst.vec(rr,ss,qvec)
      }
    }
    res = numeric(Rglobal+1)
    for (tt in 0:Rglobal) {
      for(rr in 1:Rglobal) {
        for(ss in max(0,tt+rr-Rglobal):min(tt,rr)) {
          # Statement "if(nnrrss[rr,ss+1] != 0)" to prevent occurance of NAN, frequently used later on
          if(nnrrss[rr,ss+1] != 0) {
            memory.matrix[rr,ss+1,tt+1] = nnrrss[rr,ss+1] * p[rr,ss+1,tt+1]
          }
        }
      }
      res[tt+1] = sum(memory.matrix[,,tt+1])
    }
    return(res)
  }
  
  # Calculate ntprime's for group A & B & C at the k's iteration
  ntprimeAk = ntprime.est(qtA.vec,nrsA)
  # which is essentially the original values of ntprimeA = resultA$ntprime
  
  ntprimeBk = ntprime.est(qtB.transformed,nrsB)
  ntprimeCk = ntprime.est(qtC.transformed,nrsC)
  ntprimeDk = ntprime.est(qtD.transformed,nrsD)
  # different from ntprimeB = resultB$ntprime since qtB.transformed is used
  
  
  ### Q function ###
  # component1 represents qvec
  # component2 represents ntprime
  ntlog = function(component1,component2) {
    ifelse((component2>0 & component1>0),log(component1)*component2,0)
  }
  Qfunction = sum(ntlog(qtA.vec,ntprimeAk) + ntlog(qtB.transformed,ntprimeBk) +
                    ntlog(qtC.transformed,ntprimeCk) + ntlog(qtD.transformed,ntprimeDk))
  
  ##############################
  ##### MAXIMIZATION STEP ######
  # Minorize-Maximize inserted #
  ##############################
  
  # weight / coefficient
  weight = function(TT,THETA,QVEC) {
    # weight returns a vector
    # TT,THETA is a scalar
    # QVEC is a vector
    coeff = numeric(Rglobal+1)
    ALPHA = seq(from=TT,to=Rglobal,by=1)
    numerators = 
      choose(ALPHA,TT) * (THETA^TT) * ((1-THETA)^(ALPHA-TT)) * QVEC[ALPHA+1]
    coeff[ALPHA+1] = ifelse(numerators>0,(numerators/sum(numerators)),0) 
    coeff
  }
  
  # update of theta using EM MM algorithm
  thetaupdate1 = function(THETA,QVEC,NTPRIMEBK) {
    weightedsum1 = 0 ; weightedsum2 = 0
    weightmatrix = matrix(0,nrow=(Rglobal+1),ncol=(Rglobal+1))
    for(TT in 0:Rglobal) {
      weightmatrix[TT+1,] = weight(TT,THETA,QVEC)
    }
    for(TT in 0:Rglobal) {
      alphaR = numeric(Rglobal+1)
      alphat = seq(from=TT,to=Rglobal,by=1)
      alphaR[alphat+1] = alphat
      weightedsum1 = weightedsum1 + NTPRIMEBK[TT+1]*TT
      weightedsum2 = weightedsum2 + NTPRIMEBK[TT+1]*sum(weightmatrix[TT+1,]*alphaR)
    }
    weightedsum1/weightedsum2
  }
  
  # update of theta using EM MM algorithm
  thetaupdate2 = function(THETA,QVEC,NTPRIMEBK) {
    weightedsum1 = 0 ; weightedsum2 = 0
    for(TT in 0:Rglobal) {
      alphaR = numeric(Rglobal+1)
      alphat = seq(from=TT,to=Rglobal,by=1)
      alphaR[alphat+1] = alphat
      weightedsum1 = weightedsum1 + NTPRIMEBK[TT+1]*TT
      weightedsum2 = weightedsum2 + NTPRIMEBK[TT+1]*sum(weight(TT,THETA,QVEC)*alphaR)
    }
    weightedsum1/weightedsum2
  }
  # theta,qtA.vec: from the pre kth iteration
  # ntprimeBk: expectaion from the pre kth iteration
  theta1kplus1 = thetaupdate1(theta1,qtA.vec,ntprimeBk)
  # ntprimeCk: expectaion from the pre kth iteration
  theta2kplus1 = thetaupdate1(theta2,qtA.vec,ntprimeCk)
  # ntprimeDk: expectaion from the pre kth iteration
  theta3kplus1 = thetaupdate1(theta3,qtA.vec,ntprimeDk)
  
  # update of qtA's using EM MM algorithm
  qtAupdate = function(THETA1,THETA2,THETA3,QVEC,NTPRIMEAK,NTPRIMEBK,NTPRIMECK,NTPRIMEDK) {
    NUMERATOR = numeric(Rglobal+1)
    for(BETA in 0:Rglobal) {
      rho = seq(from=0,to=BETA,by=1)
      weight1BETA = numeric(BETA+1)
      weight2BETA = numeric(BETA+1)
      weight3BETA = numeric(BETA+1)
      for(RHO in 0:BETA) {
        weight1BETA[RHO+1] = weight(RHO,THETA1,QVEC)[BETA+1]
        weight2BETA[RHO+1] = weight(RHO,THETA2,QVEC)[BETA+1]
        weight3BETA[RHO+1] = weight(RHO,THETA3,QVEC)[BETA+1]
      }
      NUMERATOR[BETA+1] = NTPRIMEAK[BETA+1] + sum(weight1BETA*(NTPRIMEBK[1:(BETA+1)])) +
        sum(weight2BETA*(NTPRIMECK[1:(BETA+1)])) +
        sum(weight3BETA*(NTPRIMEDK[1:(BETA+1)]))
    }
    NUMERATOR / sum(NUMERATOR)
  }
  
  # theta1&2,qtA.vec: from pre k iteration
  # ntprimeAk,ntprimeBk,ntprimeCk: expectaion from pre k iteration
  qtAkplus1 = qtAupdate(theta1,theta2,theta3,qtA.vec,ntprimeAk,ntprimeBk,ntprimeCk,ntprimeDk)
  
  ## Iterations of EM algorithm until convergence ##
  # Absolute convergence: difference1 & difference1 <= 1e-7
  itermax = 5000; iter = 0
  # Initial value for qtA and theta's
  theta1 = lambdaB[2] / lambdaAnotpooled[2]
  theta2 = lambdaC[2] / lambdaAnotpooled[2]
  theta3 = lambdaD[2] / lambdaAnotpooled[2]
  qtAs = qtA.vec
  difference = 10
  
  while((difference > (1e-7)) & (iter < itermax)) {
    iter = iter + 1
    theta1pre = theta1
    theta2pre = theta2
    theta3pre = theta3
    qtAspre = qtAs
    qtB.transformedpre = qBqATransform(theta1pre,qtAspre) 
    qtC.transformedpre = qBqATransform(theta2pre,qtAspre) 
    qtD.transformedpre = qBqATransform(theta3pre,qtAspre) 
    
    # expectation 
    ntprimeAkpre = ntprime.est(qtAspre,nrsA)  
    ntprimeBkpre = ntprime.est(qtB.transformedpre,nrsB)  
    ntprimeCkpre = ntprime.est(qtC.transformedpre,nrsC)  
    ntprimeDkpre = ntprime.est(qtD.transformedpre,nrsD)  
    
    # maximization (MM inserted)
    # theta1&2pre,qtAspre: from previous kth iteration
    # ntprimeAkpre,ntprimeBkpre,ntprimeCkpre: expectaion from previous kth iteration
    theta1 = thetaupdate1(theta1pre,qtAspre,ntprimeBkpre)
    theta2 = thetaupdate1(theta2pre,qtAspre,ntprimeCkpre)
    theta3 = thetaupdate1(theta3pre,qtAspre,ntprimeDkpre)
    qtAs = qtAupdate(theta1pre,theta2pre,theta3pre,
                     qtAspre,ntprimeAkpre,ntprimeBkpre,ntprimeCkpre,ntprimeDkpre)
    difference = sum(abs(qtAs-qtAspre)) + abs(theta1-theta1pre) + 
      abs(theta2-theta2pre) + abs(theta3-theta3pre)
  }
  #cat("From the EM MM Algorithm, the estimated theta1 from the model is",theta1,sep=" ","\n")
  #cat("From the EM MM Algorithm, the estimated theta2 from the model is",theta2,sep=" ","\n")
  #cat("From the EM MM Algorithm, the estimated qtAs from the model is",qtAs,sep=" ","\n")
  plot(qtAs,main="EM MM Algorithm Estimate of qtAs (Baseline)",xlab="r",ylim=c(0,1))
  
  # using optimized value of qtA to obtain corresponding lambda's for baseline (group A)
  lambdabaseline = numeric(Rglobal+1)
  for(tt in 0:Rglobal) {
    sumttjj = 0
    for(jj in 0:(Rglobal-tt)) {
      sumttjj = sumttjj + choose((Rglobal-tt),jj)*qtAs[Rglobal-jj+1]/choose(Rglobal,(Rglobal-jj))
    }
    lambdabaseline[tt+1] = sumttjj
  }
  #cat("lambdaA from model is",lambdabaseline,sep=" ","\n")
  plot(lambdabaseline,main="Lambda Comparison Semi-parametric vs Non-parametric",ylim=c(0,1),type="b")
  points(lambdaAnotpooled,col="orange",type="b")
  
  # using qtA and theta from EM optimization results to obtain lambda's for response (group B)
  lambdaresp1 = numeric(Rglobal+1)
  lambdaresp2 = numeric(Rglobal+1)
  lambdaresp3 = numeric(Rglobal+1)
  exponents = seq(from=0,to=Rglobal,by=1)
  for(tt in 0:Rglobal) {
    lambdaresp1[tt+1] = (theta1^exponents[tt+1]) * lambdabaseline[tt+1]
    lambdaresp2[tt+1] = (theta2^exponents[tt+1]) * lambdabaseline[tt+1]
    lambdaresp3[tt+1] = (theta3^exponents[tt+1]) * lambdabaseline[tt+1]
  }
  #cat("lambdaB from model is",lambdaresp1,sep=" ","\n")
  #cat("lambdaC from model is",lambdaresp2,sep=" ","\n")
  
  # the overlay plot red vs blue is comparison of model estimate with non-parametric estimate
  points(lambdaresp1,col="red",type="b")
  points(lambdaB,col="blue",type="b")
  
  points(lambdaresp2,col="green",type="b")
  points(lambdaC,col="purple",type="b")
  
  points(lambdaresp3,col="yellow",type="b")
  points(lambdaD,col="brown",type="b")
  
  legend(5,1,legend=c("lambdaA from model","lambdaA nonparametric",
                      "lambdaB from model","lambdaB nonparametric",
                      "lambdaC from model","lambdaC nonparametric",
                      "lambdaD from model","lambdaD nonparametric"),
         col=c("black","orange","red","blue","green","purple","yellow","brown"),lty=1:8)
  
  # ggplot2 version
  #lambdacomp = cbind.data.frame(rep(c("0.4%","0.2%","0.1%","0%"),each=(2*17)),
  #                              rep(c(rep("SemiParametric",times=17),rep("NonParametric",times=17)),4),
  #                              rep((0:Rglobal),4),
  #                              c(lambdabaseline,lambdaAnotpooled,lambdaresp1,lambdaB,
  #                                lambdaresp2,lambdaC,lambdaresp3,lambdaD))
  #names(lambdacomp) = c("Doselevel","Model","Lambda","Estimates")
  
  #lambdaSemi = lambdacomp[lambdacomp$Model=="SemiParametric",]
  #lambdaNon = lambdacomp[lambdacomp$Model=="NonParametric",]
   
  #Semi = ggplot(lambdaSemi, aes(x=Lambda, y=Estimates, group=Doselevel)) +
  #  geom_line(aes(color=Doselevel)) +
  #  geom_point(aes(color=Doselevel)) +
  #  ggtitle("Lambda Comparison Semi-parametric vs Non-parametric") +
  #  theme(legend.position="bottom")
  
  #Combined = Semi + 
  #  geom_line(data=lambdaNon, aes(x=Lambda, y=Estimates, group=Doselevel), linetype="dashed") +
  #  geom_point(aes(color=Doselevel))

  qtBs = qBqATransform(theta1,qtAs)
  qtCs = qBqATransform(theta2,qtAs)
  qtDs = qBqATransform(theta3,qtAs)
  # also plot for qt comparison EM MM vs Nonparametric
  plot(qtAs,main="Qt Comparison EM MM vs Nonparametric",ylim=c(0,1),type="b")
  points(resultA$q,col="orange",type="b")
  
  points(qtBs,col="red",type="b")
  points(resultB$q,col="blue",type="b")
  
  points(qtCs,col="green",type="b")
  points(resultC$q,col="purple",type="b")
  
  points(qtDs,col="yellow",type="b")
  points(resultD$q,col="brown",type="b")
  
  legend(5,1,legend=c("qtA from model","qtA nonparametric",
                      "qtB from model","qtB nonparametric",
                      "qtC from model","qtC nonparametric",
                      "qtD from model","qtD nonparametric"),
         col=c("black","orange","red","blue","green","purple","yellow","brown"),lty=1:8)
  
  list(Rglobal=Rglobal,
       meanclustersizeA = meanclustersizeA,
       meanclustersizeB = meanclustersizeB,
       meanclustersizeC = meanclustersizeC,
       meanclustersizeD = meanclustersizeD,
       
       nrsA=nrsA,nrsB=nrsB,nrsC=nrsC,nrsD=nrsD,
       
       ntpriAmodel=ntprime.est(qtAs,nrsA),
       ntpriBmodel=ntprime.est(qtBs,nrsB),
       ntpriCmodel=ntprime.est(qtCs,nrsC),
       ntpriDmodel=ntprime.est(qtDs,nrsD),
       
       ntpriAnonpar=resultA$ntprime,
       ntpriBnonpar=resultB$ntprime,
       ntpriCnonpar=resultC$ntprime,
       ntpriDnonpar=resultD$ntprime,
       
       qtAnonpar=resultA$q,
       qtBnonpar=resultB$q,
       qtCnonpar=resultC$q,
       qtDnonpar=resultD$q,
       
       qtAmodel=qtAs,qtBmodel=qtBs,qtCmodel=qtCs,qtDmodel=qtDs,
       
       theta1model=theta1,
       theta2model=theta2,
       theta3model=theta3,
       
       lambdaAmodel=lambdabaseline,
       lambdaBmodel=lambdaresp1,
       lambdaCmodel=lambdaresp2,
       lambdaDmodel=lambdaresp3,
       
       lambdaAnonpar=lambdaAnotpooled,
       lambdaBnonpar=lambdaB,
       lambdaCnonpar=lambdaC,
       lambdaDnonpar=lambdaD)
}



##################################################
###### EM MM Algorithm Variance Estimation #######
######## Using the Observed Information ##########
##################################################

# THETA123: interger indicating dose level category
# LAMBDAalpha: interger indicating order in baseline lambda estimates
VARIANCE = function(ESTIMATES) {
  Rglobal = ESTIMATES$Rglobal
  thetasmodel = c(ESTIMATES$theta1model,ESTIMATES$theta2model,ESTIMATES$theta3model)
  lambdaAmodel = ESTIMATES$lambdaAmodel
  nrsA = ESTIMATES$nrsA
  nrsB = ESTIMATES$nrsB
  nrsC = ESTIMATES$nrsC
  nrsD = ESTIMATES$nrsD
  
  p1 = function(ALPHA,RR,SS,THETA123) {
    P11 = numeric(RR-SS+1)
    P12 = numeric(RR-SS+1)
    P13 = numeric(RR-SS+1)
    P14 = numeric(RR-SS+1)
    P15 = numeric(RR-SS+1)
    for(j in 0:(RR-SS)) {
      P11[j+1] = ((-1)^j) * choose(RR-SS,j) * (lambdaAmodel[SS+j+1]) * ((thetasmodel[THETA123])^(SS+j))
      P12[j+1] = ((-1)^j) *(SS+j)*(SS+j-1)*choose(RR-SS,j)*(lambdaAmodel[SS+j+1])*((thetasmodel[THETA123])^(SS+j-2))
      P13[j+1] = ((-1)^j) *(SS+j)*choose(RR-SS,j)*(lambdaAmodel[SS+j+1])*((thetasmodel[THETA123])^(SS+j-1))
      P14[j+1] = ((-1)^j) *(ALPHA-SS-j)* choose(RR-SS,j) * (lambdaAmodel[SS+j+1]) * ((thetasmodel[THETA123])^j)
      P15[j+1] = ((-1)^j) * choose(RR-SS,j) * (lambdaAmodel[SS+j+1]) 
    }
    P11sum = sum(P11)
    P12sum = sum(P12)
    P13sum = sum(P13)
    P14sum = sum(P14)
    P15sum = sum(P15)
    list(P11sum=P11sum, P12sum=P12sum, P13sum=P13sum,P14sum=P14sum,P15sum=P15sum)
  }
  
  p2 = function(THETA123) {
    if (THETA123 == 0) {
      NRS = nrsA
    } else if (THETA123 == 1) {
      NRS = nrsB
    } else if (THETA123 == 2) {
      NRS = nrsC
    } else {
      NRS = nrsD
    }
    NRS
  }
  
  # returns a scalar
  Dthetatheta = function(THETA123) {
    part2 = p2(THETA123)
    THETAgSquare = matrix(0,nrow=Rglobal,ncol=Rglobal+1)
    
    for(rr in 1:Rglobal) {
      for(ss in 0:rr) {
        part1 = p1(2,rr,ss,THETA123)
        part11 = part1$P11sum
        part12 = part1$P12sum
        part13 = part1$P13sum
        
        THETAgSquare[rr,ss+1] = part2[rr,ss+1] * (part12*part11-(part13^2)) / (part11^2)
      }
    }
    THETAgSquaresum = sum(THETAgSquare)
    THETAgSquaresum
  }
  
  # resturns a scalar
  Dthetalambda = function(THETA123,LAMBDAalpha) {
    part2 = p2(THETA123)
    ThetaLambda = matrix(0,nrow=Rglobal,ncol=Rglobal+1)
    
    for(rr in LAMBDAalpha:Rglobal) {
      for(ss in 0:min(LAMBDAalpha,rr)) {
        part1 = p1(LAMBDAalpha,rr,ss,THETA123)
        part11 = part1$P11sum
        part14 = part1$P14sum
        
        ThetaLambda[rr,ss+1] = part2[rr,ss+1] * part14 * ((-1)^(LAMBDAalpha-ss)) *
          choose(rr-ss,LAMBDAalpha-ss) * ((thetasmodel[THETA123])^(LAMBDAalpha+ss-1)) / (part11^2)
      }
    }
    ThetaLambdasum = sum(ThetaLambda)
    ThetaLambdasum
  }
  
  p3 = function(THETA123,ALPHA,RR,SS,NRSG) {
    P31 = NRSG[RR,SS+1] * ((-1)^(2*(ALPHA-SS)+1)) * (choose(RR-SS,ALPHA-SS)^2)
    P32 = NRSG[RR,SS+1] * ((-1)^(2*(ALPHA-SS)+1)) * (choose(RR-SS,ALPHA-SS)^2) * ((thetasmodel[THETA123])^(2*ALPHA))
    list(P31=P31,P32=P32)
  }
  
  Dlambda1lambda1 = function(LAMBDAalpha) {
    LambdaSquare1 = matrix(0,nrow=Rglobal,ncol=Rglobal+1)
    LambdaSquare2 = array(0,dim=c(Rglobal,Rglobal+1,3))
    
    for(rr in LAMBDAalpha:Rglobal) {
      for(ss in 0:min(rr,LAMBDAalpha)) {
        for(gg in 1:3) {
          part1 = p1(LAMBDAalpha,rr,ss,gg)
          part11 = part1$P11sum
          part15 = part1$P15sum
          part2 = p2(gg)
          part3 = p3(gg, LAMBDAalpha, rr, ss, part2)
          part32 = part3$P32
          
          LambdaSquare2[rr,ss+1,gg] = part32 / (part11^2)
        }
        part30 = p3(1, LAMBDAalpha, rr, ss, nrsA)$P31
        LambdaSquare1[rr,ss+1] = part30 / (part15^2)
      }
    }
    LambdaSquaresum = sum(LambdaSquare1) + sum(LambdaSquare2)
    LambdaSquaresum
  }
  
  p4 = function(THETA123,ALPHA,BETA,RR,SS,NRSG) {
    P41 = NRSG[RR,SS+1] * ((-1)^(ALPHA+BETA-2*SS+1)) *choose(RR-SS,ALPHA-SS)*choose(RR-SS,BETA-SS)
    P42 = NRSG[RR,SS+1] * ((-1)^(ALPHA+BETA-2*SS+1)) *choose(RR-SS,ALPHA-SS)*choose(RR-SS,BETA-SS) * ((thetasmodel[THETA123])^(ALPHA+BETA))
    list(P41=P41,P42=P42)
  }
  
  Dlambda1lambda2 = function(LAMBDAalpha,LAMBDAbeta) {
    LambdaSquare3 = matrix(0,nrow=Rglobal,ncol=Rglobal+1)
    LambdaSquare4 = array(0,dim=c(Rglobal,Rglobal+1,3))
    
    for(rr in max(LAMBDAalpha,LAMBDAbeta):Rglobal) {
      for(ss in 0:min(rr,LAMBDAalpha,LAMBDAbeta)) {
        for(gg in 1:3) {
          part1 = p1(LAMBDAalpha,rr,ss,gg)
          part11 = part1$P11sum
          part15 = part1$P15sum
          part2 = p2(gg)
          part4 = p4(gg, LAMBDAalpha, LAMBDAbeta, rr, ss, part2)
          part42 = part4$P42
          
          LambdaSquare4[rr,ss+1,gg] = part42 / (part11^2)
        }
        part40 = p4(1, LAMBDAalpha, LAMBDAbeta, rr, ss, nrsA)$P41
        LambdaSquare3[rr,ss+1] = part40 / (part15^2)
      }
    }
    Lambda1lambda2sum = sum(LambdaSquare3) + sum(LambdaSquare4)
    Lambda1lambda2sum
  }
  
  covariance = matrix(0,nrow=(3+Rglobal),ncol=(3+Rglobal))
  for (i in 1:3) {
    covariance[i,i] = Dthetatheta(i)
  }
  
  for (i in 1:3) {
    for (j in 1:Rglobal) {
      entry = Dthetalambda(i,j)
      covariance[i,j+3] = entry
      covariance[j+3,i] = entry
    }
  }
  
  for (j in 1:Rglobal) {
    covariance[j+3,j+3] = Dlambda1lambda1(j)
  }
  
  for (j in 1:Rglobal) {
    for (k in 1:Rglobal) {
      covariance[j+3,k+3] = Dlambda1lambda2(j,k)
    } 
  }
  covariance
}

# for exmaple, the simulated dataset
dataprocessed = DATA.Reversion(GENERATE)
doselevels = DoseLevels(GENERATE)
cat("The doselevels in this dataset are",doselevels,sep=" ","\n")

# first loc is baseline A; second loc is response B; third loc is response C
Doselvlsinput = doselevels[c(1,2,3,4)]
cat("The input doselevels to the model are",Doselvlsinput,sep=" ","\n")

modelresult34 = modelestimate(dataprocessed,Doselvlsinput,SIMULATORNOT=FALSE,RMAX=1)
cat("The model estimate results are:",sep=" ","\n")
print(modelresult34)

covBAmouse = VARIANCE(modelresult34)
SecndObserved = as.data.frame(covBAmouse,row.names=c("theta1","theta2","theta3","lambda1","lambda2","lambda3",
                                                     "lambda4","lambda5","lambda6","lambda7","lambda8","lambda9",
                                                     "lambda10"))
colnames(SecndObserved) = c("theta1","theta2","theta3","lambda1","lambda2","lambda3","lambda4","lambda5",
                            "lambda6","lambda7","lambda8","lambda9","lambda10")
write.xlsx(SecndObserved,file="SecndObserved.xlsx",row.names=TRUE)

COVARIANCE = solve(-SecndObserved)
COVARIANCEObs = as.data.frame(COVARIANCE,
                              row.names=c("theta1","theta2","theta3","lambda1","lambda2","lambda3",
                                          "lambda4","lambda5","lambda6","lambda7","lambda8","lambda9",
                                          "lambda10"))
colnames(COVARIANCEObs) = c("theta1","theta2","theta3","lambda1","lambda2","lambda3","lambda4","lambda5",
                            "lambda6","lambda7","lambda8","lambda9","lambda10")
write.xlsx(COVARIANCEObs,file="COVARIANCEObs.xlsx",row.names=TRUE)



##################################################
###### EM MM Algorithm Variance Estimation #######
############### Using Bootstrap ##################
##################################################
dataprocessed = DATA.Reversion(BoricAcidMousedata_processed)
doselevels = DoseLevels(BoricAcidMousedata_processed)
cat("The doselevels in this dataset are",doselevels,sep=" ","\n")

# first loc is baseline A; second loc is response B; third loc is response C
Doselvlsinput = doselevels[c(4, 1,2,3)]
cat("The input doselevels to the model are",Doselvlsinput,sep=" ","\n")

modelresult34 = modelestimate(dataprocessed,Doselvlsinput,SIMULATORNOT=FALSE,RMAX=1)
cat("The model estimate results are:",sep=" ","\n")
print(modelresult34)



set.seed(123)
BOOTSTRAP = function(DATASET,ESTIMATES,B,DOSElvls) {
  Rglobal = ESTIMATES$Rglobal
  THETAS = matrix(0,nrow=B,ncol=3)
  LAMBDAS = matrix(0,nrow=B,ncol=Rglobal)
  
  for(b in 1:B) {
    SAMPLEINX = sample(1:(dim(DATASET)[1]), size=(dim(DATASET)[1]), replace=TRUE,
                       prob=rep(1,times=(dim(DATASET)[1])))
    SAMPLE = DATASET[SAMPLEINX,]
    SAMPLEresults = modelestimate(SAMPLE,DOSElvls,SIMULATORNOT=TRUE,RMAX=Rglobal)
    THETAS[b,] = c(SAMPLEresults$theta1model,SAMPLEresults$theta2model,SAMPLEresults$theta3model)
    LAMBDAS[b,] = as.numeric((SAMPLEresults$lambdaAmodel)[2:(Rglobal+1)])
  }
  Bootstrap = cbind(THETAS,LAMBDAS)
  Means = colMeans(Bootstrap,na.rm=FALSE,dims=1)
  Centered = sweep(Bootstrap,2,Means)
  COVARIANCE = (t(Centered) %*% Centered) / (B-1)
  COVARIANCE 
}

covBootstrap = BOOTSTRAP(dataprocessed,modelresult34,1000,Doselvlsinput) 
COVBootstrap = as.data.frame(covBootstrap,
                             row.names=c("theta1","theta2","theta3",
                                         "lambda1","lambda2","lambda3","lambda4","lambda5",
                                         "lambda6","lambda7","lambda8","lambda9","lambda10"))
colnames(COVBootstrap) = c("theta1","theta2","theta3",
                           "lambda1","lambda2","lambda3","lambda4","lambda5",
                           "lambda6","lambda7","lambda8","lambda9","lambda10")
write.xlsx(COVBootstrap,file="COVBootstrap.xlsx",row.names=TRUE)
write.xlsx(COVBootstrap,file="COVBootstrap1.xlsx",row.names=TRUE)


####################################
############ Simulation ############
####################################
set.seed(567)
MaxClusterSize = 10

SIMULATIONS = function(SimuTotal) {
  THETAS = matrix(0,nrow=SimuTotal,ncol=3)
  LAMBDAS = matrix(0,nrow=SimuTotal,ncol=MaxClusterSize)
  SES = matrix(0,nrow=SimuTotal,ncol=(3+MaxClusterSize))
  
  for (simu in 1:SimuTotal) {
    GENERATES = DataGenerat(qbaseline,thetaPopulat,MaxClusterSize)
    
    GENprocessed = DATA.Reversion(GENERATES)
    doselevels = DoseLevels(GENERATES)
    DOSElvls = doselevels[c(1,2,3,4)]
    
    SAMPLEest = modelestimate(GENprocessed,DOSElvls,SIMULATORNOT=TRUE,RMAX=MaxClusterSize)
    THETAS[simu,] = c(SAMPLEest$theta1model,SAMPLEest$theta2model,SAMPLEest$theta3model)
    LAMBDAS[simu,] = as.numeric((SAMPLEest$lambdaAmodel)[2:(MaxClusterSize+1)])
    
    COVGen = solve(-VARIANCE(SAMPLEest))
    SES[simu,] = sqrt(diag(COVGen))
  }
  list(THETAS=THETAS, LAMBDAS=LAMBDAS,SES=SES)
}

SIM = SIMULATIONS(100)
SIMest = SIMest[rowSums(is.na(SIMest))==0,]
SIMest = cbind.data.frame(SIM$THETAS,SIM$LAMBDAS,SIM$SES)
names(SIMest) = c("theta1","theta2","theta3",
                  "lambda1","lambda2","lambda3","lambda4","lambda5",
                  "lambda6","lambda7","lambda8","lambda9","lambda10",
                  "SEtheta1","SEtheta2","SEtheta3",
                  "SElambda1","SElambda2","SElambda3","SElambda4","SElambda5",
                  "SElambda6","SElambda7","SElambda8","SElambda9","SElambda10")
write.xlsx(SIMest,file="SIMest.xlsx",row.names=TRUE)

# summary statistics for lambda
SIMest.table = matrix(0,nrow=dim(SIMest)[2],ncol=8)

for (kk in 1:(dim(SIMest)[2])) {
  mean = mean(SIMest[,kk])
  median = median(SIMest[,kk])
  var = var(SIMest[,kk])
  sd = sd(SIMest[,kk])
  max = max(SIMest[,kk])
  min = min(SIMest[,kk])
  quantile25 = as.numeric(quantile(SIMest[,kk],probs=0.25))
  quantile75 = as.numeric(quantile(SIMest[,kk],probs=0.75))
  SIMest.table[kk,] = c(mean,median,var,sd,max,min,quantile25,quantile75)
}
SIMest.table = as.data.frame(SIMest.table)
names(SIMest.table) = c("mean","median","var","sd","max","min","quantile25","quantile75")
row.names(SIMest.table) = c("theta1","theta2","theta3",
                            "lambda1","lambda2","lambda3","lambda4","lambda5",
                            "lambda6","lambda7","lambda8","lambda9","lambda10",
                            "SEtheta1","SEtheta2","SEtheta3",
                            "SElambda1","SElambda2","SElambda3","SElambda4","SElambda5",
                            "SElambda6","SElambda7","SElambda8","SElambda9","SElambda10")
write.xlsx(SIMest.table,file="SIMestSummary.xlsx",row.names=TRUE)


####################################
###### Likelihood Ratio Test #######
####################################
# function calculating nrs*log(sum(qt*choose(R-r,t-s)/choose(R-t)))
nrslogsum = function(RMAX,NRS,QTVEC) {
  SUM1 = 0
  for(rrr in 1:RMAX) {
    for(sss in 0:rrr) {
      if(NRS[rrr,sss+1] != 0) {
        SUM2 = 0
        for(ttt in sss:(RMAX-rrr+sss)) {
          SUM2 = SUM2 + 
            QTVEC[ttt+1]*choose(RMAX-rrr,ttt-sss)/choose(RMAX,ttt)
        }
        SUM1 = SUM1 + NRS[rrr,sss+1]*log(SUM2)
      }
    }
  }
  SUM1
}

# log-likehood function
loglikelihood = function(RMAX,NRSA,QTVECA,NRSB,QTVECB,NRSC,QTVECC,NRSD,QTVECD) {
  nrslogsum(RMAX,NRSA,QTVECA) + nrslogsum(RMAX,NRSB,QTVECB) + 
    nrslogsum(RMAX,NRSC,QTVECC) + nrslogsum(RMAX,NRSD,QTVECD)
}

LRT = function(ModelResults) {
  Rglobal = ModelResults$Rglobal
  df.diff = 3*Rglobal - 3
  nrsA = ModelResults$nrsA
  nrsB = ModelResults$nrsB
  nrsC = ModelResults$nrsC
  nrsD = ModelResults$nrsD
  
  qtAnonpar = ModelResults$qtAnonpar
  qtBnonpar = ModelResults$qtBnonpar
  qtCnonpar = ModelResults$qtCnonpar
  qtDnonpar = ModelResults$qtDnonpar
  
  qtAmodel = ModelResults$qtAmodel
  qtBmodel = ModelResults$qtBmodel
  qtCmodel = ModelResults$qtCmodel
  qtDmodel = ModelResults$qtDmodel
  
  loglnonpar = loglikelihood(Rglobal,nrsA,qtAnonpar,nrsB,qtBnonpar,nrsC,qtCnonpar,nrsD,qtDnonpar)
  loglmodel  = loglikelihood(Rglobal,nrsA,qtAmodel,nrsB,qtBmodel,nrsC,qtCmodel,nrsD,qtDmodel)
  pvalue = pchisq(-2*(loglmodel-loglnonpar),df=df.diff,lower.tail=FALSE)
  
  list(loglnonpar=loglnonpar,loglmodel=loglmodel,pvalue=pvalue)
}



####################################
########### Simulation #############
####################################
set.seed(123)
# Method1: truncate lambda^hat; use inclusion-exclusion to retrieve qts
simulation1 = function(DoseLvls,ModelResults) {
  doseA = DoseLvls[1]              # baseline A
  doseB = DoseLvls[2]              # response B
  doseC = DoseLvls[3]              # response C
  doseD = DoseLvls[4]              # response D
  
  meanclustersizeA = ModelResults$meanclustersizeA
  meanclustersizeB = ModelResults$meanclustersizeB
  meanclustersizeC = ModelResults$meanclustersizeC
  meanclustersizeD = ModelResults$meanclustersizeD
  
  Rglobal = ModelResults$Rglobal
  lambdaAmodel = ModelResults$lambdaAmodel
  lambdaBmodel = ModelResults$lambdaBmodel
  lambdaCmodel = ModelResults$lambdaCmodel
  lambdaDmodel = ModelResults$lambdaDmodel
  
  # retrieve qt's using model estimate of lambda
  QfromLambda = function(clustersize,Lambdas) {
    QS = numeric(clustersize+1)
    for(s in 0:clustersize) {
      js = seq(from=0,to=(clustersize-s),by=1)
      QS[s+1] = choose(clustersize,s)*
        sum(ifelse(Lambdas[s+js+1]!=0,((-1)^js)*choose((clustersize-s),js)*Lambdas[s+js+1],0))
    }
    QS = ifelse(QS<0,0,QS)
    QS
  }
  
  # distribution of clustersize: binomial((samplemean/Rmax),Rmax)
  PERMUTATION1 = function(DOSELVL,MeanClusize,Lambdas) {
    sample = NULL
    # clustersizes are binomially distributed with n=Rglobal, p=MeanClusterSize/Rglobal
    Clusizes = seq(from=1,to=Rglobal,by=1)
    RandSizes = sample(Clusizes,size=10,replace=TRUE,
                       prob=dbinom(x=Clusizes,size=Rglobal,prob=(MeanClusize/Rglobal)))
    for(i in 1:length(RandSizes)) {
      qts = QfromLambda(RandSizes[i],Lambdas[1:(RandSizes[i]+1)])
      subsamp = matrix(0,nrow=(RandSizes[i]+1),ncol=3)
      subsamp[,1] = rep(RandSizes[i],times=RandSizes[i]+1)
      subsamp[,2] = seq(from=0,to=RandSizes[i],by=1)
      subsamp[,3] = as.numeric(rmultinom(n=1,size=15,prob=qts))
      sample = rbind(sample,subsamp)
    }
    sample = as.data.frame(sample)
    sample = cbind.data.frame(DOSELVL,sample)
    names(sample) = c("dose","nTotal","nResp","Freq")
    sample = sample[which(sample[,4]!=0),]
  }
  sampleA = PERMUTATION1(doseA,meanclustersizeA,lambdaAmodel)
  sampleB = PERMUTATION1(doseB,meanclustersizeB,lambdaBmodel)
  sampleC = PERMUTATION1(doseC,meanclustersizeC,lambdaCmodel)
  sampleD = PERMUTATION1(doseD,meanclustersizeD,lambdaDmodel)
  samplesABCD = rbind.data.frame(sampleA,sampleB,sampleC,sampleD)
  samplesABCD
}



# Method2: model estimate "qt" as probability parameter of multinomial distribution
# use hypergeometric for truncation
simulation2 = function(DoseLvls,ModelResults) {
  doseA = DoseLvls[1]              # baseline A
  doseB = DoseLvls[2]              # response B
  doseC = DoseLvls[3]              # response C
  doseD = DoseLvls[4]              # response D
  
  meanclustersizeA = ModelResults$meanclustersizeA
  meanclustersizeB = ModelResults$meanclustersizeB
  meanclustersizeC = ModelResults$meanclustersizeC
  meanclustersizeD = ModelResults$meanclustersizeD
  
  Rglobal = ModelResults$Rglobal
  qtAmodel = ModelResults$qtAmodel
  qtBmodel = ModelResults$qtBmodel
  qtCmodel = ModelResults$qtCmodel
  qtDmodel = ModelResults$qtDmodel
  
  PERMUTATION2 = function(DOSELVL,QTMODEL,MeanClusize) {
    multinoms = as.numeric(rmultinom(n=1,size=150,prob=QTMODEL))
    multresplongvec = numeric(0)
    for(resp in 0:Rglobal) {
      multresplongvec = c(multresplongvec,
                          rep(resp,times=multinoms[resp+1]))
    }
    # clustersizes are binomially distributed with n=Rglobal, p=MeanClusterSize/Rglobal
    Clusizes = seq(from=1,to=Rglobal,by=1)
    RandSizes = sample(Clusizes,size=150,replace=TRUE,
                       prob=dbinom(x=Clusizes,size=Rglobal,prob=(MeanClusize/Rglobal)))
    sample = matrix(0,nrow=length(RandSizes),ncol=3)
    sample[,1] = RandSizes
    sample[,3] = rep(1,times=length(RandSizes))
    for(recd in 1:length(RandSizes)) {
      sample[recd,2] = rhyper(nn=1,m=multresplongvec[recd],
                              n=(Rglobal-multresplongvec[recd]),k=RandSizes[recd])
    }
    sample = as.data.frame(sample)
    sample = cbind.data.frame(DOSELVL,sample)
    names(sample) = c("dose","nTotal","nResp","Freq")
    sample = sample[which(sample[,4]!=0),]
  }
  sampleA = PERMUTATION2(doseA,qtAmodel,meanclustersizeA)
  sampleB = PERMUTATION2(doseB,qtBmodel,meanclustersizeB)
  sampleC = PERMUTATION2(doseC,qtCmodel,meanclustersizeC)
  sampleD = PERMUTATION2(doseD,qtDmodel,meanclustersizeD)
  samplesABCD = rbind.data.frame(sampleA,sampleB,sampleC,sampleD)
  samplesABCD
}



# Method3: model estimate "qt" as probability parameter of multinomial distribution
# use hypergeometric for truncation
simulation3 = function(DoseLvls,ModelResults) {
  doseA = DoseLvls[1]              # baseline A
  doseB = DoseLvls[2]              # response B
  doseC = DoseLvls[3]              # response C
  doseD = DoseLvls[4]              # response D
  
  meanclustersizeA = ModelResults$meanclustersizeA
  meanclustersizeB = ModelResults$meanclustersizeB
  meanclustersizeC = ModelResults$meanclustersizeC
  meanclustersizeD = ModelResults$meanclustersizeD
  
  Rglobal = ModelResults$Rglobal
  qtAmodel = ModelResults$qtAmodel
  qtBmodel = ModelResults$qtBmodel
  qtCmodel = ModelResults$qtCmodel
  qtDmodel = ModelResults$qtDmodel
  
  PERMUTATION3 = function(DOSELVL,QTMODEL,MeanClusize) {
    ps.rt = array(0,dim=c(Rglobal+1,Rglobal,Rglobal+1))
    products = array(0,dim=c(Rglobal+1,Rglobal,Rglobal+1))
    ps.r = matrix(0,nrow=Rglobal+1,ncol=Rglobal)
    
    for(TrunClusize in 1:Rglobal) {      # truncated cluster size: r
      for(CompletResp in 0:Rglobal) {    # complete data nResponse/Rglobal : t
        lower = max(0,(TrunClusize+CompletResp-Rglobal))
        upper = min(TrunClusize,CompletResp)
        ps.rt[,TrunClusize,CompletResp+1] = dhyper(x=seq(from=0,to=Rglobal,by=1),
                                                   m=CompletResp,n=(Rglobal-CompletResp),
                                                   k=TrunClusize,log=FALSE)
      }
    }
    for(TrunResp in 0:Rglobal) {         # truncated nResponse: s
      for(TrunClusize in 1:Rglobal) {    # truncated cluster size: r
        products[TrunResp+1,TrunClusize,] = ps.rt[TrunResp+1,TrunClusize,]*QTMODEL
      }
    }
    for(TrunResp in 0:Rglobal) {         # truncated nResponse: s
      for(CompletResp in 0:Rglobal) {    # complete data nResponse/Rglobal : t  
        clusizepdf = dbinom(x=seq(from=1,to=Rglobal,by=1),size=Rglobal,prob=(MeanClusize/Rglobal))
        products[TrunResp+1,,CompletResp+1] = products[TrunResp+1,,CompletResp+1]*(clusizepdf/sum(clusizepdf))
      }
    }
    for(TrunResp in 0:Rglobal) {          # truncated nResponse: s
      for(TrunClusize in 1:Rglobal) {     # truncated cluster size: r
        ps.r[TrunResp+1,TrunClusize] = sum(products[TrunResp+1,TrunClusize,])
      }
    }
    prs = t(ps.r)
    prslongvec = as.numeric(prs)
    lablrs = outer(X=seq(from=1,to=Rglobal,by=1),Y=seq(from=0,to=Rglobal,by=1),FUN="paste")
    lablrslongvec = as.vector(lablrs)
    indlongvec = seq(from=1,to=Rglobal*(Rglobal+1),by=1)
    # sample the above long vectors with replacement with probabilities "prslongvec"
    SampInd = sample(indlongvec,replace=TRUE,size=150,prob=prslongvec)
    SampLablrs = lablrslongvec[SampInd] # nTotal & nResponse in the simulation dataset
    
    # construct simulation dataset & extract values in the long vectors
    sample = matrix(0,nrow=length(SampInd),ncol=3)
    sample[,3] = rep(1,times=length(SampInd))
    for (obs in 1:length(SampInd)) {
      obslabl =  as.numeric(unlist(strsplit(SampLablrs[obs],split=" ")))
      sample[obs,1] = obslabl[1]
      sample[obs,2] = obslabl[2]
    }
    sample = as.data.frame(sample)
    sample = cbind.data.frame(DOSELVL,sample)
    names(sample) = c("dose","nTotal","nResp","Freq")
    sample = sample[which(sample[,4]!=0),]
  }
  sampleA = PERMUTATION3(doseA,qtAmodel,meanclustersizeA)
  sampleB = PERMUTATION3(doseB,qtBmodel,meanclustersizeB)
  sampleC = PERMUTATION3(doseC,qtCmodel,meanclustersizeC)
  sampleD = PERMUTATION3(doseD,qtDmodel,meanclustersizeD)
  samplesABCD = rbind.data.frame(sampleA,sampleB,sampleC,sampleD)
  samplesABCD
}



SIMULANALYZE = function(DOSELVL,SIMULDATA,SIMULMETHOD,ModelResults,LRTresult34) {
  RMax = ModelResults$Rglobal
  simumodel34 = modelestimate(SIMULDATA,DOSELVL,SIMULATORNOT=TRUE,RMAX=RMax)
  cat("The method",SIMULMETHOD,"simulation results are:",sep=" ","\n")
  print(simumodel34)
  
  loglsimumodel34 = loglikelihood(ModelResults$Rglobal,
                                  ModelResults$nrsA,simumodel34$qtAmodel,
                                  ModelResults$nrsB,simumodel34$qtBmodel,
                                  ModelResults$nrsC,simumodel34$qtCmodel,
                                  ModelResults$nrsD,simumodel34$qtDmodel)
  cat("loglikelihood from simulation using method",SIMULMETHOD,"is",loglsimumodel34,
      "loglikelihood from model estimate is",LRTresult34$loglmodel,sep=" ","\n")
  
  # comparison of qt's
  plot(ModelResults$qtAmodel,main=paste0("Qt Comparison Model vs Simulation",SIMULMETHOD),
       ylim=c(0,1),type="b")
  points(simumodel34$qtAmodel,col="orange",type="b")
  
  points(ModelResults$qtBmodel,col="red",type="b")
  points(simumodel34$qtBmodel,col="blue",type="b")
  
  points(ModelResults$qtCmodel,col="green",type="b")
  points(simumodel34$qtCmodel,col="purple",type="b")
  
  points(ModelResults$qtDmodel,col="yellow",type="b")
  points(simumodel34$qtDmodel,col="brown",type="b")
  
  legend(5,1,legend=c("qtA from model",paste0("qtA from simulation",SIMULMETHOD),
                      "qtB from model",paste0("qtB from simulation",SIMULMETHOD),
                      "qtC from model",paste0("qtC from simulation",SIMULMETHOD),
                      "qtD from model",paste0("qtD from simulation",SIMULMETHOD)),
         col=c("black","orange","red","blue","green","purple","yellow","brown"),lty=1:8)
  
  # comparison of lambda's
  plot(ModelResults$lambdaAmodel,main=paste0("Lambda Comparison Model vs Simulation",SIMULMETHOD),
       ylim=c(0,1),type="b")
  points(simumodel34$lambdaAmodel,col="orange",type="b")
  
  points(ModelResults$lambdaBmodel,col="red",type="b")
  points(simumodel34$lambdaBmodel,col="blue",type="b")
  
  points(ModelResults$lambdaCmodel,col="green",type="b")
  points(simumodel34$lambdaCmodel,col="purple",type="b")
  
  points(ModelResults$lambdaDmodel,col="yellow",type="b")
  points(simumodel34$lambdaDmodel,col="brown",type="b")
  
  legend(5,0.8,legend=c("lambdaA from model",paste0("lambdaA from simulation",SIMULMETHOD),
                        "lambdaB from model",paste0("lambdaB from simulation",SIMULMETHOD),
                        "lambdaC from model",paste0("lambdaC from simulation",SIMULMETHOD),
                        "lambdaD from model",paste0("lambdaD from simulation",SIMULMETHOD)),
         col=c("black","orange","red","blue","green","purple","yellow","brown"),lty=1:8)
}



# DOSELEVELS: vector of numeric numbers
# first loc is baseline A; second loc is response B; third loc is response C; fourth loc is response D
# DATASET: dataset for analysis
WHOLEPROCESS = function(DOSELEVELS,DATASET) {
  dataprocessed = DATA.Reversion(DATASET)
  doselevels = DoseLevels(DATASET)
  cat("The doselevels in this dataset are",doselevels,sep=" ","\n")
  
  # first loc is baseline A; second loc is response B; third loc is response C
  Doselvlsinput = doselevels[DOSELEVELS]
  cat("The input doselevels to the model are",Doselvlsinput,sep=" ","\n")
  
  modelresult34 = modelestimate(dataprocessed,Doselvlsinput,SIMULATORNOT=FALSE,RMAX=1)
  cat("The model estimate results are:",sep=" ","\n")
  print(modelresult34)
  
  # LRT
  LRTresult34 = LRT(modelresult34)
  cat("The log-likelihood test results are:",sep=" ","\n")
  print(LRTresult34)
  
  # simulation using method1
  sim1result1 = simulation1(Doselvlsinput,modelresult34)
  datasimu1 = DATA.Reversion(sim1result1)
  SIMULANALYZE(Doselvlsinput,datasimu1,1,modelresult34,LRTresult34)
  
  # simulation using method2
  sim2result1 = simulation2(Doselvlsinput,modelresult34)
  datasimu2 = DATA.Reversion(sim2result1)
  SIMULANALYZE(Doselvlsinput,datasimu2,2,modelresult34,LRTresult34)
  
  # simulation using method3
  sim3result1 = simulation3(Doselvlsinput,modelresult34)
  datasimu3 = DATA.Reversion(sim3result1)
  SIMULANALYZE(Doselvlsinput,datasimu3,3,modelresult34,LRTresult34)
}



# Application to different dataset * combinations of lambdas A & B & C
# DOSELEVELS: vector of numeric numbers
# first loc is baseline A; second loc is response B; third loc is response C
WHOLEPROCESS(c(4,3,2,1),BoricAcidRatdata_processed)

WHOLEPROCESS(c(4,3,2,1),BBPdata_processed)

WHOLEPROCESS(c(4,3,2,1),BoricAcidMousedata_processed)

WHOLEPROCESS(c(1,4,3,2),BoricAcidRatdata_Reversed)

WHOLEPROCESS(c(2,4,3,1),BBPdata_Reversed)

WHOLEPROCESS(c(1,4,3,2),BoricAcidMousedata_Reversed)






