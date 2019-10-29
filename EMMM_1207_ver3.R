
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
#BoricAcidMousedata_processed <- read_csv("C:/Users/xinra/OneDrive/Documents/ReadingResearch/BoricAcidMousedata_processed.csv")
#DEHPmousedata_processed <- read_csv("C:/Users/Administrator/Desktop/DEHPmousedata_processed.csv")
#DEHPratdata_processed <- read_csv("C:/Users/Administrator/Desktop/DEHPratdata_processed.csv")
#DESdata_processed <- read_csv("C:/Users/Administrator/Desktop/DESdata_processed.csv")
#MEHPdata_processed <- read_csv("C:/Users/Administrator/Desktop/MEHPdata_processed.csv")


# 3D Scatterplot with Coloring and Vertical Lines and Regression Plane 
library(scatterplot3d) 

attach(BoricAcidRatdata_processed) 
s3d <-scatterplot3d(nTotal,nResp,Freq, pch=16, highlight.3d=TRUE,type="h", main="3D Scatterplot")
fit <- lm(mpg ~ wt+disp) 
s3d$plane3d(fit)

set.seed(123)
library(plotly)

#BoricAcidRatdata_processed$dose[which(BoricAcidRatdata_processed$dose==0)] <- '0% BA'
#BoricAcidRatdata_processed$dose[which(BoricAcidRatdata_processed$dose==0.1)] <- '0.1% BA'
#BoricAcidRatdata_processed$dose[which(BoricAcidRatdata_processed$dose==0.2)] <- '0.2% BA'
#BoricAcidRatdata_processed$dose[which(BoricAcidRatdata_processed$dose==0.4)] <- '0.4% BA'

p <- plot_ly(BoricAcidRatdata_processed, x = ~nTotal, y = ~nResp, z = ~Freq, 
             marker = list(color = ~dose, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
             add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Cluster size'),
                      yaxis = list(title = 'Number of responses'),
                      zaxis = list(title = 'Frequency')),
         annotations = list(
           x = 1.13,
           y = 1.05,
           text = 'Dose level (% BA)',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))
p

# stacked bar plot
BArat = DATA.Reversion(BoricAcidMousedata_processed)
BArat = cbind.data.frame(BArat,BArat$nTotal - BArat$nResp)
names(BArat) = c("Dose","nTotal","nResp","nFail")

BA0rat2 = BArat[which(BArat$Dose==0), ]
BA0rat3 = cbind.data.frame(rep("0% BA", times=2*nrow(BA0rat2)),
                           rep(1:nrow(BA0rat2), times=2), 
                           c(BA0rat2$nResp, BA0rat2$nFail),
                           c(rep("CombinedEndpoint", times=nrow(BA0rat2)), rep("Alive", times=nrow(BA0rat2)))) 
names(BA0rat3) = c("Dose","ClusterID", "Value","ResponseLevel")

BA1rat2 = BArat[which(BArat$Dose==0.1), ]
BA1rat3 = cbind.data.frame(rep("0.1% BA", times=2*nrow(BA1rat2)),
                           rep(1:nrow(BA1rat2), times=2), 
                           c(BA1rat2$nResp, BA1rat2$nFail),
                           c(rep("CombinedEndpoint", times=nrow(BA1rat2)), rep("Alive", times=nrow(BA1rat2)))) 
names(BA1rat3) = c("Dose","ClusterID", "Value","ResponseLevel")

BA2rat2 = BArat[which(BArat$Dose==0.2), ]
BA2rat3 = cbind.data.frame(rep("0.2% BA", times=2*nrow(BA2rat2)),
                           rep(1:nrow(BA2rat2), times=2), 
                           c(BA2rat2$nResp, BA2rat2$nFail),
                           c(rep("CombinedEndpoint", times=nrow(BA2rat2)), rep("Alive", times=nrow(BA2rat2)))) 
names(BA2rat3) = c("Dose","ClusterID", "Value","ResponseLevel")

BA4rat2 = BArat[which(BArat$Dose==0.4), ]
BA4rat3 = cbind.data.frame(rep("0.4% BA", times=2*nrow(BA4rat2)),
                           rep(1:nrow(BA4rat2), times=2), 
                           c(BA4rat2$nResp, BA4rat2$nFail),
                           c(rep("CombinedEndpoint", times=nrow(BA4rat2)), rep("Alive", times=nrow(BA4rat2)))) 
names(BA4rat3) = c("Dose","ClusterID", "Value","ResponseLevel")

BArat2 = rbind.data.frame(BA0rat3, BA1rat3, BA2rat3, BA4rat3)

# approach 1
d <- ggplot(BArat2, aes(fill=ResponseLevel, y=Value, x=ClusterID)) + 
  geom_bar(stat="identity", alpha=.5) +
  ggtitle("Responses From Clusters Under Different BA Concentrations")

library(lemon)
d3 <- d + facet_wrap(~Dose, ncol=2, nrow=3) + scale_color_discrete(guide=guide_legend(nrow=3))
reposition_legend(d3, 'right', panel='panel-2-2')

par(mfrow=c(2,2))

# approach 2
p1 <- ggplot(BA0rat3, aes(fill=ResponseLevel, y=Value, x=ClusterID)) + 
  geom_bar(stat="identity", alpha=.5) 

p2 <- ggplot(BA1rat3, aes(fill=ResponseLevel, y=Value, x=ClusterID)) + 
  geom_bar(stat="identity", alpha=.5) 

p3 <- ggplot(BA2rat3, aes(fill=ResponseLevel, y=Value, x=ClusterID)) + 
  geom_bar(stat="identity", alpha=.5) 

p4 <- ggplot(BA4rat3, aes(fill=ResponseLevel, y=Value, x=ClusterID)) + 
  geom_bar(stat="identity", alpha=.5) 

grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 2, nrow = 2, position='right')


# transform cbdata to GEE-like form
DATA.Expand1 = function(cbdata) {
  attach(cbdata)
  datareversion = function(columnname,reptimes) {
    repeatterm = rep(columnname,times=reptimes)
  }
  dataprocessed = cbind.data.frame(datareversion(simul,Freq),
                                   datareversion(dose,Freq),
                                   datareversion(nTotal,Freq),
                                   datareversion(nResp,Freq))
  names(dataprocessed) = c("simul","dose","nTotal","nResp")
  detach(cbdata)
  dataprocessed
}

DATA.Expand2 = function(cbdata) {
  attach(cbdata)
  datareversion = function(columnname,reptimes) {
    repeatterm = rep(columnname,times=reptimes)
  }
  dataprocessed = cbind.data.frame(datareversion(simul,freq),
                                   datareversion(case,freq),
                                   datareversion(dose,freq),
                                   datareversion(resp,freq))
  detach(cbdata)
  dataprocessed
}

GEETrans = function(cbdata) {
  cbdataTrans = DATA.Expand1(cbdata)
  cbdataTrans = cbind.data.frame(rep(1:(4*250),times=1000),cbdataTrans)
  names(cbdataTrans) = c("case","simul","dose","nTotal","nResp")
  # cbdata for GEE
  cbdataGEE = matrix(0,nrow=2*dim(cbdataTrans)[1],ncol=5)
  cbdataGEE[,4] = rep(c(1,0),times=dim(cbdataTrans)[1])
  for (obs in 1:(2*dim(cbdataTrans)[1])) {
    if (obs %% 2 == 1) {
      cbdataGEE[obs,1] = cbdataTrans[(obs+1)/2,2]
      cbdataGEE[obs,2] = cbdataTrans[(obs+1)/2,1]
      cbdataGEE[obs,3] = cbdataTrans[(obs+1)/2,3]
      cbdataGEE[obs,5] = cbdataTrans[(obs+1)/2,5]
    } else {
      cbdataGEE[obs,1] = cbdataTrans[obs/2,2]
      cbdataGEE[obs,2] = cbdataTrans[obs/2,1]
      cbdataGEE[obs,3] = cbdataTrans[obs/2,3]
      cbdataGEE[obs,5] = cbdataTrans[obs/2,4] - cbdataTrans[obs/2,5]
    }
  }
  cbdataGEE = as.data.frame(cbdataGEE)
  names(cbdataGEE) =  c("simul","case","dose","resp","freq")
  cbdataGEE = DATA.Expand2(cbdataGEE)
  names(cbdataGEE) =  c("simul","case","dose","resp")
  write.csv(cbdataGEE,file="C:/Users/xinqi/Desktop/cbdataGEE.csv",row.names=FALSE)
}
GEETrans(GENERATE)



# generate cbdata
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

GENERATE = NULL
for(simulat in 1:1000) {
  SimulData = cbind.data.frame(simulat,DataGenerat(qbaseline,thetaPopulat,MaxClusterSize)) 
  GENERATE = rbind.data.frame(GENERATE, SimulData)
}
names(GENERATE) = c("simul","dose","nTotal","nResp","Freq")



### Reverse coding of responses ###
ReverseCoding = function(cbdata) {
  cbdata$nResp =  cbdata$nTotal - cbdata$nResp
  cbdata
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
DATA.Reversion = function(cbdata) {
  attach(cbdata)
  datareversion = function(columnname,reptimes) {
    repeatterm = rep(columnname,times=reptimes)
  }
  dataprocessed = cbind.data.frame(datareversion(dose,Freq),
                                   datareversion(nTotal,Freq),
                                   datareversion(nResp,Freq))
  names(dataprocessed) = c("dose","nTotal","nResp")
  detach(cbdata)
  dataprocessed
}



# doselvls are values of solution concentrations
DATA.Subset = function(cbdata,doselvls) {
  attach(cbdata)
  dosesubset = NULL
  # Clusters partitioning with respect to dose levels
  for (lvls in 1:length(doselvls)) {
    dosesubset = rbind.data.frame(dosesubset,
                                  subset(cbdata,dose==doselvls[lvls]))
  }
  names(dosesubset) = c("dose","nTotal","nResp")
  detach(cbdata)
  dosesubset
}



# returns a vector of solution concentrations
DoseLevels = function(cbdata) {
  sort(unique(cbdata$dose),decreasing=FALSE)
}



# Function to estimate qk
# where qk represents: for the complete data setting, given all cluster sizes are Rglobal,
# the probability of achieving k successes/responses
# doselevel are values of solution concentrations
qestimate = function(cbdata,doselevel,RMAX,itermax,eps) {
  BARat = DATA.Subset(cbdata,doselevel)
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
  #NrsS = xtabs(~nTotal+nResp, data=BARat)
  
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
  while((difference > eps) && (iter < itermax)) {
    iter = iter + 1
    qnext = maximization(q0)
    difference = sum(abs(qnext-q0))
    q0 = qnext
  }
  # cat("The estimated q0 is",q0, "\n")
  
  # Calculation of p(r)
  sequence = seq(0:Rglobal)
  probability = ifelse(q0!=0,q0/choose(Rglobal,sequence),0) 
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
  
  list(q=q0,
       p=probability,
       lambda=lambdaest,
       ntprime=ntprime,
       nrs=nrs.matrix,
       meanclustersize = meanclustersize)
}
 

  
QestGenerate1 = qestimate(GENERATE,1,MaxClusterSize,5000,1e-7)
QestGenerate2 = qestimate(GENERATE,2,MaxClusterSize,5000,1e-7)
QestGenerate3 = qestimate(GENERATE,3,MaxClusterSize,5000,1e-7)
QestGenerate0 = qestimate(GENERATE,0,MaxClusterSize,5000,1e-7)

plot(QestGenerate1$lambda,col="black",ylim=c(0,1))
lines(QestGenerate2$lambda,col="red")
lines(QestGenerate3$lambda,col="green")
lines(QestGenerate0$lambda,col="blue")

plot(QestGenerate1$q,col="black",ylim=c(0,0.35))
lines(QestGenerate2$q,col="red")
lines(QestGenerate3$q,col="green")
lines(QestGenerate0$q,col="blue")



# inputs: cbdata & baseline dose level & Maximum cluster size
#         maximum iterations & convergence criterion 
#         if we want to estimate the covariance or not
SP.est = function(cbdata,baseline=NULL,MaxClusSize=NULL,itermax=5000,eps=1e-7,est.var=FALSE) {
  
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
  # prst.vec = function(rr,ss,qvec) {
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
  
  # Maximum cluster size within the extracted doselevels
  # by default, if MaxClusSize is not specified (NA), Rglobal = max(cbdata$nTotal)
  # otherwise, Rglobal = user specified 
  Rglobal = ifelse(is.null(MaxClusSize),max(cbdata$nTotal),MaxClusSize)
  
  # doselvlsinput: real concentrations of solutions
  # 1st is baseline A, followed by dose levels in an decreasing order of concentrations
  # by default, if baseline is not specified (NULL), then the highest concentration is selected as the baseline
  concentrations = unique(sort(cbdata$dose,decreasing=TRUE))
  Dlvls = length(concentrations)
  
  # Intial Parameter storage
  InitMeanCluSizeS = numeric(Dlvls)
  InitLambdaNonParS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  InitQNonParS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  InitNrsS = matrix(0,nrow=(Dlvls*Rglobal),ncol=(Rglobal+1))
  InitNtprimeNonParS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  
  for (dlvl in 1:Dlvls) {
    Result = qestimate(cbdata,concentrations[dlvl],Rglobal,itermax,eps)  # response
    Meanclustersize = Result$meanclustersize
    Lambda = Result$lambda   # nonparametric estimate for comparison
    Qt.vec = Result$q                
    Nrs = Result$nrs 
    Ntprime = Result$ntprime

    # Parameter storasge for each dose level
    InitMeanCluSizeS[dlvl] = Meanclustersize
    InitLambdaNonParS[dlvl,] = Lambda
    InitQNonParS[dlvl,] = Qt.vec
    InitNrsS[((1:Rglobal)+(dlvl-1)*Rglobal),] = Nrs
    InitNtprimeNonParS[dlvl,] = Ntprime
  }

  Baseline = ifelse(is.null(baseline),concentrations[which.max(InitLambdaNonParS[,2])],baseline)
  doselvlsinput = c(Baseline, concentrations[-which(concentrations==Baseline)])
  Pointer = c(which(concentrations==Baseline),(1:Dlvls)[-which(concentrations==Baseline)])
  
  # assuming uniform effectiveness within two specifc doselevels
  # initial value of lambda for the baseline
  # baseline, higher values in this situation 
  #          resultpooled = qestimate(cbdata,doselvlsinput,Rglobal,itermax,eps)
  #          lambdaApooled = resultpooled$lambda
  # reason that qtA also pooled over all doselevels: lambdaA never used in EM algorithm
  #          qtApooled = resultpooled$q          # baseline
  # notice the order of input doselevels
 
  meanclustersizeA = InitMeanCluSizeS[Pointer[1]]
  lambdaAnotpooled = InitLambdaNonParS[Pointer[1],]  # baseline A, nonparametric estimate for comparison
  # Initial values for qt's from group A & B
  qtA.notpooled = InitQNonParS[Pointer[1],]          # baseline A, nonparametric estimate for comparison
  qtA.vec = qtA.notpooled            # baseline qtA.notpooled  # qtApooled before
  nrsA = InitNrsS[((1:Rglobal)+(Pointer[1]-1)*Rglobal),]       # baseline A
  # Calculate ntprime's for group A based on non-parametric model
  ntprimeAk = InitNtprimeNonParS[Pointer[1],]
  # which is essentially the original values of ntprimeA = resultA$ntprime
  
  # Parameter storage
  MeanCluSizeS = numeric(Dlvls)
  LambdaModelS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  LambdaNonParS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  
  QAModelInitS = numeric(Rglobal+1)
  QAModelTraceS = numeric(Rglobal+1)
  QNonParS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  QModelS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  
  NrsS = matrix(0,nrow=(Dlvls*Rglobal),ncol=(Rglobal+1))
  
  ThetaInitS = numeric(Dlvls-1)
  ThetaTraceS = numeric(Dlvls-1)
  
  NtprimeInitS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  NtprimeTraceS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  NtprimeModelS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  NtprimeNonParS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
  
  MeanCluSizeS[1] = meanclustersizeA
  LambdaNonParS[1,] = lambdaAnotpooled
  QNonParS[1,] = qtA.notpooled
  QAModelInitS = qtA.notpooled
  NrsS[1:Rglobal,] = nrsA
  NtprimeNonParS[1,] = ntprimeAk
  NtprimeInitS[1,] = ntprimeAk
  
  for (dlvl in 2:Dlvls) {
    resultB = qestimate(cbdata,doselvlsinput[dlvl],Rglobal,itermax,eps)  # response
    meanclustersizeB = InitMeanCluSizeS[Pointer[dlvl]]
    lambdaB = InitLambdaNonParS[Pointer[dlvl],]   # nonparametric estimate for comparison
    qtB.vec = InitQNonParS[Pointer[dlvl],]                
    nrsB = InitNrsS[((1:Rglobal)+(Pointer[dlvl]-1)*Rglobal),]
    
    # Initial value for theta's
    theta1 = lambdaB[2] / lambdaAnotpooled[2]
    # response initial value according to seme-parametric model
    qtB.transformed = qBqATransform(theta1,QAModelInitS) 
    ntprimeBk = ntprime.est(qtB.transformed,nrsB)
    # different from ntprimeB = resultB$ntprime since qtB.transformed is used
    
    # Parameter storasge fir each dose level
    MeanCluSizeS[dlvl] = meanclustersizeB
    LambdaNonParS[dlvl,] = lambdaB
    QNonParS[dlvl,] = qtB.vec
    NrsS[((1:Rglobal)+(dlvl-1)*Rglobal),] = nrsB
    ThetaInitS[dlvl-1] = theta1
    NtprimeNonParS[dlvl,] = ntprime.est(qtB.vec,nrsB)
    NtprimeInitS[dlvl,] = ntprimeBk
  }
  
  ### Q function ###
  # component1 represents qvec
  # component2 represents ntprime
  ntlog = function(component1,component2) {
    ifelse(component2>0,log(component1)*component2,0)
  }
  # Qfunction = sum(ntlog(qtA.vec,ntprimeAk) + ntlog(qtB.transformed,ntprimeBk) +
  #                   ntlog(qtC.transformed,ntprimeCk) + ntlog(qtD.transformed,ntprimeDk))
  
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
  
  # update of qtA's using EM MM algorithm
  qtAupdate = function(THETAS,QVEC,NTPRIMEModels) {
    DlvlMinus1 = length(THETAS)
    ProductSum = matrix(0,nrow=DlvlMinus1,ncol=Rglobal+1)
    NUMERATOR = numeric(Rglobal+1)
    
    for(BETA in 0:Rglobal) {
      weightBETAS = matrix(0,nrow=DlvlMinus1,ncol=(BETA+1))
      for(RHO in 0:BETA) {
        for(dlvlMinus in 1:DlvlMinus1) {
          weightBETAS[dlvlMinus,RHO+1] = weight(RHO,THETAS[dlvlMinus],QVEC)[BETA+1]
        }
      }
      for(dlvls in 1:DlvlMinus1) {
        ProductSum[dlvls,BETA+1] = sum(weightBETAS[dlvls,] * (NTPRIMEModels[(dlvls+1),1:(BETA+1)]))
      }
      NUMERATOR[BETA+1] = NTPRIMEModels[1,BETA+1] + sum(ProductSum[,BETA+1])
    }
    NUMERATOR / sum(NUMERATOR)
  }

  ## Iterations of EM algorithm until convergence ##
  # Absolute convergence: difference1 & difference1 <= 1e-7
  iter = 0
  difference = 10
  
  while((difference > eps) & (iter < itermax)) {
    iter = iter + 1
    # Initial value for qtA and theta's
    ThetaTraceS = ThetaInitS
    QAModelTraceS = QAModelInitS
    
    ###########################################################
    ##################### EM MM ALGORITHM #####################
    ######## Expectation Maximization Minorize-Maximize #######
    ###########################################################
    
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
    QsModelEst = matrix(0,nrow=(Dlvls-1),ncol=(Rglobal+1))
    for(dlvlMinus in 1:(Dlvls-1)) {
      QsModelEst[dlvlMinus,] = qBqATransform(ThetaTraceS[dlvlMinus],QAModelTraceS) 
    }
    
    # expectation 
    NtprimeTraceS = matrix(0,nrow=Dlvls,ncol=(Rglobal+1))
    NtprimeTraceS[1,] = ntprime.est(QAModelTraceS,nrsA)  
    for(dlvlMinus in 1:(Dlvls-1)) {
      NtprimeTraceS[dlvlMinus+1,] = ntprime.est(QsModelEst[dlvlMinus,],NrsS[((1:Rglobal)+dlvlMinus*Rglobal),])  
    }
    
    ##############################
    ##### MAXIMIZATION STEP ######
    # Minorize-Maximize inserted #
    ##############################
    # theta,qtA.vec: from the pre kth iteration
    # ntprimeBk: expectaion from the pre kth iteration
    # maximization (MM inserted)
    # theta1&2pre,qtAspre: from previous kth iteration
    # ntprimeAkpre,ntprimeBkpre,ntprimeCkpre: expectaion from previous kth iteration
    for(dlvlMinus in 1:(Dlvls-1)) {
      ThetaInitS[dlvlMinus] = thetaupdate1(ThetaTraceS[dlvlMinus],QAModelTraceS,NtprimeTraceS[dlvlMinus+1,])
    }
    QAModelInitS = qtAupdate(ThetaTraceS,QAModelTraceS,NtprimeTraceS)
    difference = sum(abs(QAModelInitS-QAModelTraceS)) + sum(abs(ThetaInitS-ThetaTraceS)) 
  }
 
  # using optimized value of qtA to obtain corresponding lambda's for baseline (group A)
  lambdabaseline = numeric(Rglobal+1)
  for(tt in 0:Rglobal) {
    sumttjj = 0
    for(jj in 0:(Rglobal-tt)) {
      sumttjj = sumttjj + choose((Rglobal-tt),jj)*QAModelInitS[Rglobal-jj+1]/choose(Rglobal,(Rglobal-jj))
    }
    lambdabaseline[tt+1] = sumttjj
  }
  #cat("lambdaA from model is",lambdabaseline,sep=" ","\n")
  #plot(lambdabaseline,main="Lambda Comparison Semi-parametric vs Non-parametric",ylim=c(0,1),type="b")
  #points(LambdaNonParS[1,],col="orange",type="b")
  
  LambdaModelS[1,] = lambdabaseline
  # using qtA and theta from EM optimization results to obtain lambda's for response
  exponents = seq(from=0,to=Rglobal,by=1)
  
  for(tt in 0:Rglobal) {
    for(LVLs in 2:Dlvls) {
      LambdaModelS[LVLs,tt+1] = ((ThetaInitS[LVLs-1])^(exponents[tt+1])) * lambdabaseline[tt+1]
    }
  }
  
  QModelS[1,] = QAModelInitS
  for(LVLs in 2:Dlvls) {
    QModelS[LVLs,] = qBqATransform(ThetaInitS[LVLs-1],QAModelInitS)
  }
  
  if (!est.var) {
    list(Rglobal=Rglobal,
         doselvlsinput = doselvlsinput,
         MeanCluSizeS = MeanCluSizeS,
         NrsS = NrsS,
         QNonParS = QNonParS,
         QModelS = QModelS,
         ThetaEst = ThetaInitS,
         LambdaNonParS = LambdaNonParS,
         LambdaModelS = LambdaModelS,
         QBaseline = QAModelInitS)
  } else {
    ##################################################
    ###### EM MM Algorithm Variance Estimation #######
    ######## Using the Observed Information ##########
    ##################################################
    # THETA123: interger indicating dose level category
    # LAMBDAalpha: interger indicating order in baseline lambda estimates
    lambdaAmodel = LambdaModelS[1,]
    p1 = function(ALPHA,RR,SS,THETA123) {
      P11 = numeric(RR-SS+1)
      P12 = numeric(RR-SS+1)
      P13 = numeric(RR-SS+1)
      P14 = numeric(RR-SS+1)
      P15 = numeric(RR-SS+1)
      for(j in 0:(RR-SS)) {
        P11[j+1] = ((-1)^j) * choose(RR-SS,j) * (lambdaAmodel[SS+j+1]) * ((ThetaInitS[THETA123])^(SS+j))
        P12[j+1] = ((-1)^j) *(SS+j)*(SS+j-1)*choose(RR-SS,j)*(lambdaAmodel[SS+j+1])*((ThetaInitS[THETA123])^(SS+j-2))
        P13[j+1] = ((-1)^j) *(SS+j)*choose(RR-SS,j)*(lambdaAmodel[SS+j+1])*((ThetaInitS[THETA123])^(SS+j-1))
        P14[j+1] = ((-1)^j) *(ALPHA-SS-j)* choose(RR-SS,j) * (lambdaAmodel[SS+j+1]) * ((ThetaInitS[THETA123])^j)
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
      NRS = NrsS[(1:Rglobal)+(THETA123*Rglobal), ]
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
          
          THETAgSquare[rr,ss+1] = ifelse(part11==0, 0, part2[rr,ss+1] * (part12*part11-(part13^2)) / (part11^2) )
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
          
          ThetaLambda[rr,ss+1] = ifelse(part11==0, 0, 
                                        part2[rr,ss+1] * part14 * ((-1)^(LAMBDAalpha-ss)) * choose(rr-ss,LAMBDAalpha-ss) * ((ThetaInitS[THETA123])^(LAMBDAalpha+ss-1)) / (part11^2) )
        }
      }
      ThetaLambdasum = sum(ThetaLambda)
      ThetaLambdasum
    }
    
    p3 = function(THETA123,ALPHA,RR,SS,NRSG) {
      P31 = NRSG[RR,SS+1] * ((-1)^(2*(ALPHA-SS)+1)) * (choose(RR-SS,ALPHA-SS)^2)
      P32 = NRSG[RR,SS+1] * ((-1)^(2*(ALPHA-SS)+1)) * (choose(RR-SS,ALPHA-SS)^2) * ((ThetaInitS[THETA123])^(2*ALPHA))
      list(P31=P31,P32=P32)
    }
    
    Dlambda1lambda1 = function(LAMBDAalpha) {
      LambdaSquare1 = matrix(0,nrow=Rglobal,ncol=Rglobal+1)
      LambdaSquare2 = array(0,dim=c(Rglobal,Rglobal+1,3))
      
      for(rr in LAMBDAalpha:Rglobal) {
        for(ss in 0:min(rr,LAMBDAalpha)) {
          for(gg in 1:length(ThetaInitS)) {
            part1 = p1(LAMBDAalpha,rr,ss,gg)
            part11 = part1$P11sum
            part15 = part1$P15sum
            part2 = p2(gg)
            part3 = p3(gg, LAMBDAalpha, rr, ss, part2)
            part32 = part3$P32
            
            LambdaSquare2[rr,ss+1,gg] = ifelse(part11==0, 0, part32 / (part11^2) )
          }
          part30 = p3(1, LAMBDAalpha, rr, ss, nrsA)$P31
          LambdaSquare1[rr,ss+1] = ifelse(part15==0, 0, part30 / (part15^2) )
        }
      }
      LambdaSquaresum = sum(LambdaSquare1) + sum(LambdaSquare2)
      LambdaSquaresum
    }
    
    p4 = function(THETA123,ALPHA,BETA,RR,SS,NRSG) {
      P41 = NRSG[RR,SS+1] * ((-1)^(ALPHA+BETA-2*SS+1)) *choose(RR-SS,ALPHA-SS)*choose(RR-SS,BETA-SS)
      P42 = NRSG[RR,SS+1] * ((-1)^(ALPHA+BETA-2*SS+1)) *choose(RR-SS,ALPHA-SS)*choose(RR-SS,BETA-SS) * ((ThetaInitS[THETA123])^(ALPHA+BETA))
      list(P41=P41,P42=P42)
    }
    
    Dlambda1lambda2 = function(LAMBDAalpha,LAMBDAbeta) {
      LambdaSquare3 = matrix(0,nrow=Rglobal,ncol=Rglobal+1)
      LambdaSquare4 = array(0,dim=c(Rglobal,Rglobal+1,length(ThetaInitS)))
      
      for(rr in max(LAMBDAalpha,LAMBDAbeta):Rglobal) {
        for(ss in 0:min(rr,LAMBDAalpha,LAMBDAbeta)) {
          for(gg in 1:length(ThetaInitS)) {
            part1 = p1(LAMBDAalpha,rr,ss,gg)
            part11 = part1$P11sum
            part15 = part1$P15sum
            part2 = p2(gg)
            part4 = p4(gg, LAMBDAalpha, LAMBDAbeta, rr, ss, part2)
            part42 = part4$P42
            
            LambdaSquare4[rr,ss+1,gg] = ifelse(part11==0, 0, part42 / (part11^2) )
          }
          part40 = p4(1, LAMBDAalpha, LAMBDAbeta, rr, ss, nrsA)$P41
          LambdaSquare3[rr,ss+1] = ifelse(part15==0, 0, part40 / (part15^2) )
        }
      }
      Lambda1lambda2sum = sum(LambdaSquare3) + sum(LambdaSquare4)
      Lambda1lambda2sum
    }
    
    SecDerivat = matrix(0,nrow=(length(ThetaInitS)+Rglobal),ncol=(length(ThetaInitS)+Rglobal))
    for (i in 1:length(ThetaInitS)) {
      SecDerivat[i,i] = Dthetatheta(i)
    }
    
    for (i in 1:length(ThetaInitS)) {
      for (j in 1:Rglobal) {
        entry = Dthetalambda(i,j)
        SecDerivat[i,j+length(ThetaInitS)] = entry
        SecDerivat[j+length(ThetaInitS),i] = entry
      }
    }
    
    for (j in 1:Rglobal) {
      SecDerivat[j+length(ThetaInitS),j+length(ThetaInitS)] = Dlambda1lambda1(j)
    }
    
    for (j in 1:Rglobal) {
      for (k in 1:Rglobal) {
        SecDerivat[j+length(ThetaInitS),k+length(ThetaInitS)] = Dlambda1lambda2(j,k)
      } 
    }
    CovMatr = solve(-SecDerivat)
    
    RowNam = c(paste0("theta", 1:length(ThetaInitS)), paste0("lambda", 1:(length(lambdaAmodel)-1) ))
    CovMatr = as.data.frame(CovMatr,row.names=RowNam)
    colnames(CovMatr) = RowNam
    
    list(Rglobal=Rglobal,
         doselvlsinput = doselvlsinput,
         MeanCluSizeS = MeanCluSizeS,
         NrsS = NrsS,
         QNonParS = QNonParS,
         QModelS = QModelS,
         ThetaEst = ThetaInitS,
         LambdaNonParS = LambdaNonParS,
         LambdaModelS = LambdaModelS,
         QBaseline = QAModelInitS,
         CovMatr = CovMatr)
    }
}






### Application ###
# for example, the simulated cbdata
Thetas = matrix(0,nrow=1000,ncol=3)
SEs = matrix(0,nrow=1000,ncol=3)
for (i in 1:1000) {
  CBData = DATA.Reversion(GENERATE[((1:(4*250))+(i-1)*(4*250)), ])
  fit = SP.est(CBData,baseline=NULL,MaxClusSize=NULL,itermax=5000,eps=1e-7,est.var=TRUE) 
  Thetas[i,] = fit$ThetaEst
  SEs[i,] = sqrt(as.numeric(diag(as.matrix(fit$CovMatr))[1:3]))
}
colnames(Thetas) = c("theta1","theta2","theta3")
colnames(SEs) = c("se1","se2","se3")
write.csv(Thetas,file="C:/Users/xinqi/Desktop/Thetas0719.csv",row.names=FALSE)
write.csv(SEs,file="C:/Users/xinqi/Desktop/SEs0719.csv",row.names=FALSE)

# comparison with GEE estimates
library(haven)
GEEest = read_sas("C:/Users/xinqi/Desktop/SimulationGEE/geeest.sas7bdat",NULL)[,1:5]

GEEBetas = matrix(0,nrow=100,ncol=3)
GEEThetas = matrix(0,nrow=100,ncol=3)
GEESEs = matrix(0,nrow=100,ncol=3)
GEEThetaSEs = matrix(0,nrow=100,ncol=3)

attach(GEEest)
for (i in 1:dim(GEEest)[1]) {
  index = (i %/% 5) + 1
  if (i %% 5 == 2) {
    GEEBetas[index,1] = Estimate[i]
    GEESEs[index,1] = Stderr[i]
  } 
  if (i %% 5 == 3) {
    GEEBetas[index,2] = Estimate[i]
    GEESEs[index,2] = Stderr[i]
  } 
  if (i %% 5 == 4) {
    GEEBetas[index,3] = Estimate[i]
    GEESEs[index,3] = Stderr[i]
  }
}
detach(GEEest)
GEEThetas = exp(GEEBetas)[ ,3:1]
GEEThetaSEs = GEEThetas * (GEESEs[,3:1])
colnames(GEEThetas) = c("theta1","theta2","theta3")
colnames(GEEThetaSEs) = c("se1","se2","se3")
write.csv(GEEThetas,file="C:/Users/xinqi/Desktop/GEEThetas.csv",row.names=FALSE)
write.csv(GEEThetaSEs,file="C:/Users/xinqi/Desktop/GEEThetaSEs.csv",row.names=FALSE)



# real data analysis
CBData = DATA.Reversion(BoricAcidMousedata_processed)
doselevels = DoseLevels(BoricAcidMousedata_processed)

#CBData = DATA.Reversion(BoricAcidRatdata_processed)
#doselevels = DoseLevels(BoricAcidRatdata_processed)
cat("The doselevels in this cbdata are",doselevels,sep=" ","\n")

fit = SP.est(CBData,baseline=NULL,MaxClusSize=NULL,itermax=5000,eps=1e-7,est.var=FALSE) 
cat("The model estimate results are:",sep=" ","\n")
print(fit)



#######################################################
### use bootstrap for parameter variance estimation ###
#######################################################

set.seed(123)

BOOTSTRAP = function(DATASET,ESTIMATES,B,DOSElvls) {
  Rglobal = ESTIMATES$Rglobal
  THETAS = matrix(0,nrow=B,ncol=3)
  LAMBDAS = matrix(0,nrow=B,ncol=Rglobal)
  
  for(b in 1:B) {
    SAMPLEINX = sample(1:(dim(DATASET)[1]), size=(dim(DATASET)[1]), replace=TRUE,
                       prob=rep(1,times=(dim(DATASET)[1])))
    SAMPLE = DATASET[SAMPLEINX,]
    SAMPLEresults = SP.est(SAMPLE,baseline=0.4,MaxClusSize=Rglobal,itermax=5000,eps=1e-4,est.var=FALSE) 
    
    THETAS[b,] = SAMPLEresults$ThetaEst
    LAMBDAS[b,] = as.numeric((SAMPLEresults$LambdaModelS)[1,2:(Rglobal+1)])
  }
  Bootstrap = cbind(THETAS,LAMBDAS)
  
  list(Bootstrap = Bootstrap,
       THETAS = THETAS,
       LAMBDAS = LAMBDAS)
}

covBootstrap = BOOTSTRAP(CBData,fit,1000,Doselvlsinput) 

Means = colMeans(covBootstrap$Bootstrap,na.rm=FALSE,dims=1)
Centered = sweep(covBootstrap$Bootstrap,2,Means)
COVARIANCE = (t(Centered) %*% Centered) / (B-1)
COVARIANCE 

COVBootstrap = as.data.frame(COVARIANCE,
                             row.names=c("theta1","theta2","theta3",
                                         "lambda1","lambda2","lambda3","lambda4","lambda5",
                                         "lambda6","lambda7","lambda8","lambda9","lambda10"))
colnames(COVBootstrap) = c("theta1","theta2","theta3",
                           "lambda1","lambda2","lambda3","lambda4","lambda5",
                           "lambda6","lambda7","lambda8","lambda9","lambda10")
write.xlsx(COVBootstrap,file="COVBootstrap.xlsx",row.names=TRUE)
write.xlsx(COVBootstrap,file="COVBootstrap1.xlsx",row.names=TRUE)



