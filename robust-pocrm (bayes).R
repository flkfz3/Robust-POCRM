#############使用贝叶斯方法的Robust-POCRM/BMA-POCRM
#rgetwm函数：获取包含所研究骨架、顺序信息的列表alpha，用于rpocrmb.imp、rpocrmb.sim，见使用最大似然方法的Robust-POCRM/BMA-POCRM
#rpocrmb.imp函数：用使用贝叶斯方法的BMA-POCRM或Robust-POCRM方法推荐下一批患者的治疗剂量
#rpocrmb.imp(alpha,模型先验矩阵,目标毒性率,患者DLT数据,患者剂量数据,方法BMA/BMS,参数a预设方差)
rpocrmb.imp <-
  function(alpha,prior.o,theta,y,combos,method,scale){
    data<-as.matrix(table(combos,y))
    level<-as.numeric(row.names(data))
    nontox<-as.numeric(data[,1])
    tox<-as.numeric(data[,2])
    denpred<-apred<-lik<-matrix(0,nrow=nrow(prior.o),ncol=ncol(prior.o))
    crmh<-function(a,salpha,level,tox,nontox,scale){
      la<-exp(-a^2/2/scale^2) 
      for(i in level){
        index<-match(i,level)
        la<-la*((salpha[i]**exp(a))**tox[index])*((1-salpha[i]**exp(a))**nontox[index])
      }
      la
    }
    crmht<-function(a,salpha,level,tox,nontox,scale){
      la<-a * exp(-a^2/2/scale^2) 
      for(i in level){
        index<-match(i,level)
        la<-la*((salpha[i]**exp(a))**tox[index])*((1-salpha[i]**exp(a))**nontox[index])
      }
      la
    }
    for(j in 1:length(alpha)){
      for(k in 1:nrow(alpha[[j]])){			
        salpha<-alpha[[j]][k,]
        denpred[j,k]<-integrate(crmh,-Inf,Inf,salpha,level,tox,nontox,scale,abs.tol=0)[[1]]
        apred[j,k]<-integrate(crmht,-10,10,salpha,level,tox,nontox,scale,abs.tol=0)[[1]] / denpred[j,k]
      }
    }
    pord<-(denpred*prior.o)/sum(denpred*prior.o)
    if(method=="BMS"){
      ske<-which(pord==max(pord),arr.ind=T)[1]
      ord<-which(pord==max(pord),arr.ind=T)[2]
      ahat<-apred[ske,ord]
      rpred<-alpha[[ske]][ord,]**exp(ahat)
      next.lev<-which.is.max(-(abs(rpred-theta)))
      out<-list(ord.prob=round(pord,3),skeleton.est=ske,order.est=ord,a.est=round(ahat,3),ptox.est=round(rpred,3),dose.rec=next.lev)
    }
    else if(method=="BMA"){
      spord<-apply(pord,1,sum)
      opord<-apply(pord,2,sum)
      ord<-which.is.max(opord)
      ahat<-apred[,ord]
      srpred<-matrix(0,ncol=ncol(alpha[[1]]),nrow=length(alpha))
      for (i in 1:length(alpha)){
        srpred[i,]<-alpha[[i]][ord,]**exp(ahat[i])
      }
      rpred<-apply(srpred*spord,2,sum)
      next.lev<-which.is.max(-(abs(rpred-theta)))
      out<-list(ord.prob=round(pord,3),order.est=ord,a.est=round(ahat,3),ptox.est=round(rpred,3),dose.rec=next.lev)
    }
    else{
      print("Incorrect parameter of method!")
    }
  }
#rpocrmb.sim函数：使用贝叶斯方法的BMA-POCRM或Robust-POCRM方法的模拟试验
#rpocrmb.sim(真实毒性率,alpha,模型先验矩阵,第一阶段升剂量方案,停止样本量,最大样本量,目标毒性率,模拟次数,PCS区间半宽,方法BMA/BMS,参数a预设方差)
rpocrmb.sim<-function(r,alpha,prior.o,x0,stop,n,theta,nsim,tox.range,method,scale){
  sim <- sim1 <- apred <- lik <- pord <-ske <- ord <- ahat <- rpred <- next.lev <- n1 <- N <- NULL
  d<-ncol(alpha[[1]])
  s<-nrow(alpha[[1]])
  if(nsim==1){
    ###LOAD FUNCTION 'crm'
    #crm函数：输出一次BMS/BMA-POCRM的推荐结果
    crm<-function(obs,alpha,prior.o,theta,method,scale){
      sim<<-table(obs$level,obs$tox)
      ifelse(dim(sim)[2]==1,sim1<<-data.frame(as.numeric(row.names(sim)),sim[,1],1-sim[,1]),sim1<<-data.frame(as.numeric(row.names(sim)),sim[,1],sim[,2]))
      names(sim1)<<-c('level','nontox','tox')
      denpred<<-matrix(0,nrow=nrow(prior.o),ncol=ncol(prior.o))
      apred<<-matrix(0,nrow=nrow(prior.o),ncol=ncol(prior.o))
      lik<<-matrix(0,nrow=nrow(prior.o),ncol=ncol(prior.o))
      crmh<-function(a,salpha,level,tox,nontox,scale){
        la<-exp(-a^2/2/scale^2) 
        for(i in sim1$level){
          index<-match(i,sim1$level)
          la<-la*(salpha[i]**exp(a)**sim1$tox[index])*((1-salpha[i]**exp(a))**sim1$nontox[index])
        }
        la
      }
      crmht<-function(a,salpha,level,tox,nontox,scale){
        la<-a * exp(-a^2/2/scale^2) 
        for(i in sim1$level){
          index<-match(i,sim1$level)
          la<-la*(salpha[i]**exp(a)**sim1$tox[index])*((1-salpha[i]**exp(a))**sim1$nontox[index])
        }
        la
      }
      for(j in 1:length(alpha)){
        for(k in 1:s){			
          salpha<-alpha[[j]][k,]
          denpred[j,k]<-integrate(crmh,-Inf,Inf,salpha,sim1$level,sim1$tox,sim1$nontox,scale,abs.tol=0)[[1]]
          apred[j,k]<-integrate(crmht,-10,10,salpha,sim1$level,sim1$tox,sim1$nontox,scale,abs.tol=0)[[1]] / denpred[j,k]
        }
      }
      pord<-(denpred*prior.o)/sum(denpred*prior.o)
      if(method=="BMS"){
        ske<<-which(pord==max(pord),arr.ind=T)[1]
        ord<<-which(pord==max(pord),arr.ind=T)[2]
        ahat<<-apred[ske,ord]
        rpred<<-alpha[[ske]][ord,]**exp(ahat)
        next.lev<<-which.is.max(-(abs(rpred-theta)))
        next.lev
      }
      else if(method=="BMA"){
        spord<<-apply(pord,1,sum)
        opord<<-apply(pord,2,sum)
        ord<<-which.is.max(opord)
        ahat<<-apred[,ord]
        srpred<<-matrix(0,ncol=ncol(alpha[[1]]),nrow=length(alpha))
        for (i in 1:length(alpha)){
          srpred[i,]<<-alpha[[i]][ord,]**exp(ahat[i])
        }
        rpred<<-apply(srpred*spord,2,sum)
        next.lev<<-which.is.max(-(abs(rpred-theta)))
        next.lev
      }
      else{
        print("Incorrect parameter of method!")
      }
    }
    ###'crm' ENDS HERE
    ###LOAD FUNCTION 'twostgcrm1'
    #twostgcrm1函数：两阶段CRM框架下用BMS-POCRM推荐下一剂量
    twostgcrm1<-function(r,x0,stop,n,theta){
      n1<<-n+1
      obs<-data.frame(cbind(1:n1,rep(0,n1),rep(0,n1),rep(0,n1),rep(0,n1),rep(0,n1)))
      names(obs)<-c('patient','level','tox','a','skeleton','order')
      
      #######Beginning at dose level 1
      ###1st Stage:  Up and down scheme with cohort size 1
      i<-1
      #x0<-lapply(zones,ff)
      ##'initial.scheme' is a vector indicating the Stage I escalation scheme
      initial.scheme<-x0
      initial.scheme<-c(initial.scheme,rep(ncol(alpha[[1]]),n1-length(initial.scheme)))
      while(i < n1){
        obs$skeleton[i]<-0
        obs$order[i]<-0
        obs$level[1]<-initial.scheme[1]
        p<-runif(1)
        #number of tox in 1st patient
        index<-p<=r[obs$level[i]] ##determines any toxicities
        obs$tox[i]<-obs$tox[i]+as.numeric(index)
        if(any(obs$tox[1:i]==1) & any(obs$tox[1:i]==0)){
          q<-2
          break
        }
        if(all(obs$tox[1:i]==1)){
          i<-i+1
          obs$level[i]<-initial.scheme[1]
        }
        if(all(obs$tox[1:i]==0)){
          i<-i+1
          obs$level[i]<-initial.scheme[i]
        }
        
        if(length(obs$level[obs$level==d])==stop+1){
          MTD<-d
          break
        }
      }
      ##2nd stage
      N<<-table(obs$level>0)[2]+1
      if(any(obs$tox>0)){
        level<-crm(obs[1:(N-1),],alpha,prior.o,theta,method="BMS",scale)
        obs$a[N-1]<-ahat
        obs$skeleton[N-1]<-ske
        obs$order[N-1]<-ord
        for(j in N:n1){
          ##assigment for remaining patients
          obs$level[j]<-level
          if(obs$level[n1]>0){
            MTD<-obs$level[n1]
            break
          }
          if(length(obs$level[obs$level==level])==stop+1){
            MTD<-level
            break
          }
          index<-runif(1)<=r[obs$level[j]]
          if(index){obs$tox[j]<-1}
          level<-crm(obs[1:j,],alpha,prior.o,theta,method="BMS",scale)
          obs$a[j]<-ahat
          obs$skeleton[j]<-ske
          obs$order[j]<-ord
          ##crm dose recommendation for Nth patient
        }
      } else
        MTD<-d
      out<-list(trial=obs[obs$level>0,],MTD.selection=MTD)
    }
    ###'twostgcrm1' ENDS HERE
    ###LOAD FUNCTION 'twostgcrm2'
    #twostgcrm2函数：两阶段CRM框架下用BMA-POCRM推荐下一剂量
    twostgcrm2<-function(r,x0,stop,n,theta){
      
      n1<<-n+1
      obs<-data.frame(cbind(1:n1,rep(0,n1),rep(0,n1),rep(0,n1)))
      names(obs)<-c('patient','level','tox','order')
      #######Beginning at dose level 1
      ###1st Stage:  Up and down scheme with cohort size 1
      i<-1
      #x0<-lapply(zones,ff)
      ##'initial.scheme' is a vector indicating the Stage I escalation scheme
      initial.scheme<-x0
      initial.scheme<-c(initial.scheme,rep(ncol(alpha[[1]]),n1-length(initial.scheme)))
      while(i < n1){
        obs$order[i]<-0
        obs$level[1]<-initial.scheme[1]
        p<-runif(1)
        #number of tox in 1st patient
        index<-p<=r[obs$level[i]] ##determines any toxicities
        obs$tox[i]<-obs$tox[i]+as.numeric(index)
        if(any(obs$tox[1:i]==1) & any(obs$tox[1:i]==0)){
          q<-2
          break
        }
        if(all(obs$tox[1:i]==1)){
          i<-i+1
          obs$level[i]<-initial.scheme[1]
        }
        if(all(obs$tox[1:i]==0)){
          i<-i+1
          obs$level[i]<-initial.scheme[i]
        }
        
        if(length(obs$level[obs$level==d])==stop+1){
          MTD<-d
          break
        }
      }
      ##2nd stage
      N<<-table(obs$level>0)[2]+1
      if(any(obs$tox>0)){
        level<-crm(obs[1:(N-1),],alpha,prior.o,theta,method="BMA",scale)
        obs$order[N-1]<-ord
        for(j in N:n1){
          ##assigment for remaining patients
          obs$level[j]<-level
          if(obs$level[n1]>0){
            MTD<-obs$level[n1]
            break
          }
          
          if(length(obs$level[obs$level==level])==stop+1){
            MTD<-level
            break
          }
          index<-runif(1)<=r[obs$level[j]]
          if(index){obs$tox[j]<-1}
          level<-crm(obs[1:j,],alpha,prior.o,theta,method="BMA",scale)
          obs$order[j]<-ord
          ##crm dose recommendation for Nth patient
        }
      } else
        MTD<-d
      out<-list(trial=obs[obs$level>0,],MTD.selection=MTD)
    }
    ###'twostgcrm2' ENDS HERE
  }
  
  if(nsim>1){
    ###Load the function 'lpocrm' 
    #lpocrm函数：nsim>1时在两阶段CRM框架下用BMS/BMA-POCRM推荐下一剂量
    lpocrm<-function(r,alpha,prior.o,x0,stop,n,theta,method,scale){
      
      nord.tox = nrow(alpha[[1]])
      nske.tox = length(alpha)
      mprior.tox = prior.o
      ###Load the function 'bcrmh' 
      bcrmh<-function(a,p1,y,n,s){
        lik=exp(-a^2/2/s^2) 
        for(j in 1:length(p1)){
          lik=lik*((p1[j]**exp(a))**y[j])*((1-p1[j]**exp(a))**(n[j]-y[j]));
        }
        return(lik);
      }
      ###'bcrmh' ENDS HERE
      
      ###Load the function 'bcrmht' 
      bcrmht<-function(a,p1,y,n,s){
        lik=a * exp(-a^2/2/s^2) 
        for(j in 1:length(p1)){
          lik=lik*((p1[j]**exp(a))**y[j])*((1-p1[j]**exp(a))**(n[j]-y[j]));
        }
        return(lik);
      }
      ###'bcrmht' ENDS HERE
      ncomb = ncol(alpha[[1]])
      y=npts=ptox.hat=comb.select=numeric(ncomb)
      comb.curr = x0[1]
      stoprule=0
      i=1
      stage1<-c(x0,rep(ncol(alpha[[1]]),n-length(x0)))
      ##Stage 1
      while(i <= n){
        y[comb.curr] = y[comb.curr] + rbinom(1,1,r[comb.curr]);
        npts[comb.curr] = npts[comb.curr] + 1;
        
        if(sum(y)==sum(npts)){
          comb.curr<-ifelse(comb.curr==1,comb.curr,comb.curr-1)
        } else if(sum(y)==0){
          comb.curr<-ifelse(comb.curr==ncomb,comb.curr,stage1[i+1])
        } else {
          break
        }
        if(any(npts>stop)){
          stoprule<-0
          break
        }
        i=i+1
      }
      #Stage 2
      while(sum(npts) <= n)
      {
        if(sum(y)==0){
          stop=0
          break
        } else{
          like.tox= est.tox=matrix(0,nrow=nrow(prior.o),ncol=ncol(prior.o))
          for(m in 1:nske.tox)
          {
            for(k in 1:nord.tox)
            {
              salpha<-alpha[[m]][k,]
              like.tox[m,k]<-integrate(bcrmh,-Inf,Inf,p1=salpha,y=y,n=npts,s=scale,abs.tol=0)[[1]]
              est.tox[m,k]<-integrate(bcrmht,-10,10,p1=salpha,y=y,n=npts,s=scale,abs.tol=0)[[1]] /like.tox[m,k]
            }		
          }
          postprob.tox = (like.tox*mprior.tox)/sum(like.tox*mprior.tox)
          if(method=="BMS"){
            if(nske.tox>1){ 
              mtoxske.sel = which(postprob.tox==max(postprob.tox),arr.ind=T)[1]; 
            } else{
              mtoxske.sel = 1;
            }
            
            if(nord.tox>1){ 
              mtoxord.sel = which(postprob.tox==max(postprob.tox),arr.ind=T)[2]; 
            } else{
              mtoxord.sel = 1;
            }
            ptox.hat=alpha[[mtoxske.sel]][mtoxord.sel,]**exp(est.tox[mtoxske.sel,mtoxord.sel])
          }
          else if(method=="BMA"){
            spord.tox<-apply(postprob.tox,1,sum)
            opord.tox<-apply(postprob.tox,2,sum)
            mtoxord.sel<-which.is.max(opord.tox)
            a.hat<-est.tox[,mtoxord.sel]
            sptox.hat<-matrix(0,ncol=ncol(alpha[[1]]),nrow=length(alpha))
            for (l in 1:length(alpha)){
              sptox.hat[l,]<-alpha[[l]][mtoxord.sel,]**exp(a.hat[l])
            }
            ptox.hat<-apply(sptox.hat*spord.tox,2,sum)
          }
          else{
            print("Incorrect parameter of method!")
          }
          loss=abs(ptox.hat-theta)
          comb.curr=which.is.max(-loss)
          if(npts[comb.curr]==stop){
            stoprule<-0
            break
          }
          if(sum(npts)==n){
            stoprule=0
            break
          } else{
            y[comb.curr] = y[comb.curr] + rbinom(1,1,r[comb.curr]);
            npts[comb.curr] = npts[comb.curr] + 1;
          }
        }
      }
      if(stoprule==0){
        comb.select[comb.curr]=comb.select[comb.curr]+1;
      }
      return(list(MTD.selection=comb.select,tox.data=y,patient.allocation=npts))
    }
    ##########'lpocrm' ENDS HERE
  }
  ###Load the function 'lpocrm.sim' 
  ##lpocrm.sim函数：调用lpocrm函数做nsim次模拟
  lpocrm.sim<-function(nsim){
    ncomb=length(r)
    comb.select<-y<-npts<-matrix(nrow=nsim,ncol=ncomb)
    trialsize<-rep(0,nsim)
    nstop=0
    for(i in 1:nsim){
      result<-lpocrm(r,alpha,prior.o,x0,stop,n,theta,method,scale)
      comb.select[i,]=result$MTD.selection
      y[i,]=result$tox.data
      npts[i,]=result$patient.allocation
      trialsize[i]=sum(result$patient.allocation)
    }
    return(list(method=method,true.prob=r,MTD.selection=round(colMeans(comb.select),2),patient.allocation=round(colMeans(npts)/mean(trialsize),2),percent.DLT=sum(colMeans(y))/mean(trialsize),mean.n=mean(trialsize),acceptable=sum(colMeans(comb.select)[which(round(abs(r-theta),2)<=tox.range)])))
  }
  ##########'lpocrm.sim' end here
  if(nsim==1){
    if(method=="BMS"){
      twostgcrm1(r,x0,stop,n,theta)
    }
    else if(method=="BMA"){
      twostgcrm2(r,x0,stop,n,theta)
    }
    else{
      print("Incorrect parameter of method!")
    }
  } else{
    lpocrm.sim(nsim)
  }
}