optTrain = function(geno, cand, n.train, subpop=NULL, test=NULL,target=TRUE, method="rScore", min.iter=NULL){  
  if(!method%in%c("rScore","MSPE")){stop("Method not found. Please choose one from (rScore,MSPE)")}
  n=n.train; N=nrow(geno); Nc=length(cand)
  geno = as.matrix(geno)
  if(target==TRUE){
    subpop=as.character(subpop)
    if(method=="MSPE"){
      ## pop target MSPE
      if(is.null(min.iter)){min.iter=round((Nc-n)*n+1)}
      sol=matrix(NA, n, 30); 
      for(i in 1:30){
        pops = names(table(subpop))
        pop.ratio = round(n*as.numeric(table(subpop[test]))/sum(as.numeric(table(subpop[test]))))
        stop=0
        if(pop.ratio[which.min(pop.ratio)]==1){pop.ratio[which.min(pop.ratio)]=2}
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else if(sum(pop.ratio)<n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]+1
          }else{
            stop=1
          }
        }
        a=c()
        for(j in 1:length(pops)){
          a = c(a, sample(cand[subpop[cand]==pops[j]],pop.ratio[j]))
        }
        sol[,i] = sort(a)}
      score=rep(NA,30); for(i in 1:30){score[i]=mspe(geno[sol[,i],], geno[test,])}
      iter.score=score; top.score=min(score)
      stop=0; iter=1
      while(stop==0){
        elite = which(rank(1/score)>27)
        del = sample(seq(30)[-elite], 3, prob=((score[-elite])/sum(score[-elite])))
        for(i in 1:length(del)){
          par = sample(seq(30)[-del],2)
          b=c()
          for(j in 1:length(pops)){
            b=c(b,sample(unique(c(sol[,par[1]],sol[,par[2]]))[subpop[unique(c(sol[,par[1]],sol[,par[2]]))]==pops[j]],pop.ratio[j]))
            sol[,del[i]] = sort(b)}
        }
        for(i in 1:30){
          new.sol=sol[,i]
          p = sample(length(pops),1,prob=pop.ratio)
          if(i %in% del){
            sol[,i][sample(which(subpop[sol[,i]]==pops[p]),1)]=cand[!cand%in%sol[,i]][sample(which(subpop[cand[!cand%in%sol[,i]]]==pops[p]),1)]
            score[i] = mspe(geno[sol[,i],], geno[test,])
          }else{
            new.sol[sample(which(subpop[new.sol]==pops[p]),1)]=cand[!cand%in%new.sol][sample(which(subpop[cand[!cand%in%new.sol]]==pops[p]),1)]
            new.score = mspe(geno[new.sol,], geno[test,])
            if(new.score < score[i]){sol[,i] = new.sol; score[i]=new.score}
          }
        }
        iter.score=c(iter.score,mean(score)); top.score=c(top.score, min(score))
        cat(iter, "..", sep=""); iter = iter + 1
        if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-(Nc-n)*n])==0){stop=1}}
      }
      sol = sol[,which(score==min(score))[1]]
    }else if(method=="rScore"){
      ## pop target r-score
      if(is.null(min.iter)){min.iter=round((Nc-n)*n+1)}
      sol=matrix(NA, n, 30); 
      for(i in 1:30){
        pops = names(table(subpop))
        pop.ratio = round(n*as.numeric(table(subpop[test]))/sum(as.numeric(table(subpop[test]))))
        stop=0
        if(pop.ratio[which.min(pop.ratio)]==1){pop.ratio[which.min(pop.ratio)]=2}
        while(stop==0){
          if(sum(pop.ratio)>n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
          }else if(sum(pop.ratio)<n){
            pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]+1
          }else{
            stop=1
          }
        }
        a=c()
        for(j in 1:length(pops)){
          a = c(a, sample(cand[subpop[cand]==pops[j]],pop.ratio[j]))
        }
        sol[,i] = sort(a)}
      score=rep(NA,30); for(i in 1:30){score[i]=r_score(geno[sol[,i],], geno[test,])}
      iter.score=score; top.score=max(score)
      stop=0; iter=1
      while(stop==0){
        elite = which(rank(score)>27)
        del = sample(seq(30)[-elite], 3, prob=((1/score[-elite])/sum(1/score[-elite])))
        for(i in 1:length(del)){
          par = sample(seq(30)[-del],2)
          b=c()
          for(j in 1:length(pops)){
            b=c(b,sample(unique(c(sol[,par[1]],sol[,par[2]]))[subpop[unique(c(sol[,par[1]],sol[,par[2]]))]==pops[j]],pop.ratio[j]))
            sol[,del[i]] = sort(b)}
        }
        for(i in 1:30){
          new.sol=sol[,i]
          p = sample(length(pops),1,prob=pop.ratio)
          if(i %in% del){
            sol[,i][sample(which(subpop[sol[,i]]==pops[p]),1)]=cand[!cand%in%sol[,i]][sample(which(subpop[cand[!cand%in%sol[,i]]]==pops[p]),1)]
            score[i] = r_score(geno[sol[,i],], geno[test,])
          }else{
            new.sol[sample(which(subpop[sol[,i]]==pops[p]),1)]=cand[!cand%in%new.sol][sample(which(subpop[cand[!cand%in%new.sol]]==pops[p]),1)]
            new.score = r_score(geno[new.sol,], geno[test,])
            if(new.score > score[i]){sol[,i] = new.sol; score[i]=new.score}
          }
        }
        iter.score=c(iter.score,mean(score)); top.score=c(top.score, max(score))
        cat(iter, "..", sep=""); iter = iter + 1
        if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-(Nc-n)*n])==0){stop=1}}
      }
      sol = sol[,which(score==max(score))[1]]
    }}else{
      subpop=as.character(subpop)
      if(method=="MSPE"){
        ## pop untarget MSPE
        if(is.null(min.iter)){min.iter=round((Nc-n)*n+1)}
        sol=matrix(NA, n, 30); 
        for(i in 1:30){
          pops = names(table(subpop))
          pop.ratio = round(n*as.numeric(table(subpop[cand]))/sum(as.numeric(table(subpop[cand]))))
          stop=0
          if(pop.ratio[which.min(pop.ratio)]==1){pop.ratio[which.min(pop.ratio)]=2}
          while(stop==0){
            if(sum(pop.ratio)>n){
              pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
            }else if(sum(pop.ratio)<n){
              pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]+1
            }else{
              stop=1
            }
          }
          a=c()
          for(j in 1:length(pops)){
            a = c(a, sample(cand[subpop[cand]==pops[j]],pop.ratio[j]))
          }
          sol[,i] = sort(a)}
        score=rep(NA,30); for(i in 1:30){score[i]=mspe(geno[sol[,i],], geno[cand[!cand%in%sol[,i]],])}
        iter.score=score; top.score=min(score)
        stop=0; iter=1
        while(stop==0){
          elite = which(rank(1/score)>27)
          del = sample(seq(30)[-elite], 3, prob=((score[-elite])/sum(score[-elite])))
          for(i in 1:length(del)){
            par = sample(seq(30)[-del],2)
            b=c()
            for(j in 1:length(pops)){
              b=c(b,sample(unique(c(sol[,par[1]],sol[,par[2]]))[subpop[unique(c(sol[,par[1]],sol[,par[2]]))]==pops[j]],pop.ratio[j]))
              sol[,del[i]] = sort(b)}
          }
          for(i in 1:30){
            new.sol=sol[,i]
            p = sample(length(pops),1,prob=pop.ratio)
            if(i %in% del){
              sol[,i][sample(which(subpop[sol[,i]]==pops[p]),1)]=cand[!cand%in%sol[,i]][sample(which(subpop[cand[!cand%in%sol[,i]]]==pops[p]),1)]
              score[i] = mspe(geno[sol[,i],], geno[cand[!cand%in%sol[,i]],])
            }else{
              new.sol[sample(which(subpop[new.sol]==pops[p]),1)]=cand[!cand%in%new.sol][sample(which(subpop[cand[!cand%in%new.sol]]==pops[p]),1)]
              new.score = mspe(geno[new.sol,], geno[cand[!cand%in%new.sol],])
              if(new.score < score[i]){sol[,i] = new.sol; score[i]=new.score}
            }
          }
          iter.score=c(iter.score,mean(score)); top.score=c(top.score, min(score))
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-(Nc-n)*n])==0){stop=1}}
        }
        sol = sol[,which(score==min(score))[1]]
      }else if(method=="rScore"){
        ## pop untarget r-score
        if(is.null(min.iter)){min.iter=round((Nc-n)*n+1)}
        sol=matrix(NA, n, 30); 
        for(i in 1:30){
          pops = names(table(subpop))
          pop.ratio = round(n*as.numeric(table(subpop[cand]))/sum(as.numeric(table(subpop[cand]))))
          stop=0
          if(pop.ratio[which.min(pop.ratio)]==1){pop.ratio[which.min(pop.ratio)]=2}
          while(stop==0){
            if(sum(pop.ratio)>n){
              pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
            }else if(sum(pop.ratio)<n){
              pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]+1
            }else{
              stop=1
            }
          }
          a=c()
          for(j in 1:length(pops)){
            a = c(a, sample(cand[subpop[cand]==pops[j]],pop.ratio[j]))
          }
          sol[,i] = sort(a)}
        score=rep(NA,30); for(i in 1:30){score[i]=r_score(geno[sol[,i],], geno[cand[!cand%in%sol[,i]],])}
        iter.score=score; top.score=max(score)
        stop=0; iter=1
        while(stop==0){
          elite = which(rank(score)>27)
          del = sample(seq(30)[-elite], 3, prob=((1/score[-elite])/sum(1/score[-elite])))
          for(i in 1:length(del)){
            par = sample(seq(30)[-del],2)
            b=c()
            for(j in 1:length(pops)){
              b=c(b,sample(unique(c(sol[,par[1]],sol[,par[2]]))[subpop[unique(c(sol[,par[1]],sol[,par[2]]))]==pops[j]],pop.ratio[j]))
              sol[,del[i]] = sort(b)}
          }
          for(i in 1:30){
            new.sol=sol[,i]
            p = sample(length(pops),1,prob=pop.ratio)
            if(i %in% del){
              sol[,i][sample(which(subpop[sol[,i]]==pops[p]),1)]=cand[!cand%in%sol[,i]][sample(which(subpop[cand[!cand%in%sol[,i]]]==pops[p]),1)]
              score[i] = r_score(geno[sol[,i],], geno[cand[!cand%in%sol[,i]],])
            }else{
              new.sol[sample(which(subpop[sol[,i]]==pops[p]),1)]=cand[!cand%in%new.sol][sample(which(subpop[cand[!cand%in%new.sol]]==pops[p]),1)]
              new.score = r_score(geno[new.sol,], geno[cand[!cand%in%new.sol],])
              if(new.score > score[i]){sol[,i] = new.sol; score[i]=new.score}
            }
          }
          iter.score=c(iter.score,mean(score)); top.score=c(top.score, max(score))
          cat(iter, "..", sep=""); iter = iter + 1
          if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-(Nc-n)*n])==0){stop=1}}
        }
        sol = sol[,which(score==max(score))[1]]
      }}
  ret = list(
    OPTtrain=as.numeric(sort(sol)),
    TOPscore=as.numeric(top.score[-1]),
    ITERscore=as.numeric(iter.score[-1])
  )
  return(ret)
}

optTrain = function(geno, cand, n.train, subpop=NULL, test=NULL,target=TRUE, method="rScore", min.iter=NULL){  
  if(!method%in%c("rScore","MSPE")){stop("Method not found. Please choose one from (rScore,MSPE)")}
  n=n.train; N=nrow(geno); Nc=length(cand)
  geno = as.matrix(geno)
  if(target==TRUE){
    subpop=as.character(subpop)
    if(method=="MSPE"){
      ## pop target MSPE
      if(is.null(min.iter)){min.iter=round((Nc-n)*n+1)}
      pops = names(table(subpop))
      pop.ratio = round(n*as.numeric(table(subpop[test]))/sum(as.numeric(table(subpop[test]))))
      stop=0
      while(stop==0){
        if(sum(pop.ratio)>n){
          pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
        }else if(sum(pop.ratio)<n){
          pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]+1
        }else{
          stop=1
        }
      }
      sol=c()
      for(i in 1:length(pops)){
        sol = c(sol, sample(cand[subpop[cand]==pops[i]],pop.ratio[i]))
      }
      score=mspe(geno[sol,], geno[test,])
      iter.score=score; top.score=score
      stop=0; iter=1
      while(stop==0){
        new.sol=sol
        p = sample(length(pops),1,prob=pop.ratio)
        new.sol[sample(which(subpop[sol]==pops[p]),1)]=cand[!cand%in%new.sol][sample(which(subpop[cand[!cand%in%new.sol]]==pops[p]),1)]
        new.score=mspe(geno[new.sol,],geno[test,])
        if(new.score<score){
          sol=new.sol; score=new.score
          iter.score=c(iter.score,score)
          top.score=c(top.score,score)
        }else{
          iter.score=c(iter.score,new.score)
          top.score=c(top.score,score)
        }
        cat(iter, "..", sep=""); iter = iter + 1
        if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-(Nc-n)*n])==0){stop=1}}
        
      }
    }else if(method=="rScore"){
      ## pop target r-score
      if(is.null(min.iter)){min.iter=round((Nc-n)*n+1)}
      pops = names(table(subpop))
      pop.ratio = round(n*as.numeric(table(subpop[test]))/sum(as.numeric(table(subpop[test]))))
      stop=0
      while(stop==0){
        if(sum(pop.ratio)>n){
          pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
        }else if(sum(pop.ratio)<n){
          pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]+1
        }else{
          stop=1
        }
      }
      sol=c()
      for(i in 1:length(pops)){
        sol = c(sol, sample(cand[subpop[cand]==pops[i]],pop.ratio[i]))
      }
      score=r_score(geno[sol,], geno[test,])
      iter.score=score; top.score=score
      stop=0; iter=1
      while(stop==0){
        new.sol=sol
        p = sample(length(pops),1,prob=pop.ratio)
        new.sol[sample(which(subpop[sol]==pops[p]),1)]=cand[!cand%in%new.sol][sample(which(subpop[cand[!cand%in%new.sol]]==pops[p]),1)]
        new.score=r_score(geno[new.sol,],geno[test,])
        if(new.score>score){
          sol=new.sol; score=new.score
          iter.score=c(iter.score,score)
          top.score=c(top.score,score)
        }else{
          iter.score=c(iter.score,new.score)
          top.score=c(top.score,score)
        }
        cat(iter, "..", sep=""); iter = iter + 1
        if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-(Nc-n)*n])==0){stop=1}}
        
      }
    }
  }else{
    subpop=as.character(subpop)
    if(method=="MSPE"){
      ## pop untarget MSPE
      if(is.null(min.iter)){min.iter=round((Nc-n)*n+1)}
      pops = names(table(subpop))
      pop.ratio = round(n*as.numeric(table(subpop[cand]))/sum(as.numeric(table(subpop[cand]))))
      stop=0
      while(stop==0){
        if(sum(pop.ratio)>n){
          pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
        }else if(sum(pop.ratio)<n){
          pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]+1
        }else{
          stop=1
        }
      }
      sol=c()
      for(i in 1:length(pops)){
        sol = c(sol, sample(cand[subpop[cand]==pops[i]],pop.ratio[i]))
      }
      score=mspe(geno[sol,], geno[cand[!cand%in%sol],])
      iter.score=score; top.score=score
      stop=0; iter=1
      while(stop==0){
        new.sol=sol
        p = sample(length(pops),1,prob=pop.ratio)
        new.sol[sample(which(subpop[sol]==pops[p]),1)]=cand[!cand%in%new.sol][sample(which(subpop[cand[!cand%in%new.sol]]==pops[p]),1)]
        new.score=mspe(geno[new.sol,],geno[cand[!cand%in%new.sol],])
        if(new.score<score){
          sol=new.sol; score=new.score
          iter.score=c(iter.score,score)
          top.score=c(top.score,score)
        }else{
          iter.score=c(iter.score,new.score)
          top.score=c(top.score,score)
        }
        cat(iter, "..", sep=""); iter = iter + 1
        if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-(Nc-n)*n])==0){stop=1}}
        
      }
    }else if(method=="rScore"){
      ## pop untarget r-score
      if(is.null(min.iter)){min.iter=round((Nc-n)*n+1)}
      pops = names(table(subpop))
      pop.ratio = round(n*as.numeric(table(subpop[cand]))/sum(as.numeric(table(subpop[cand]))))
      stop=0
      while(stop==0){
        if(sum(pop.ratio)>n){
          pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]-1
        }else if(sum(pop.ratio)<n){
          pop.ratio[which(pop.ratio==max(pop.ratio))[1]]=pop.ratio[which(pop.ratio==max(pop.ratio))[1]]+1
        }else{
          stop=1
        }
      }
      sol=c()
      for(i in 1:length(pops)){
        sol = c(sol, sample(cand[subpop[cand]==pops[i]],pop.ratio[i]))
      }
      score=r_score(geno[sol,], geno[cand[!cand%in%sol],])
      iter.score=score; top.score=score
      stop=0; iter=1
      while(stop==0){
        new.sol=sol
        p = sample(length(pops),1,prob=pop.ratio)
        new.sol[sample(which(subpop[sol]==pops[p]),1)]=cand[!cand%in%new.sol][sample(which(subpop[cand[!cand%in%new.sol]]==pops[p]),1)]
        new.score=r_score(geno[new.sol,],geno[cand[!cand%in%new.sol],])
        if(new.score>score){
          sol=new.sol; score=new.score
          iter.score=c(iter.score,score)
          top.score=c(top.score,score)
        }else{
          iter.score=c(iter.score,new.score)
          top.score=c(top.score,score)
        }
        cat(iter, "..", sep=""); iter = iter + 1
        if(iter > min.iter){if(abs(top.score[iter]-top.score[iter-(Nc-n)*n])==0){stop=1}}
      }
    }
  }
  ret = list(
    OPTtrain=as.numeric(sort(sol)),
    TOPscore=as.numeric(top.score[-1]),
    ITERscore=as.numeric(iter.score[-1])
  )
  return(ret)
}

library(compiler)
r_score<-function(x,x0){
  n0=nrow(x0);n=nrow(x);I=diag(1,n0,n0);J=matrix(1/n0,n0,n0);
  A=t(x)%*%solve(x%*%t(x)+diag(1,n,n));X0AX=x0%*%A%*%x
  q1=(n0+1)+sum(diag(t(x0)%*%(I-J)%*%x0))
  q2=sum(diag(t(A)%*%t(x0)%*%(I+J)%*%x0%*%A))+sum(diag(t(X0AX)%*%(I-J)%*%X0AX))
  q12=sum(diag(t(x0)%*%(I-J)%*%X0AX))
  return(q12/sqrt(q1*q2))
}
mspe<-function(x,x0){
  n0=nrow(x0);n=nrow(x);I=diag(1,n0,n0);J=matrix(1/n0,n0,n0);
  A=t(x)%*%solve(x%*%t(x)+diag(1,n,n));X0AX=x0%*%A%*%x
  return(1+(1/n0)*(sum(diag(x0%*%A%*%t(A)%*%t(x0)))+sum(diag(t(x0-X0AX)%*%(x0-X0AX)))))
}
r_score<-cmpfun(r_score);mspe<-cmpfun(mspe)
library(foreach);library(doParallel)
cpu.cores<-detectCores()
cl=makeCluster(cpu.cores)
registerDoParallel(cl)

