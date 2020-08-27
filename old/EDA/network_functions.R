Heatmap=function(matrix, bound, filename=NA, numbers=FALSE, width=7, height=7,
                 show_colnames=FALSE, show_rownames=FALSE){
  Tmp=matrix
  Tmp[lower.tri(Tmp, diag = T)]=0
  if (numbers){
    TmpNumbers=round(matrix, digits=2)
    TmpNumbers[lower.tri(TmpNumbers, diag = T)]=""
    pheatmap(Tmp, cluster_rows = FALSE, cluster_cols = FALSE, width = width, height = height,
             breaks=seq(-bound,bound,length.out=100), border=NA,
             show_rownames = show_rownames, show_colnames = show_colnames, filename = filename, display_numbers = TmpNumbers,
             color=colorRampPalette(c('navy', 'blue', 'skyblue', 'white', 'gold', 'red', 'darkred'))(100))
  } else {
    pheatmap(Tmp, cluster_rows = FALSE, cluster_cols = FALSE, width = width, height = height,
             breaks=seq(-bound,bound,length.out=100), border=NA,
             show_rownames = show_rownames, show_colnames = show_colnames, filename = filename,
             color=colorRampPalette(c('navy', 'blue', 'skyblue', 'white', 'gold', 'red', 'darkred'))(100))
  }
}


LogtoLog10=function(x){
  return(x*log10(exp(1)))
}


EstimationSelection=function(data, y=NULL, q, implementation="pcor"){
  if (!implementation%in%c("cor", "pcor", "pcor.shrink")){
    stop("Argument implementation must be one of 'cor', 'pcor', 'pcor.shrink'.")
  }
  
  if (is.null(y)){
    if (implementation=="pcor"){
      metric=suppressWarnings(pcor(data)$estimate)
    }
    if (implementation=="pcor.shrink"){
      metric=pcor.shrink(data, verbose=FALSE)
    }
  } else {
    if (implementation=="pcor"){
      # print(data[y==0,])
      metric0=suppressWarnings(pcor(data[y==0,])$estimate)
      metric1=suppressWarnings(pcor(data[y==1,])$estimate)
    }
    if (implementation=="pcor.shrink"){
      # print(data[y==0,])
      # print(data[y==1,])
      metric0=pcor.shrink(data[y==0,], verbose=FALSE)
      metric1=pcor.shrink(data[y==1,], verbose=FALSE)
    }
    if (implementation=="cor"){
      metric0=cor(data[y==0,], verbose=FALSE)
      metric1=cor(data[y==1,], verbose=FALSE)
    }
    metric=metric1-metric0
  }
  
  rank_est=rank(1-abs(metric[upper.tri(metric)]))
  A=matrix(0, nrow=nrow(metric), ncol=nrow(metric))
  A[upper.tri(A)]=rank_est
  A=A+t(A)
  A=ifelse(A<=q,1,0)
  diag(A)=0
  colnames(A)=rownames(A)=colnames(data)
  
  return(A)
}


ComputePFER=function(q, pi, N){
  ## from MB
  return(1/(2*pi-1)*q^2/N)
}


GetSelectionProportions=function(data, y=NULL, lambda=NULL, q=NULL, tau=0.5, K=100, seed=1, folds=NULL, 
                                 method="GraphicalLasso", implementation="glasso"){
  # if folds is NULL, subsampling is performed with subsampling proportion of tau
  # if folds is set to e.g. 5, 5-fold is performed and repeated for K iterations. In this case, tau is ignored. 
  # method is in c("GraphicalLasso", "EstimationSelection")
  
  method=tolower(method)
  
  ### TO DO ###
  # Include stop/error messages if (i) lambda is not provided for GraphicalLasso, (ii) q is not provided for EstimationSelection
  # Include implementation for estimation-selection with choice between shrinkage, pseudo-inverse and others?
  # Find a way to re-include huge 
  
  if (!method%in%c("graphicallasso", "estimationselection")){
    stop("Argument method must be one of 'graphicallasso', 'estimationselection'.")
  }
  
  set.seed(seed)
  stab_iter=matrix(0,ncol=ncol(data),nrow=ncol(data))
  
  if (method=="graphicallasso"){ ###### TO DO 
    if (!implementation%in%c("glasso", "QUIC")){
      stop("For method='graphicallasso', implementation must be one of 'glasso', 'QUIC'.")
    }
    if (is.null(folds)){
      for (i in 1:K){
        s=sample(1:nrow(data), size=tau*nrow(data))
        data_sub=data[s,]
        cov_sub=cov(data_sub)
        if (implementation=="glasso"){
          g_sub=glasso(s=cov_sub, rho=lambda)
          omega=g_sub$wi
        }
        if (implementation=="QUIC"){ 
          g_sub=QUIC(S=cov_sub, rho=lambda, msg=0)
          omega=g_sub$X
        }
        # if (implementation=="huge"){ 
        #   g_sub=huge(x=cov_sub, lambda=lambda, method="glasso", verbose=FALSE)
        #   omega=g_sub$icov[[1]]
        # }
        A=ifelse(omega!=0, yes=1, no=0)
        A=A+t(A)
        A=ifelse(A!=0, yes=1, no=0)
        stab_iter=stab_iter+A
      }
      stab_iter=stab_iter/K
      diag(stab_iter)=0
    } else {
      for (i in 1:K){
        kfolds=split(sample(1:nrow(data), size=nrow(data)), f=ceiling(seq(10^-16,folds,length.out=nrow(data))))
        for (k in 1:folds){
          s=kfolds[[k]]
          data_sub=data[s,]
          cov_sub=cov(data_sub)
          if (implementation=="glasso"){
            g_sub=glasso(s=cov_sub, rho=lambda)
            omega=g_sub$wi
          }
          if (implementation=="QUIC"){ 
            g_sub=QUIC(S=cov_sub, rho=lambda, msg=0)
            omega=g_sub$X
          }
          A=ifelse(omega!=0, yes=1, no=0)
          A=A+t(A)
          A=ifelse(A!=0, yes=1, no=0)
          stab_iter=stab_iter+A
        }
      }
      stab_iter=stab_iter/(K*folds)
      diag(stab_iter)=0
    }
  }
  
  if (method=="estimationselection"){
    if (!implementation%in%c("cor", "pcor", "pcor.shrink")){
      stop("For method='estimationselection', implementation must be one of 'cor', 'pcor', 'pcor.shrink'.")
    }
    if (is.null(folds)){
      for (i in 1:K){
        s=sample(1:nrow(data), size=tau*nrow(data))
        data_sub=data[s,]
        y_sub=y[s]
        A=EstimationSelection(data_sub, y=y_sub, q, implementation=implementation)
        stab_iter=stab_iter+A
      }
      stab_iter=stab_iter/K
      diag(stab_iter)=0
    } else {
      for (i in 1:K){
        kfolds=split(sample(1:nrow(data), size=nrow(data)), f=ceiling(seq(10^-16,folds,length.out=nrow(data))))
        for (k in 1:folds){
          s=kfolds[[k]]
          data_sub=data[s,]
          y_sub=y[s]
          A=EstimationSelection(data_sub, y=y_sub, q, implementation=implementation)
          stab_iter=stab_iter+A
        }
      }
      stab_iter=stab_iter/(K*folds)
      diag(stab_iter)=0
    }
  }
  
  return(stab_iter)
}


# GetSelectionProportions=function(data, lambda=NULL, q=NULL, tau=0.5, K=100, seed=1, folds=NULL, 
#                                  method="GraphicalLasso", implementation="glasso"){
#   # if folds is NULL, subsampling is performed with subsampling proportion of tau
#   # if folds is set to e.g. 5, 5-fold is performed and repeated for K iterations. In this case, tau is ignored. 
#   # method is in c("GraphicalLasso", "EstimationSelection")
#   
#   method=tolower(method)
#   
#   ### TO DO ###
#   # Include stop/error messages if (i) lambda is not provided for GraphicalLasso, (ii) q is not provided for EstimationSelection
#   # Include implementation for estimation-selection with choice between shrinkage, pseudo-inverse and others?
#   # Find a way to re-include huge 
#   
#   
#   if (!method%in%c("graphicallasso", "estimationselection")){
#     stop("Argument method must be one of 'graphicallasso', 'estimationselection'.")
#   }
#   
#   set.seed(seed)
#   stab_iter=matrix(0,ncol=ncol(data),nrow=ncol(data))
#   
#   if (method=="graphicallasso"){
#     if (!implementation%in%c("glasso", "QUIC")){
#       stop("For method='graphicallasso', implementation must be one of 'glasso', 'QUIC'.")
#     }
#     if (is.null(folds)){
#       for (i in 1:K){
#         s=sample(1:nrow(data), size=tau*nrow(data))
#         data_sub=data[s,]
#         cov_sub=cov(data_sub)
#         if (implementation=="glasso"){
#           g_sub=glasso(s=cov_sub, rho=lambda)
#           omega=g_sub$wi
#         }
#         if (implementation=="QUIC"){ 
#           g_sub=QUIC(S=cov_sub, rho=lambda, msg=0)
#           omega=g_sub$X
#         }
#         # if (implementation=="huge"){ 
#         #   g_sub=huge(x=cov_sub, lambda=lambda, method="glasso", verbose=FALSE)
#         #   omega=g_sub$icov[[1]]
#         # }
#         A=ifelse(omega!=0, yes=1, no=0)
#         A=A+t(A)
#         A=ifelse(A!=0, yes=1, no=0)
#         stab_iter=stab_iter+A
#       }
#       stab_iter=stab_iter/K
#       diag(stab_iter)=0
#     } else {
#       for (i in 1:K){
#         kfolds=split(sample(1:nrow(data), size=nrow(data)), f=ceiling(seq(10^-16,folds,length.out=nrow(data))))
#         for (k in 1:folds){
#           s=kfolds[[k]]
#           data_sub=data[s,]
#           cov_sub=cov(data_sub)
#           if (implementation=="glasso"){
#             g_sub=glasso(s=cov_sub, rho=lambda)
#             omega=g_sub$wi
#           }
#           if (implementation=="QUIC"){ 
#             g_sub=QUIC(S=cov_sub, rho=lambda, msg=0)
#             omega=g_sub$X
#           }
#           # if (implementation=="huge"){ 
#           #   g_sub=huge(x=cov_sub, lambda=lambda, method="glasso", verbose=FALSE)
#           #   omega=g_sub$icov[[1]]
#           # }
#           A=ifelse(omega!=0, yes=1, no=0)
#           A=A+t(A)
#           A=ifelse(A!=0, yes=1, no=0)
#           stab_iter=stab_iter+A
#         }
#       }
#       stab_iter=stab_iter/(K*folds)
#       diag(stab_iter)=0
#     }
#   }
#   
#   if (method=="estimationselection"){
#     if (!implementation%in%c("cor", "pcor", "pcor.shrink")){
#       stop("For method='estimationselection', implementation must be one of 'cor', 'pcor', 'pcor.shrink'.")
#     }
#     if (is.null(folds)){
#       for (i in 1:K){
#         s=sample(1:nrow(data), size=tau*nrow(data))
#         data_sub=data[s,]
#         A=EstimationSelection(data_sub, q, implementation=implementation)
#         stab_iter=stab_iter+A
#       }
#       stab_iter=stab_iter/K
#       diag(stab_iter)=0
#     } else {
#       for (i in 1:K){
#         kfolds=split(sample(1:nrow(data), size=nrow(data)), f=ceiling(seq(10^-16,folds,length.out=nrow(data))))
#         for (k in 1:folds){
#           s=kfolds[[k]]
#           data_sub=data[s,]
#           A=EstimationSelection(data_sub, q, implementation=implementation)
#           stab_iter=stab_iter+A
#         }
#       }
#       stab_iter=stab_iter/(K*folds)
#       diag(stab_iter)=0
#     }
#   }
#   
#   return(stab_iter)
# }


GetBinomialProbabilities=function(q, N, pi, K){
  p_1=pbinom(K*pi, size=K, prob=q/N, log.p=TRUE) # proportion <= pi
  p_2=log(pbinom(K*(1-pi)-1, size=K, prob=q/N)-pbinom(K*pi, size=K, prob=q/N)) # pi < proportion < (1-pi)
  if (is.infinite(p_2)){
    p_2=0
    for (i in seq(K*(1-pi)-1, K*pi+1)){
      p_2=p_2+dbinom(round(i), size=K, prob=q/N)
    }
    p_2=log(p_2)
  }
  p_3=pbinom(K*(1-pi)-1, size=K, prob=q/N, lower.tail=FALSE, log.p=TRUE) # proportion >= pi 
  # p2=log(1-exp(p_1)-exp(p_3))
  
  if (abs(exp(p_1)+exp(p_2)+exp(p_3)-1)>1e-10){
    stop(paste0("Probabilities do not sum to 1 (Binomial distribution) \n p_1+p_2+p_3=", exp(p_1)+exp(p_2)+exp(p_3)))
  }
  return(list(p_1=p_1, p_2=p_2, p_3=p_3))
}


GetBinomialScore=function(stab_iter, q=NULL, N=NULL, pi, K){
  if (is.matrix(stab_iter)){
    stab_iter=stab_iter[upper.tri(stab_iter)]
  }
  
  if (is.null(q)){
    q=round(sum(stab_iter))
  }
  
  # if (is.null(N)){
  #   p=ncol(stab_iter)
  #   N=p*(p-1)/2
  # }
  
  p_vect=GetBinomialProbabilities(q, N, pi, K)
  
  S_0=sum(stab_iter<=pi)
  S_1=sum(stab_iter>=(1-pi))
  U=sum((stab_iter<(1-pi))&(stab_iter>pi))
  
  if (S_0+S_1+U!=N){
    stop(paste0("Inconsistency in number of edges \n S_0+S_1+U=", S_0+S_1+U, " instead of ", N))
  }
  
  l=S_0*p_vect$p_1+U*p_vect$p_2+S_1*p_vect$p_3
  
  if (is.infinite(l)){
    l=NA
  }
  
  # if (S_1==0){
  #   l=NA # to prevent from having EMPTY networks
  # }
  
  return(l)
}


GetPFERPairs=function(q_list, pi_list, K, N){
  if (is.null(q_list)){
    q_list=1:N
  }
  
  PFER=matrix(NA, nrow=length(q_list), ncol=length(pi_list))
  for (k in 1:length(q_list)){
    q=q_list[k]
    for (j in 1:length(pi_list)){
      pi=pi_list[j]
      # PFER[k,j]=stabsel_parameters(p=N, cutoff=pi, q=q, B=K, sampling.type="MB")$PFER
      PFER[k,j]=ComputePFER(q=q,pi=pi,N=N)
    }
  }
  
  return(PFER)
}


# PFER=GetPFERPairs(q_list=1:N, pi_list, K=100, N)
# pheatmap(PFER, cluster_rows = FALSE, cluster_cols = FALSE)


GetBlockMatrix=function(pk){
  nblocks=sum(upper.tri(matrix(NA, ncol=length(pk), nrow=length(pk)), diag=TRUE))
  blocks=matrix(NA, nrow=length(pk), ncol=length(pk))
  blocks[upper.tri(blocks, diag=TRUE)]=1:nblocks
  
  mybreaks=c(0,cumsum(pk))
  bigblocks=matrix(ncol=sum(pk), nrow=sum(pk))
  row_id_start=matrix(mybreaks[row(blocks)],ncol=length(pk))+1
  row_id_end=matrix(mybreaks[row(blocks)+1],ncol=length(pk))
  col_id_start=matrix(mybreaks[col(blocks)],ncol=length(pk))+1
  col_id_end=matrix(mybreaks[col(blocks)+1],ncol=length(pk))
  
  row_id_start=row_id_start[upper.tri(row_id_start, diag=TRUE)]
  row_id_end=row_id_end[upper.tri(row_id_end, diag=TRUE)]
  col_id_start=col_id_start[upper.tri(col_id_start, diag=TRUE)]
  col_id_end=col_id_end[upper.tri(col_id_end, diag=TRUE)]
  
  for (block_id in blocks[upper.tri(blocks, diag=TRUE)]){
    ids=rbind(expand.grid(row_id_start[block_id]:row_id_end[block_id], 
                          col_id_start[block_id]:col_id_end[block_id]),
              expand.grid(col_id_start[block_id]:col_id_end[block_id], 
                          row_id_start[block_id]:row_id_end[block_id]))
    bigblocks[as.matrix(ids)]=block_id
  }
  
  return(bigblocks)
}


# stab=NULL; data=NULL; pi_list=seq(0.6,0.9,by=0.01); cov=NULL;
# K=100; tau=0.5; seed=1; Lambda=NULL; Q=NULL; PFER_thr=Inf;
# method="EstimationSelection"; pk=NULL; implementation="pcor"; blocks=NULL
# data=X$X
# mycov=cov(X$X)
# lmax=getMaxCov(mycov)
# Q=seq(1,N,length.out=N_it)
# calib_method="grid"
# # Lambda=getLamPath(lmax,lmax*0.001,len=N_it,log=TRUE)
# K=10
# # pk=c(5,15)


BinomCalib=function(stab=NULL, data=NULL, y=NULL, pi_list=seq(0.6,0.9,by=0.01), cov=NULL,
                    K=100, B=NULL, tau=0.5, seed=1, Lambda=NULL, Q=NULL, PFER_thr=Inf, 
                    method="GraphicalLasso", pk=NULL, implementation=NULL, calib_method="grid", blocks=NULL){
  # tau should be 0.5 as advised by stabs 
  # If PFER_thr is Inf, problem is unconstrained, otherwise it is constrained
  # method is in c("GraphicalLasso", "EstimationSelection")
  # implementation is in c("glasso", "QUIC") for the Graphical LASSO
  # and in c("pcor", "pcor.shrink") for estimation selection
  # blocks: block ids as given by GetBlockMatrix() of the block to calibrate
  # (e.g. can be used to parallelise M-O calibration or if the user is only interested in one of the blocks, 
  # like only bipartite edges, or only transcripts but conditionally on the other variables)
  # if blocks is not specified, the default is to calibrate all blocks
  
  ### TO DO: ###
  # If stab is provided, K, tau, seed, Lambda, bootstrap, data are ignored
  # Lambda has to be provided if stab is not
  # If stab is not provided, the subsampling procedure is run (data, K, tau, seed, Lambda and bootstrap are used). 
  # Otherwise selection proportions from stab are used.
  
  method=tolower(method)
  
  if (is.null(implementation)){
    if (method=="graphicallasso"){
      implementation="glasso"
      Q=NULL
    } 
    if (method=="estimationselection"){
      implementation="pcor.shrink"
    }
  }
  
  p=ncol(data)
  N=p*(p-1)/2
  
  if (method=="graphicallasso"){
    if (is.null(Lambda)){
      stop("Lambda must be provided for the graphical lasso")
    }
    if (!is.null(stab)){
      stab$selprop=stab$selprop[,,sort.list(stab$Lambda, decreasing=TRUE)]
      stab$Lambda=sort(stab$Lambda, decreasing=TRUE)
      Lambda=stab$Lambda
    } else {
      if (length(dim(Lambda))==0){
        Lambda=matrix(Lambda,ncol=1)
      }
      Lambda=apply(Lambda, 2, sort, decreasing=TRUE)
      bigstab=array(NA,dim=c(ncol(data),ncol(data),nrow(Lambda)),dimnames=list(colnames(data),colnames(data),NULL))
    }
    if (is.null(cov)){
      cov=cov(data)
    }
  }
  
  if (is.null(pk)){
    pk=p
  } 
  
  nblocks=sum(upper.tri(matrix(NA, ncol=length(pk), nrow=length(pk)), diag=TRUE))
  bigblocks=GetBlockMatrix(pk)
  if (length(PFER_thr)==1){
    PFER_thr_blocks=ceiling(PFER_thr*table(bigblocks[upper.tri(bigblocks)])/N)
  } else {
    if (length(PFER_thr)!=nblocks){
      stop(paste0("Please provide a single value for argument PFER_thr or a vector with as many values as there are blocks (i.e.", nblocks," in this case)."))
    }
    PFER_thr_blocks=ceiling(PFER_thr)
  }
  
  if (method=="estimationselection"){
    Lambda=NULL
    if (is.null(Q)){
      stop("Q must be provided for Estimation-Selection")
    }
    p=ncol(data)
    N=p*(p-1)/2
    if (length(dim(Q))==0){
      Q=matrix(Q,ncol=1)
    }
    if ((ncol(Q)!=nblocks)&(ncol(Q)!=1)){
      stop(paste0("Matrix Q must have as many columns as there are blocks (i.e. ", nblocks, " in this instance)"))
    }
    if ((ncol(Q)==1)&(nblocks>1)){
      Q=matrix(rep(Q,nblocks), ncol=nblocks)
    }
    bigstab=array(NA,dim=c(ncol(data),ncol(data),nrow(Q)),dimnames=list(colnames(data),colnames(data),NULL))
  }
  
  if (method=="graphicallasso"){
    Q=matrix(NA,nrow=nrow(Lambda),ncol=nblocks)
    if ((ncol(Lambda)!=nblocks)&(ncol(Lambda)!=1)){
      stop(paste0("Matrix Lambda must have as many columns as there are blocks (i.e. ", nblocks, " in this instance)"))
    }
    if ((ncol(Lambda)==1)&(nblocks>1)){
      Lambda=matrix(rep(Lambda, nblocks),ncol=nblocks)
    }
  }
  
  if (is.null(blocks)){
    blocks=1:nblocks
  }
  
  for (k in blocks){
    assign(paste0("loglik", k), matrix(NA,nrow=nrow(Q),ncol=length(pi_list)))
    assign(paste0("PFERs", k), matrix(NA,nrow=nrow(Q),ncol=length(pi_list)))
  }
  
  if (K>1){
    S=P=matrix(NA, nrow=nrow(Q), ncol=nblocks)
  }
  
  N=rep(NA, nblocks)
  
  for (block_id in blocks){
    # print(block_id)
    k=min_PFER=previous_S=max_S=0
    while ((k<nrow(Q))&(min_PFER<=PFER_thr_blocks[block_id])&(max_S>=previous_S)){
      k=k+1
      previous_S=max_S
      # print(k)
      
      if (method=="graphicallasso"){
        lambda=Lambda[k,block_id]
        lambdamat=bigblocks
        lambdamat[bigblocks==block_id]=lambda
        lambdamat[bigblocks!=block_id]=max(Lambda[,block_id])*10
        
        if (is.null(stab)){
          stab_iter=GetSelectionProportions(data=data, y=y, lambda=lambdamat, tau=tau, K=K, seed=seed, method=method, implementation=implementation)
          bigstab[cbind(which(bigblocks==block_id, arr.ind = TRUE), k)]=stab_iter[which(bigblocks==block_id, arr.ind = TRUE)]
        } else {
          stab_iter=stab$selprop[,,k]
        }
        q=round(sum(stab_iter[upper.tri(stab_iter)])) # estimate of the average number of selected edges
        Q[k,block_id]=q
      }
      
      if (method=="estimationselection"){ #### TO DO: block-calibration with qmat
        q=Q[k,block_id]
        if (is.null(stab)){
          stab_iter=GetSelectionProportions(data=data, y=y, q=q, tau=tau, K=K, seed=seed, method=method, implementation=implementation)
          bigstab[,,k]=stab_iter
        } else {
          stab_iter=stab$selprop[,,k]
        }
      }
      
      stab_iter_block=stab_iter[(bigblocks==block_id)&(upper.tri(bigblocks))]
      q_block=round(sum(stab_iter_block))
      N_block=length(stab_iter_block)
      N[block_id]=N_block
      
      loglik_tmp=eval(parse(text=paste0("loglik",block_id)))
      PFERs_tmp=eval(parse(text=paste0("PFERs",block_id)))
      for (j in 1:length(pi_list)){
        pi=1-pi_list[j]
        PFERs_tmp[k,j]=ComputePFER(q=q_block,pi=1-pi,N=N_block)
        if ((PFERs_tmp[k,j]<=PFER_thr_blocks[block_id])&(K>1)){
          loglik_tmp[k,j]=GetBinomialScore(stab_iter_block, q=q_block, N=N_block, pi=pi, K=K)
        }
      }
      min_PFER=min(PFERs_tmp[k,], na.rm=TRUE)
      if (calib_method=="localmax"){
        if (sum(!is.na(loglik_tmp[k,]))>0){
          max_S=max(-loglik_tmp[k,], na.rm=TRUE)
        }
        print(max_S)
      }
      assign(paste0("PFERs",block_id), PFERs_tmp)
      
      if (K>1) {
        assign(paste0("loglik",block_id), loglik_tmp)
        
        if (sum(is.na(loglik_tmp[k,]))==ncol(loglik_tmp)){
          hat_pi=min(pi_list)
        } else {
          hat_pi=pi_list[which.min(loglik_tmp[k,])]
        }
        
        P[k,block_id]=hat_pi
        S[k,block_id]=sum(stab_iter[upper.tri(stab_iter)]>=hat_pi)
      }
    }
  }
  
  if (nblocks==1){
    PFERs=PFERs1
  } else {
    PFERs=list()
    for (block_id in blocks){
      PFERs=c(PFERs, list(eval(parse(text=paste0("PFERs",block_id,"=PFERs",block_id)))))
    }
  }
  
  if (K==1){ ## Screening
    return(list(Q=Q, PFERs=PFERs, N=N, seed=seed, K=K, tau=tau, Lambda=Lambda))
  } else {
    if (is.null(stab)){
      stab=list(seed=seed, K=K, tau=tau, Lambda=Lambda, selprop=bigstab)
    }
    
    if (nblocks==1){
      loglik=loglik1
    } else {
      loglik=list()
      for (block_id in blocks){
        loglik=c(loglik, list(eval(parse(text=paste0("loglik",block_id,"=loglik",block_id)))))
      }
    }
    
    return(list(loglik=loglik, N=N,
                PFER=PFERs, Q=Q, S=S, P=P,
                pi_list=pi_list,
                stab=stab))
  } 
}


# stab=NULL; data=NULL; pi_list=seq(0.6,0.9,by=0.01); cov=NULL;
# K=100; tau=0.5; seed=1; Lambda=NULL; Q=NULL; PFER_thr=10;
# method="GraphicalLasso"; pk=NULL; implementation="glasso"; blocks=NULL
# calib_method="grid"
# data=X$X
# delta=+Inf
# param_card1=100; param_card2=10; myfactor=0.001


Calibrate=function(data=NULL, y=NULL, pi_list=seq(0.6,0.9,by=0.01), 
                   K=100, B=NULL, tau=0.5, seed=1, PFER_thr=10, 
                   method="GraphicalLasso", pk=NULL, implementation=NULL, calib_method="grid", blocks=NULL, 
                   param_card1=100, param_card2=10, myfactor=0.0001, 
                   factor_range=0.7, delta=+Inf){
  # delta=+Inf to perform full calibration, delta=0 for screening only 
  
  ## TO DO: for estimation-selection 
  ## TO DO: for M-O
  
  if (!is.null(pk)){
    if (sum(pk)!=ncol(data)){
      stop("Argument pk is not consistent with the number of variables in data. Please make sure that sum(pk) is equal to ncol(data).")
    }
  }
  
  method=tolower(method)
  
  print("PFER threshold(s):")
  print(PFER_thr)
  print(paste("Cardinal 1:", param_card1))
  print(paste("Cardinal 2:", param_card2))
  
  # if (method=="graphicallasso"){
  #   mycov=cov(data)
  #   lmax=getMaxCov(mycov)
  # } 
  p=ncol(data)
  N=p*(p-1)/2
  
  if (is.null(pk)){
    pk=p
  }
  
  nblocks=sum(upper.tri(matrix(NA, ncol=length(pk), nrow=length(pk)), diag=TRUE))
  bigblocks=GetBlockMatrix(pk)
  
  if (length(PFER_thr)==1){
    PFER_thr_blocks=ceiling(PFER_thr*table(bigblocks[upper.tri(bigblocks)])/N)
  } else {
    if (length(PFER_thr)!=nblocks){
      stop(paste0("Please provide a single value for argument PFER_thr or a vector with as many values as there are blocks (i.e.", nblocks," in this case)."))
    }
    PFER_thr_blocks=ceiling(PFER_thr)
  }
  
  if (method=="graphicallasso"){
    print(paste("Factor:", myfactor))
    print("Defining first set Lambda...")
    mycov=cov(data)
    lmax=getMaxCov(mycov)
    if (is.infinite(PFER_thr[1])){
      max_q=0
      while (max_q<N) {
        # myfactor=myfactor*0.001
        Lambda=getLamPath(lmax,lmax*myfactor,len=param_card1,log=TRUE)
        myscreen=BinomCalib(stab=NULL, data=data, y=y, pi_list=1, cov=NULL, K=1, tau=tau, seed=seed, Lambda=Lambda, Q=NULL, 
                            PFER_thr=PFER_thr, method=method, pk=pk, implementation=implementation, blocks=blocks, calib_method=calib_method)
        max_q=max(myscreen$Q, na.rm=TRUE)
      }
    } else {
      man_PFER=0
      while (all(man_PFER<PFER_thr)){
        Lambda=getLamPath(lmax,lmax*myfactor,len=param_card1,log=TRUE)
        myscreen=BinomCalib(stab=NULL, data=data, y=y, pi_list=1, cov=NULL, K=1, tau=tau, seed=seed, Lambda=Lambda, Q=NULL, 
                            PFER_thr=PFER_thr, method=method, pk=pk, implementation=implementation, blocks=blocks, calib_method=calib_method)
        if (length(pk)==1){
          tmpPFER=myscreen$PFERs
          tmpPFER=tmpPFER[apply(tmpPFER, 1, FUN=function(x){sum(!is.na(x))})>0,,drop=FALSE]
          man_PFER=max(apply(tmpPFER, 1, min, na.rm=TRUE), na.rm=TRUE)
        } else {
          man_PFER=NULL
          for (k in 1:nblocks){
            tmpPFER=myscreen$PFER[[k]]
            tmpPFER=tmpPFER[apply(tmpPFER, 1, FUN=function(x){sum(!is.na(x))})>0,,drop=FALSE]
            man_PFER=c(man_PFER,max(apply(tmpPFER, 1, min, na.rm=TRUE), na.rm=TRUE))
          }
        }
      }
    }
    print("Corresponding Q:")
    # print(myscreen$Q[!is.na(myscreen$Q)])
    print(myscreen$Q[apply(myscreen$Q, 1, FUN=function(x){sum(!is.na(x))})>=1,])
    Lambda=matrix(NA, ncol=nblocks, nrow=param_card2)
    for (k in 1:nblocks){
      myrange=range(myscreen$Lambda[!is.na(myscreen$Q[,k]),k])
      Lambda[,k]=getLamPath(max(myrange),min(myrange)*factor_range,len=param_card2,log=TRUE) # *0.7 to be less stringent on PFER in case
    }
    # myrange=range(myscreen$Lambda[apply(myscreen$Q, 1, FUN=function(x){sum(!is.na(x))})>1,,drop=FALSE])
    # Lambda=getLamPath(max(myrange),min(myrange)*factor_range,len=param_card2,log=TRUE) # *0.7 to be less stringent on PFER in case
    Q=NULL
    if (delta==0){
      out_full=list(method=method, Lambda=Lambda)
    } 
  } 
  
  ### TO DO: check that it works for M-O network with estimation-selection
  if (method=="estimationselection"){
    Q=seq(1,N,length.out=param_card1)
    tmpPFER=GetPFERPairs(q_list=Q, pi_list=1, K=K, N=N)
    Q=round(seq(1, max(Q[tmpPFER<PFER_thr]), length.out=param_card2))
    Q=matrix(rep(Q,nblocks), ncol=nblocks)
    if (delta==0){
      out_full=list(method=method, Q=Q)
    } 
  }
  
  k=0
  while (TRUE%in%(delta>1)){
    k=k+1
    print(paste("Refining the grid of parameters... Iteration", k))
    out=BinomCalib(stab=NULL, data=data, y=y, pi_list=pi_list, cov=NULL, K=K, tau=tau, seed=seed, Lambda=Lambda, Q=Q, 
                   PFER_thr=PFER_thr, method=method, pk=pk, implementation=implementation, blocks=blocks, calib_method=calib_method)
    Q=out$Q
    tmpQ=out$Q[apply(out$Q, 1, FUN=function(x){sum(!is.na(x))})>=1,,drop=FALSE]
    print("Corresponding Q:")
    if (nblocks==1){
      print(as.vector(tmpQ))
    } else {  
      print(tmpQ)
    }
    delta=apply(tmpQ, 2, FUN=function(x){round((max(x,na.rm=TRUE)-min(x,na.rm=TRUE))/length(x))})
    if (nblocks==1){
      myS=-out$loglik
      myS[is.na(myS)]=0
      id=which.max(apply(myS, 1, max, na.rm=TRUE))
      if (method=="graphicallasso"){
        lmax=Lambda[ifelse((id-1)>=1, yes=id-1, no=1)]
        lmin=Lambda[ifelse((id+1)<=length(Lambda), yes=id+1, no=length(Lambda))]
        Lambda=getLamPath(lmax,lmin,len=param_card2,log=TRUE)
      }
      if (method=="estimationselection"){
        qmin=Q[ifelse((id-1)>=1, yes=id-1, no=1)]
        qmax=Q[ifelse((id+1)<=length(Q), yes=id+1, no=length(Q))]
        Q=unique(round(seq(qmin, qmax, length.out=param_card2)))
      }
    } else {
      for (block_id in 1:nblocks){
        myS=-out$loglik[[block_id]]
        myS[is.na(myS)]=0
        id=which.max(apply(myS, 1, max, na.rm=TRUE))
        if (method=="graphicallasso"){
          lmax=Lambda[ifelse((id-1)>=1, yes=id-1, no=1),block_id]
          lmin=Lambda[ifelse((id+1)<=nrow(Lambda), yes=id+1, no=nrow(Lambda)),block_id]
          Lambda[,block_id]=getLamPath(lmax,lmin,len=param_card2,log=TRUE)
        }
        if (method=="estimationselection"){
          qmin=Q[ifelse((id-1)>=1, yes=id-1, no=1),block_id]
          qmax=Q[ifelse((id+1)<=nrow(Q), yes=id+1, no=nrow(Q)),block_id]
          Q[,block_id]=(round(seq(qmin, qmax, length.out=param_card2)))
        }
      }
    }
    if (k==1){
      out_full=out
    } else {
      out_full=MergeOutput(out_full, out)
    }
  }
  
  if (ncol(out_full$Q)==1){
  rownames(out_full$loglik)=out_full$Q
  colnames(out_full$loglik)=out_full$pi_list
  } else {
    for (k in 1:ncol(out_full$Q)){
      rownames(out_full$loglik[[k]])=out_full$Q[,k]
      colnames(out_full$loglik[[k]])=out_full$pi_list
    }
  }
  return(out_full)
}


Screen=function(stab=NULL, data=NULL, pi_list=seq(0.6,0.9,by=0.01), 
                K=100, B=NULL, tau=0.5, seed=1, PFER_thr=10, 
                method="GraphicalLasso", pk=NULL, implementation=NULL, calib_method="grid", blocks=NULL, 
                param_card1=100, param_card2=10, myfactor=0.001, 
                factor_range=0.7, delta=+Inf){
  return(Calibrate(data=data, pi_list=pi_list, 
                   K=K, B=B, tau=tau, seed=seed, PFER_thr=PFER_thr, 
                   method=method, pk=pk, implementation=implementation, calib_method=calib_method, blocks=blocks, 
                   param_card1=param_card1, param_card2=param_card2, myfactor=myfactor, 
                   factor_range=factor_range, delta=0))
}


MergeOutput=function(output1, output2, order=TRUE){
  if (output1$stab$K!=output2$stab$K){
    stop("Both objects must have the same number of subsampling iterations (K)")
  }
  if (output1$stab$tau!=output2$stab$tau){
    stop("Both objects must have the same subsampling proportion (tau)")
  }
  if (all(output1$pi_list!=output2$pi_list)){
    stop("Both objects must have the same selection proportion parameters (pi_list)")
  }
  if (all(output1$stab$seed!=output2$stab$seed)){
    stop("Both objects must have been created with the same seed")
  }
  output=output1
  if (length(dim(output$loglik[[1]]))==0){
    output$loglik=rbind(output$loglik, output2$loglik)
    output$PFER=rbind(output$PFER, output2$PFER)
  } else {
    for (block_id in 1:length(output$loglik)){
      output$loglik[[block_id]]=rbind(output$loglik[[block_id]], output2$loglik[[block_id]])
      # print(nrow(output2$loglik[[block_id]]))
      output$PFER[[block_id]]=rbind(output$PFER[[block_id]], output2$PFER[[block_id]])
    }
  }
  output$Q=rbind(output$Q, output2$Q)
  output$S=rbind(output$S, output2$S)
  output$P=rbind(output$P, output2$P)
  output$stab$Lambda=rbind(output$stab$Lambda, output2$stab$Lambda)
  # print(paste("Lambda:", length(output2$stab$Lambda)))
  output$stab$selprop=abind(output$stab$selprop, output2$stab$selprop, along=3)
  
  if (order){
    output=OrderOutput(output)
  }
  
  return(output)
}


OrderOutput=function(output){
  if (!is.null(output$stab$Lambda)){
    if (length(dim(output$loglik[[1]]))==0){
      output$loglik=output$loglik[sort.list(output$stab$Lambda[,1], decreasing=TRUE),]
      output$PFER=output$PFER[sort.list(output$stab$Lambda[,1], decreasing=TRUE),]
    } else {
      for (block_id in 1:length(output$loglik)){
        output$loglik[[block_id]]=output$loglik[[block_id]][sort.list(output$stab$Lambda[,1], decreasing=TRUE),]
        output$PFER[[block_id]]=output$PFER[[block_id]][sort.list(output$stab$Lambda[,1], decreasing=TRUE),]
      }
    }
    output$Q=output$Q[sort.list(output$stab$Lambda[,1], decreasing=TRUE),,drop=FALSE]
    output$S=output$S[sort.list(output$stab$Lambda[,1], decreasing=TRUE),,drop=FALSE]
    output$P=output$P[sort.list(output$stab$Lambda[,1], decreasing=TRUE),,drop=FALSE]
    output$stab$selprop=output$stab$selprop[,,sort.list(output$stab$Lambda[,1], decreasing=TRUE)]
    output$stab$Lambda=output$stab$Lambda[sort.list(output$stab$Lambda[,1], decreasing=TRUE),,drop=FALSE]
  } else {
    if (length(dim(output$loglik[[1]]))==0){
      myQ=output$Q[,1]
      output$loglik=output$loglik[sort.list(myQ, decreasing=FALSE),]
      output$PFER=output$PFER[sort.list(myQ, decreasing=FALSE),]
      output$S=output$S[sort.list(myQ, decreasing=FALSE),,drop=FALSE]
      output$P=output$P[sort.list(myQ, decreasing=FALSE),,drop=FALSE]
      output$stab$selprop=output$stab$selprop[,,sort.list(myQ, decreasing=FALSE)]
      output$Q=output$Q[sort.list(myQ, decreasing=FALSE),,drop=FALSE]
    } else {
      for (block_id in 1:length(output$loglik)){
        myQ=output$Q[,block_id]
        output$loglik[[block_id]]=output$loglik[[block_id]][sort.list(myQ, decreasing=FALSE),]
        output$PFER[[block_id]]=output$PFER[[block_id]][sort.list(myQ, decreasing=FALSE),]
        output$S=output$S[sort.list(myQ, decreasing=FALSE),,drop=FALSE]
        output$P=output$P[sort.list(myQ, decreasing=FALSE),,drop=FALSE]
        output$stab$selprop=output$stab$selprop[,,sort.list(myQ, decreasing=FALSE)]
        output$Q=output$Q[sort.list(myQ, decreasing=FALSE),,drop=FALSE]
      }
    }
  }
  
  return(output)
}


ConstrainCalib=function(out, pk=NULL, PFER_thr=10){
  p=ncol(out$stab$selprop)
  N=p*(p-1)/2
  if (is.null(pk)){
    PFER_pairs=GetPFERPairs(q_list=out$Q, pi_list=out$pi_list, K=out$stab$K, N=N)
    constr_loglik=ifelse(PFER_pairs<PFER_thr, yes=out$loglik, no=NA)
    out$loglik=constr_loglik
  } else {
    nblocks=sum(upper.tri(matrix(NA, ncol=length(pk), nrow=length(pk)), diag=TRUE))
    bigblocks=GetBlockMatrix(pk)
    if (length(PFER_thr)<nblocks){
      PFER_thr_blocks=ceiling(PFER_thr*table(bigblocks[upper.tri(bigblocks)])/N)
    } else {
      PFER_thr_blocks=ceiling(PFER_thr)
    }
    for (block_id in 1:nblocks){
      constr_loglik=ifelse(out$PFER[[block_id]]<PFER_thr_blocks[block_id], yes=out$loglik[[block_id]],no=NA)
      out$loglik[[block_id]]=constr_loglik
    }
  }
  out$PFER_thr=PFER_thr
  return(out)
}


GetHatTheta=function(out){
  hat_theta=arrayInd(which.min(out$loglik), .dim=c(nrow(out$loglik), ncol(out$loglik)))
  hat_theta=c(out$Q[hat_theta[1]], out$P[hat_theta[1]])
  names(hat_theta)=c("q", "pi")
  return(hat_theta)
}


GetHatTheta=function(out, pk=NULL){
  if (is.null(pk)){
    pk=ncol(out$loglik)
  }
  
  if (length(pk)>1){
    bigblocks=GetBlockMatrix(pk)
    hat_Q=hat_pi=NULL
    for (block_id in 1:length(out$loglik)){
      mat=-out$loglik[[block_id]]
      mylambda=out$stab$Lambda[,block_id]
      arrid=which(mat==max(mat,na.rm=TRUE))
      arrid=arrayInd(arrid, .dim=c(nrow(mat), ncol(mat)))
      hat_theta=arrid[which.max(mylambda[arrid[,1]]),]
      hat_Q=c(hat_Q, out$Q[hat_theta[1], block_id])
      hat_pi=c(hat_pi, out$pi_list[hat_theta[2]])
    }
  } else {
    hat_theta=arrayInd(which.min(out$loglik), .dim=c(nrow(out$loglik), ncol(out$loglik)))
    hat_Q=out$Q[hat_theta[1]]
    hat_pi=out$pi_list[hat_theta[2]]
  }
  
  return(list(hat_Q=hat_Q, hat_pi=hat_pi))
}


GetEdgeSelectionProportion=function(out, pk=NULL){
  if (is.null(pk)){
    pk=ncol(out$loglik)
  }
  
  if (length(pk)>1){
    bigblocks=GetBlockMatrix(pk)
    A=matrix(0,nrow=nrow(bigblocks),ncol=ncol(bigblocks))
    for (block_id in 1:length(out$loglik)){
      mat=-out$loglik[[block_id]]
      mylambda=out$stab$Lambda[,block_id]
      arrid=which(mat==max(mat,na.rm=TRUE))
      arrid=arrayInd(arrid, .dim=c(nrow(mat), ncol(mat)))
      hat_theta=arrid[which.max(mylambda[arrid[,1]]),]
      A_block=out$stab$selprop[,,hat_theta[1]]
      A_block[bigblocks!=block_id]=0
      A=A+A_block
    }
  } else {
    hat_theta=arrayInd(which.min(out$loglik), .dim=c(nrow(out$loglik), ncol(out$loglik)))
    A=out$stab$selprop[,,hat_theta[1]]
  }
  
  return(A)
}


GetAdjacency=function(out, pk=NULL, hat_theta=NULL){
  if (!is.null(hat_theta)){
    if (length(pk)>1){
      warning("hat_theta is ignored")
    }
  }
  
  if (is.null(pk)){
    pk=ncol(out$loglik)
  }
  
  if (length(pk)>1){
    bigblocks=GetBlockMatrix(pk)
    A=matrix(0,nrow=nrow(bigblocks),ncol=ncol(bigblocks))
    for (block_id in 1:length(out$loglik)){
      mat=-out$loglik[[block_id]]
      mylambda=out$stab$Lambda[,block_id]
      if (is.null(mylambda)){
        mylambda=out$Q[,block_id]
        mylambda=-mylambda
      }
      arrid=which(mat==max(mat,na.rm=TRUE))
      arrid=arrayInd(arrid, .dim=c(nrow(mat), ncol(mat)))
      hat_theta=arrid[which.max(mylambda[arrid[,1]]),]
      # hat_theta=arrayInd(which.min(out$loglik[[block_id]]), .dim=c(nrow(out$loglik[[block_id]]), ncol(out$loglik[[block_id]])))
      A_block=ifelse(out$stab$selprop[,,hat_theta[1]]>=out$pi_list[hat_theta[2]],1,0)
      A_block[bigblocks!=block_id]=0
      A=A+A_block
    }
  } else {
    if (is.null(hat_theta)){
      hat_theta=arrayInd(which.min(out$loglik), .dim=c(nrow(out$loglik), ncol(out$loglik)))
    }
    A=ifelse(out$stab$selprop[,,hat_theta[1]]>=out$pi_list[hat_theta[2]],1,0)
  }
  
  # print(paste("lambda:",out$stab$Lambda[hat_theta[1]]))
  # print(paste("pi:",out$pi_list[hat_theta[2]]))
  
  print(all(A==t(A))) # check symmetry
  print(sum(A[upper.tri(A)])) # number of selected edges
  rownames(A)=colnames(A)=colnames(out$stab$selprop[,,1])
  
  return(A)
}


# GetGraph=function(Data, adjacency, weighted=NULL, pal=NULL){
#   if (is.null(pal)){
#     pal=c(brewer.pal(12,name="Set3"), brewer.pal(n=8, name="Set2"), brewer.pal(n=4, name="Set1"))
#   }
#   
#   gg=graph_from_adjacency_matrix(adjacency, mode = 'undirected', weighted = weighted)
#   V(gg)$size=as.numeric(as.character(cut(degree(gg), breaks = 3, labels = 1:3)))*3
#   V(gg)$color=pal[as.numeric(Data$Annot$Vertex_color[V(gg)$name])]
#   V(gg)$frame.color=V(gg)$color
#   V(gg)$label.family="sans"
#   E(gg)$color="black"
#   V(gg)$label.cex=as.numeric(as.character(cut(degree(gg), breaks = 3, labels = c(0.4, 0.6, 0.7))))
#   V(gg)$label.color="grey30"
#   E(gg)$width=0.5
#   if (!is.null(weighted)){
#     E(gg)$color=c('red', 'blue', 'forestgreen')[E(gg)$weight]
#   }
#   
#   return(gg)
# }


GetGraph=function(out=NULL, adjacency=NULL, pk=NULL, node_label=NULL, node_color=NULL, 
                  weighted=NULL, satellites=FALSE){
  # either out or adjacency have to be provided
  
  if (is.null(adjacency)){
    if (is.null(out)){
      stop("Either 'out' or 'adjacency' needs to be provided.")
    }
    adjacency=GetAdjacency(out=out, pk=pk)
  }
  
  if (is.null(node_color)){
    node_color=rep("skyblue", ncol(adjacency))
  }
  
  if (is.null(node_label)){
    node_label=colnames(adjacency)
  }
  
  names(node_color)=colnames(adjacency)
  names(node_label)=colnames(adjacency)
  
  mygraph=graph_from_adjacency_matrix(adjacency, mode = 'undirected', weighted = weighted)
  
  if (!satellites){
    mygraph=delete.vertices(mygraph, v=names(degree(mygraph))[degree(mygraph)==0])
  }
  
  V(mygraph)$size=as.numeric(as.character(cut(degree(mygraph), breaks = 3, labels = 1:3)))*3
  V(mygraph)$label=node_label[V(mygraph)$name]
  V(mygraph)$color=node_color[V(mygraph)$name]
  V(mygraph)$frame.color=V(mygraph)$color
  V(mygraph)$label.family="sans"
  E(mygraph)$color="grey40"
  V(mygraph)$label.cex=as.numeric(as.character(cut(degree(mygraph), breaks = 3, labels = c(0.4, 0.6, 0.7))))
  V(mygraph)$label.color="grey30"
  E(mygraph)$width=0.5
  
  if (!is.null(weighted)){
    E(mygraph)$color=c('red', 'blue', 'forestgreen')[E(mygraph)$weight]
  }
  
  return(mygraph)
}


GetPerformance=function(A){
  
  if (is.matrix(A)){
    p=ncol(A)
    N=p*(p-1)/2
    A=A[upper.tri(A)]
  } else {
    N=length(A)
  }
  
  TP=sum(A==3)
  FN=sum(A==2)
  FP=sum(A==1)
  TN=sum(A==0)
  
  sensitivity=TP/(TP+FN)
  specificity=TN/(TN+FP)
  accuracy=(TP+TN)/N
  if (TP+FP>0){
    precision=TP/(TP+FP)
  } else {
    precision=0
  }
  recall=TP/(TP+FN)
  return(list(TP=TP, FN=FN, FP=FP, TN=TN, 
              sensitivity=sensitivity, specificity=specificity, 
              accuracy=accuracy, precision=precision, recall=recall))
}


GetAUC <- function(TPR, FPR){
  # inputs already sorted, best scores first 
  # function from https://www.r-bloggers.com/calculating-auc-the-area-under-a-roc-curve/
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}


GetROC=function(perf, show_AUC=TRUE, col="navy"){
  plot(1-perf$specificity, perf$sensitivity, col=col, pch=19, cex=0.5,
       xlab="FPR", ylab="TPR", las=1, xlim=c(0,1), ylim=c(0,1))
  abline(0,1,lty=2,col="darkgrey")
  if (show_AUC){
    legend("bottomright", bty="n",
           legend=paste("AUC:", round(GetAUC(perf$sensitivity, 1-perf$specificity), digits=3)))
  }
}


GetPrecisionRecall=function(perf, col="navy"){
  plot(perf$recall, perf$precision, col=col, pch=19, cex=0.5,
       xlab="Recall", ylab="Precision", las=1, xlim=c(0,1), ylim=c(0,1))
}


zTransform=function(r){
  return(1/2*log((1+r)/(1-r)))
}


Zpval=function(z, n){
  z=abs(z) # true because a correlation is scaled
  pval=pnorm(z, mean=0, sd=1/sqrt(n-3), lower.tail=FALSE)/2
  if(pval==0|pval==1){
    pval=pnorm(z, mean=0, sd=1/sqrt(n-3), lower.tail=FALSE, log.p=TRUE)-log(2)  
  }
  return(pval)
}


BlockSimulation=function(n,pk,v,u=0.1,prob=0.1,graph="random"){
  d=sum(pk)
  myhuge=huge.generator(n=n, d=d, prob=prob, graph=graph)
  M=GetBlockMatrix(pk)
  if(length(unique(as.vector(M)))!=length(v)){
    stop("The length of v must be the number of blocks.")
  }
  theta=as.matrix(myhuge$theta)
  diag(theta) = 0
  omega=theta
  omega[M==1] = theta[M==1] * v[1]
  omega[M==2] = theta[M==2] * v[2]
  omega[M==3] = theta[M==3] * v[3]
  diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
  sigma = cov2cor(solve(omega))
  omega = solve(sigma)
  x = mvrnorm(n, rep(0, d), sigma)
  sigmahat = cor(x)
  sim = list(data = x, sigma = sigma, sigmahat = sigmahat, 
             omega = omega, theta = Matrix(theta, sparse = TRUE), 
             sparsity = sum(theta)/(d * (d - 1)), graph.type = graph)
  class(sim) = "sim"
  return(sim)
}




