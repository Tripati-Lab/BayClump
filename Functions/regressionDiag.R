library("superdiag")
library("ggmcmc")

##Notes: Create function to retrieve model diagnostics


regressionDiagnostics<- function(CompleteModelFit,
         interceptPrior=c(0.068, 0.323), distIntercept='Uniform',
         slopePrior=c(0.033, 0.046), distSlope='Uniform', export_ggmcmc=F, 
         prefix='Diagnostics_', Material=T){
  
  if(Material){
    
    smPost<- lapply(1:3 , function(q){
      BM<-CompleteModelFit[[q+4]]
      n=BM$n.iter
      interceptPrior=if(distIntercept =='Uniform' ){runif(n,interceptPrior[1], interceptPrior[2])}else{
        rnorm(n,interceptPrior[2],interceptPrior[2])
      }
      slopePrior=if(distSlope =='Uniform' ){runif(n, slopePrior[1],slopePrior[2])}else{
        rnorm(n,slopePrior[1],slopePrior[2])
      }
      
      namcols<-colnames(BM$BUGSoutput$sims.matrix)
      basic_no<-length( grep('[',namcols , fixed = T))
      
      priors<-if(basic_no ==0 ){
        rbind(cbind.data.frame(AllID='alpha',AllParams=interceptPrior, group='Prior'),
              cbind.data.frame(AllID='beta',AllParams=slopePrior, group='Prior'))
        
      }else{
        
        priors<- do.call(rbind.data.frame,lapply(seq_along(head(namcols,-2)), function(x){
          if(length(grep('beta',namcols[x]))==0){
            cbind.data.frame(AllID=head(namcols,-2)[x],AllParams=interceptPrior, group='Prior')
          }else{
            cbind.data.frame(AllID=head(namcols,-2)[x],AllParams=slopePrior, group='Prior')
          }
        }))
        
      }
      
      BM_MCMC<-as.mcmc(BM)
      if( export_ggmcmc ==T ){
        ggmcmc(ggs(BM_MCMC),file = paste0(prefix,'Model_',q, '.pdf'))
      }
      
      MDD<-capture.output(superdiag(BM_MCMC, burnin = 0))
      
      outMCMC<-list(chains=MyBUGSChains(BM$BUGSoutput, head(namcols,-2)),
                    priorpost=MyBUGSHist(BM$BUGSoutput, head(namcols,-2), priors=priors),
                    NumDiagnostics=MDD)
      
      
      
      return(outMCMC)
    })
    
    names(smPost)<-names(CompleteModelFit)[-c(1:4)]
    smPost
    
  }else{

      BM<-CompleteModelFit[[3]]
      n=BM$n.iter
      interceptPrior=if(distIntercept =='Uniform' ){runif(n,interceptPrior[1], interceptPrior[2])}else{
        rnorm(n,interceptPrior[2],interceptPrior[2])
      }
      slopePrior=if(distSlope =='Uniform' ){runif(n, slopePrior[1],slopePrior[2])}else{
        rnorm(n,slopePrior[1],slopePrior[2])
      }
      
      namcols<-colnames(BM$BUGSoutput$sims.matrix)
      basic_no<-length( grep('[',namcols , fixed = T))
      
      priors<-if(basic_no ==0 ){
        rbind(cbind.data.frame(AllID='alpha',AllParams=interceptPrior, group='Prior'),
              cbind.data.frame(AllID='beta',AllParams=slopePrior, group='Prior'))
        
      }else{
        
        priors<- do.call(rbind.data.frame,lapply(seq_along(head(namcols,-2)), function(x){
          if(length(grep('beta',namcols[x]))==0){
            cbind.data.frame(AllID=head(namcols,-2)[x],AllParams=interceptPrior, group='Prior')
          }else{
            cbind.data.frame(AllID=head(namcols,-2)[x],AllParams=slopePrior, group='Prior')
          }
        }))
        
      }
      
      BM_MCMC<-as.mcmc(BM)
      if( export_ggmcmc ==T ){
        ggmcmc(ggs(BM_MCMC),file = paste0(prefix,'BL_Model.pdf'))
      }
      
      MDD<-capture.output(superdiag(BM_MCMC, burnin = 0))
      
      outMCMC<-list(chains=MyBUGSChains(BM$BUGSoutput, head(namcols,-2)),
                    priorpost=MyBUGSHist(BM$BUGSoutput, head(namcols,-2), priors=priors),
                    NumDiagnostics=MDD)
      
      
      outMCMC
    
  }
  
 
}

