plotRegressions<-function(regressionParameters, selectedModels='all', colors=NULL){
  
  colors_models<-if(is.null(colors)){
    c( "#00798c","#d1495b", "#edae49","#66a182","#2e4057", "#8d96a3",'green')}else{
      colors
    }
  if(selectedModels != 'ALL'){
    regressionParameters<-regressionParameters[regressionParameters$Model %in% selectedModels,] 
    }
  
  names(regressionParameters)[c(3,4,6,7)]<-c('I_low', 'I_upper', 'S_low', 'S_upper')
  
  p<-ggplot(data=attr(regressionParameters, "data"), aes(x=T2,y=D47)) + 
    geom_point(shape = 1, alpha=.2) +
    geom_errorbarh(aes(xmin = T2 -Temp_Error ,xmax = T2 + Temp_Error,y = D47,height = 0.01), alpha=.2) +
    geom_errorbar(aes(ymin = D47 - D47_SD,ymax = D47 + D47_SD,x = T2), alpha=.2)
  
  for(i in seq_along(unique(regressionParameters$Model))){
    p<-p+geom_abline(data=regressionParameters[i,], aes(intercept = Intercept, slope =Slope, colour= Model))
    p<-p+geom_abline(data=regressionParameters[i,], aes(intercept = I_low, slope =S_low, colour= Model), linetype='dashed')
    p<-p+geom_abline(data=regressionParameters[i,], aes(intercept =I_upper, slope =S_upper, colour= Model), linetype='dashed')
  }
  
  p+ scale_color_manual(values= colors_models)+theme_bw()
  
}
