#File    : modelAnalysis.R
#Author  : William Hammond
#Date    : 7/30/2015
#Modified: 9/08/2015

conductance_names <- c('GNa', 'gNal', 'gto', 'PCa', 'GKr','GKs','GK1','Gncx','Pnak',
                       'bt')
feature_names <- c('v_apd','RMP','v_peak','v_amp','v_notch','V_maxvelocity','T_Triangluation',
                   'V_triangulation','v90','v50','diastolic','I_peak','I_amp','I90','I50','T1_I90')

genAccuracyGraphs <- function(accuracy, path, normalized = FALSE, forward = FALSE){
  for (i in 1:ncol(accuracy)){
    if (!forward){ # We are predicting conductances
      name = conductance_names[i]
      title = paste(name, "Accuracy") 
    }else{ # We are predicting APD features 
      name = feature_names[i]
      title = paste(name, "Forward_Accuracy") 
    }
    jpeg (paste(path,'/',title,'.jpeg', sep =''))
    plot(accuracy[,i], ylab = name, main = title)
    dev.off()
  }
}

genResidualGraphs <- function(path, model, normalized = FALSE ){
  for(i in 1:10){
    name = conductance_names[i]
    title = paste(name, 'residual')
    jpeg(file = paste(path,title,'.jpeg',sep = ''))
    plot (model_nnet$residuals[,i], xlab = 'index', ylab = name, main = title )
    dev.off()
  }
}

genComparisonGraphs <- function (path, test_input, predicted, normalized = FALSE){
  for (i in 1:10){
    name = conductance_names[i]
    title = paste(name, 'comparison normalized')
    
    predicted_cond <- as.data.frame(cbind(predicted[,i],seq(1,100,1)))
    comparison <- as.data.frame(cbind(test_input[,i], seq(1,100,1)))
    
    jpeg(paste(path,title,'.jpeg', sep =''))
    plot(test_input[,i], type = "p", pch = 22, col = "red", xaxt = "n", yaxt = "n", xlab = "", ylab="")
    title(main = title, sub = NULL, xlab = "Instances", ylab = name, line = 2, outer = FALSE)
    axis(2, 1:4, LETTERS[1:4])
    par(new=TRUE)
    plot(predicted[,i], type = "p", pch = 19, col = "blue", xlab= "", ylab="")
    dev.off()
  }
}