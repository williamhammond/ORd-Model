#File    : modelAnalysis.R
#Author  : William Hammond
#Date    : 7/30/2015
#Modified: 7/30/2015

conductance_names <- c('GNa', 'gNal', 'gto', 'PCa', 'GKr','GKs','GK1','Gncx','Pnak',
                       'bt')

genAccuracyGraphs <- function(accuracy, path, parameters, normalized = FALSE){
  for (i in 1:10){
    name = conductance_names[i]
    title = paste(name, "Accuracy") 
    jpeg (paste(path,title,'.jpeg', sep =''))
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