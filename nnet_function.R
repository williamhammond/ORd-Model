# Example of using Neural Networks (nnet) for Prediction Purposes

library(nnet)
library(neuralnet)
library(clusterSim)


#setwd("/Users/mme4362/R_Project_Will/Neural_Net/")

conductance_names <- c('GNa', 'gNal', 'gto', 'PCa', 'GKr','GKs','GK1','Gncx','Pnak',
                       'bt')

results_N <- read.table("results_normal.txt")
results_HF <- read.table("results_HF.txt")
distro_N <- read.table("distro_normal.txt")
distro_HF <- read.table("distro_HF.txt")

features <- do.call(rbind, list(results_N,results_HF))
# features[,1] <- NULL
# features[,51] <- NULL
distroMatrix <- do.call(cbind, list(distro_N,distro_HF))


x <- distroMatrix
x <- t(x)
# x <- x1[sample(1:nrow(x1)),]

y <- features
#y <- t(y)



# Create training and testing data (Note that the data is in random order)
# Usually the taining set is 2/3 and the testing set is 1/3
# train_input <- rbind(x[1:400,],x[501:900,])
# test_input <- rbind(x[401:500,],x[901:1000,])
train_input <- x[501:900,]
test_input <- x[901:1000,]

# train_output <- rbind(y[1:400,],y[501:900,])
# test_output <- rbind(y[401:500,],y[901:1000,])
train_output <- y[501:900,]
test_output <- y[901:1000,] 

training <- cbind(train_input, train_output)
train <- training[sample(nrow(training)),]
#train <- normalize(train)
testing <- cbind(test_input, test_output)
test <- testing[sample(nrow(testing)),]


# Build your Artificial Neural Network Model
model_nnet <- nnet(x=train[,11:26], hidden = 1000000000000000, 
                   y=train[,1:10], size=35,maxit = 10000, linout = TRUE, 
                   trace = TRUE)
x=train[,11:26]
y=train[,1:10]
model_nnet <- neuralnet(AND+OR~, x, hidden=1, threshold=0.01, 
                        stepmax=1000, rep=1, startweights=NULL, algorithm = "rprop+", 
                        learningrate.limit = NULL, 
                        learningrate.factor = list(minus = 0.5, plus = 1.2), 
                        learningrate=NULL, lifesign = "none",lifesign.step = 1000,
                        err.fct = "sse", act.fct = "logistic", 
                        linear.output = TRUE, exclude = NULL, constant.weights = NULL, 
                        likelihood = FALSE)

# Print the summary of the neural net model you built
model_nnet
for(i in 1:10){
  name = conductance_names[i]
  title = paste(name, 'residual')
  jpeg(paste(title,'.jpeg',sep = ''))
  plot (model_nnet$residuals[,i], xlab = 'index', ylab = name, main = title )
  dev.off()
}

# Use the model for prediction pruposes
predicted=predict(model_nnet, test_output, type = "raw")


accuracy=((predicted-test_input)/test_input)*100
# for (i in 1:10){
#   name = condunctance_names[i]
#   title = paste(name, "Accuracy") 
#   jpeg (paste(title,'.jpeg', sep =''))
#   plot(accuracy[,i], ylab = name, main = title)
#   dev.off()
# }

test <- as.data.frame(cbind(test_input[,1],predicted[,1],seq(1,100,1)))

# for (i in 1:10){
#   name = conductance_names[i]
#   title = paste(name, 'comparison normalized')
#   
#   predicted_cond <- as.data.frame(cbind(predicted[,i],seq(1,100,1)))
#   comparison <- as.data.frame(cbind(test_input[,i], seq(1,100,1)))
#   
#   jpeg(paste(title,'.jpeg', sep =''))
#   plot(test_input[,i], type = "p", pch = 22, col = "red", xaxt = "n", yaxt = "n", xlab = "", ylab="")
#   title(main = title, sub = NULL, xlab = "Instances", ylab = name, line = 2, outer = FALSE)
#   axis(2, 1:4, LETTERS[1:4])
#   par(new=TRUE)
#   plot(predicted[,i], type = "p", pch = 19, col = "blue", xlab= "", ylab="")
#   dev.off()
# }
# Print the confusion matrix to see how accurate your model 
# table(test_input, predicted)


