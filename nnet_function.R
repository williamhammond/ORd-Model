library(nnet)
library(neuralnet)
# library(clusterSim)


### Read in Result Data ###
results_N <- read.table("Results/results_normal.txt")
results_HF <- read.table("Results/results_hf.txt")
distro_N <- read.table("Results/distro_normal.txt")
distro_HF <- read.table("Results/distro_HF.txt")

features <- do.call(rbind, list(results_N,results_HF))
distroMatrix <- do.call(cbind, list(distro_N,distro_HF))


x <- distroMatrix
x <- t(x)
# x <- x1[sample(1:nrow(x1)),]

y <- features
#y <- t(y)

### Get # of runs to handle variable sample size ###
combinedSize <- dim(y)[1]
numRuns <- combinedSize / 2 # Divide By 2 b/c y has normal and HF


trainingSize <- (numRuns - numRuns * .2)
testSize <- (numRuns - numRuns * .8)



# Create training and testing data (Note that the data is in random order)
# Usually the taining set is 2/3 and the testing set is 1/3
HF_Start <- numRuns + 1
# train_input <- rbind(x[1:400,],x[501:900,])
# test_input <- rbind(x[401:500,],x[901:1000,])
train_input <- x[HF_Start:(numRuns + trainingSize),]
test_input <- x[(HF_Start + trainingSize):combinedSize,]

# train_output <- rbind(y[1:400,],y[501:900,])
# test_output <- rbind(y[401:500,],y[901:1000,])
train_output <- y[HF_Start:(numRuns + trainingSize),]
test_output <- y[(HF_Start + trainingSize):combinedSize,] 

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
# model_nnet <- neuralnet(y[,1]~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+
#                           V13+V14+V15+V16,x,hidden = 10, err.fct = 'sse',
#                           rep = 2)
# plot(compute(model_nnet, x[,1:16])$net.result)

# Use the model for prediction pruposes
predicted=predict(model_nnet, test_output, type = "raw")

accuracy=((predicted-test_input)/test_input)*100

test <- as.data.frame(cbind(test_input[,1],predicted[,1],seq(1,100,1)))

# Forward Model
target = train[,11:26]
n.col = 16
# forward_predicted = list()
# forward_accuracy = list()
for (i in 1:n.col){
  title <- paste('Forward Model',i,sep='')
  forward_nnet <- nnet(x=train[,1:10], hidden = 1000000000000000, 
                       y=target[,i], size=35,maxit = 10000, linout = TRUE, 
                       trace = TRUE)
  forward_predicted=predict(forward_nnet, test_input, type = "raw")
  forward_accuracy= ((forward_predicted-test_output[,i])/test_output[,i])*100
  
  jpeg (paste(title,'.jpeg', sep =''))
  plot(forward_accuracy)
  dev.off()
}
