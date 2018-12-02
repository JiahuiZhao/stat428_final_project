Interest <- read.csv("/Users/zhouhaowen/Desktop/428pro/Interest.csv", header = TRUE)
Bond <- read.csv("/Users/zhouhaowen/Desktop/428pro/Bond.csv", header = TRUE)

Index <- which(Bond[,1] == Bond[,1])

Bond <- Bond[Index,2]
Interest <- Interest[Index,2]

delta_Bond <- numeric(length(Bond))
delta_Interest <- numeric(length(Interest))

for(i in 2 : length(Bond) ){
  delta_Bond[i-1] <- Bond[i] - Bond[i-1]
  delta_Interest[i-1] <- Interest[i] - Interest[i-1]
  
}

Sigma <- cov(data.frame(delta_Bond,delta_Interest))

S_vector <- matrix(0, nrow = length(delta_Bond), ncol = 2)




delta <- matrix(1, nrow = 2, ncol = 1)
Delta <- diag(1, nrow=2, ncol = 2)


Delta_V <- S_vector %*% delta 
