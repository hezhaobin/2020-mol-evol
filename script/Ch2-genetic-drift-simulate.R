# title: simulate genetic drift using wright-fisher scheme
# author: Bin He
# date: 2020-09-21

# simulate one trial at a time. this code may be easier to understand
gen_drift <- function(popsize = 100, gen = 100, p0 = 0.5){
    # function to simulate genetic drift using the binomial distribution
	# popsize: population size
	# gen: how many generations to simulate
	# p0: initial frequency
	freq = numeric(gen)
	freq[1] = p0
	for(i in 1:(gen-1)){
		# the key to the function below is to understand how "rbinom()" works. By specifying rbinom(n, size, prob), you are asking the computer to generate n binomially distributed random
		# variable, each one of which is the result of grabbing "size" number of balls from a bag where the probability of getting a black ball is "prob" (if we call black 1 and white 0).
		# so rbinom(1, popsize, freq[i]) would generate one value that corresponds to sampling "popsize" balls from a bag where the frequency of black balls is given by the proportion of
		# black balls in the previous generation (freq[i]). Dividing this number by "popsize" will then give you the frequency of black balls in the currrent generation.
		freq[i+1] = rbinom(1, popsize, freq[i]) / popsize
	}
	plot(1:gen, freq, type = "l", ylim = c(0,1), xlab = "generation", ylab = "allele frequency", main = paste("N =", popsize, sep = " "))
}

# simulate multiple trials ("replay" the evolutionary tape multiple times)
gen_drift1 <- function(popsize = 100, gen = 100, p0 = 0.5, trial = 1){
    # function to simulate genetic drift using the binomial distribution
	# popsize: population size
	# gen: how many generations to simulate
	# p0: initial frequency
	# trial: how many times to simulate the evolution
	freq = matrix(0, ncol = trial, nrow = gen)
	freq[1,] = p0
	for(k in 1:trial){
		for(i in 1:(gen-1)){
			freq[i+1, k] = rbinom(1, popsize, freq[i, k]) / popsize
		}
	}
	matplot(freq, type = c("l"), ylim = c(0,1), col = 1, xlab = "generation", ylab = "allele frequency", main = paste("N =", popsize, sep = " "))
}

