# title: recreate part of Figure 2.5 in Chapter 2
# author: Bin He
# date: 2020-09-30

# the goal of this figure is to compute and compare the distribution of allele frequencies in a diploid population of N=10, with an initial frequency of allele A1=0.3

# -----------
# theoretical
# -----------
# set up the system with 10 diploid individuals, with an initial frequency of 0.3 for the A1 allele
N = 10; n = 2*N; p0 = 0.3
# we can easily calculate the probabiliy mass at each possible frequency of A1 alleles (i = 0-20)
pmf1 = dbinom(0:n, n, p0) # dbinom(x, size, prob) calculates the binomial probability for P[X1 = x] = (n choose x) p0^x (1-p0)^(n-x)
# ---- alternative ---
# you can also get the above with:
# i = 1:n; theoretical = choose(n, i)*p0^i*(1-p0)^(n-i)
# for t = 2, the idea is to sum over the possibilities, i.e. Prob(X2 = x) = sum_{k=1 to n}( Prob(X1 = k) * dbinom(x, n, k/n) )
# --- end alternative ---

# function for calculating the distribution at an arbitrary generation
mc.binom <- function(gen, pmf = pmf1){
	# gen: which generation to calculate the distribution for
	# pmf: probablity mass function for the previous generation 
	if(gen == 1){
		return(pmf)
	}
	else{
		i = 0:n # all possible states
		# note that pmf[i+1] stores the frequency for the state i, e.g. pmf[1] stores
		# the frequency of i = 0 A1 alleles
		pmf.new = sapply(i, function(x) sum(pmf[i+1] * dbinom(x, n, i/n)))
		gen = gen - 1
		mc.binom(gen, pmf.new)
		# --- below is the for loop solution ---
		# it works but is slower than the sapply()
		#pmf.new = numeric(n)
		#for(j in 0:n){
			#	pmf.new[j] = sum(pmf[i+1] * dbinom(j, n, i/n))
		#}
		# ------  end for loop solution  -------
	}
}

# ----------
# Simulation
# ----------

# for 1 replicate, calculate the distribution at generation gen
simulate_1 <- function(gen, p = p0){
	# gen: number of generations
	# p: frequency of A1 allele in the previous generation
	if(gen == 1 || p == 0 || p == 1){
		return(rbinom(1, n, p))
	}
	else{
		simulate_1(gen-1, rbinom(1, n, p)/n)
	}
}

simulate_1k <- function(gen, p = p0){
	n.A1 = numeric(1000)
	for( i in 1:1000 ){
		n.A1[i] = simulate_1(gen, p)
	}
	pmf = proportions(tabulate(bin = n.A1+1,nbins =  21))
	# explanation of the command above:
	#  tabulate(n.A1+1, nbins=21) takes the 1000 replicates of simulated n.A1 and count the frequency
	#  of 0, 1, ..., 20. Note that tabulate() assumes the vector contains only integers, and counts the
	#  frequency in bin = 1, 2, ..., nbins. To make it work for 0, 1, 2, ..., 20, I need to add 1. Then
	#  I used proportions() to calculate the normalized frequency for each bin
	return(pmf)
}

# -----------------
# plot and compare
# -----------------

compare <- function(gen){
	# gen: the number of generations to compare the results at
	theory = mc.binom(gen)
	simulate = simulate_1k(gen)
	barplot(rbind(theory, simulate), names.arg = 0:20, beside = TRUE, ylim = c(0, 0.75), 
	legend.text = c("theoretical", "simulation (1000 replicates)"))
}

print("To use this script, type compare(gen), where gen is the number of generations you need to compare at")
