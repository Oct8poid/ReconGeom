# function for building local Hamiltonians

# prs is a list of index pairs to be coupled by the local hamiltonian
# d is a vector of matrix factor dimensions, currently restricted to the value 2
# css is a 3x3 matrix of factors multiplying terms in the local hamiltonian

buildH <- function(prs, d, ccs){

	source("tIncl.r")
	
	H <- 0
	
	# Pauli matrices
	Sx <- matrix(c(0,1,1,0),2,2) + 1i*0
	Sy <- matrix(c(0,1i,-1i,0),2,2)
	Sz <- matrix(c(1,0,0,-1),2,2) + 1i*0
	
	S <- list(Sx, Sy, Sz)
	
	for(i in 1:length(prs)){
		for(k in 1:3){
			for(l in 1:3){
				H <- H + ccs[k,l]*tIncl(S[[k]], prs[[i]][1], d) %*% tIncl(S[[l]], prs[[i]][2], d)
			}
		}
	}
	return(H)
}
