# function for building local Hamiltonians

buildH <- function(prs, d, ccs){

	source("tIncl3.r")
	
	H <- 0
	
	# Pauli matrices
	Sx <- matrix(c(0,1,1,0),2,2) + 1i*0
	Sy <- matrix(c(0,1i,-1i,0),2,2)
	Sz <- matrix(c(1,0,0,-1),2,2) + 1i*0
	
	S <- list(Sx, Sy, Sz)
	
	for(i in 1:length(prs)){
		for(k in 1:3){
			for(l in 1:3){
				H <- H + ccs[k,l]*tIncl3(S[[k]], prs[[i]][1], d) %*% tIncl3(S[[l]], prs[[i]][2], d)
			}
		}
	}
	return(H)
}
