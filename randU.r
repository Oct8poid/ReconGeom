# function that returns a random unitary matrix of the given dimension
# NOTE: the distribution is not uniform w/r to Haar measure

randU <- function(d,radMin,radMax){
	HRe <- matrix(runif(d^2,radMin,radMax),d,d)
	HIm <- matrix(runif(d^2,radMin,radMax),d,d)
	
	H <- HRe + t(HRe) + 1i*(HIm - t(HIm))
	
	Heig <- eigen(H)
	U <- Heig$vectors %*% diag(exp(1i*Heig$values)) %*% t(Conj(Heig$vectors))
	return(U)
}
