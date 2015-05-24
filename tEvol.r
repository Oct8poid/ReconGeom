# function that takes a Hamiltonian and a time t, and returns the unitary operator giving the time evolution

tEvol <- function(H,t){
	Heig <- eigen(H)
	U <- Heig$vectors %*% diag(exp(1i*t*Heig$values)) %*% t(Conj(Heig$vectors))
	return(U)
}
