# function tCorr computes the correlation between two tensor product factors at time t
# NOTE: recommended to choose as the index i the one with larger dimension for better performance

tCorr <- function(H, t, i, j, d, normalization = FALSE){

	source("tIncl3.r")

	if((i > length(d)) || (j > length(d))){
	print("tCorr error: position of factor out of bounds!")
	return()
	}

	# compute the evolution matrix
	source("tEvol.r")
	U <- tEvol(H,t)
	Ure <- as(Re(U), "dgCMatrix")
	Uim <- as(Im(U), "dgCMatrix")
	
	idim <- d[i]
	idimtri <- ((idim+1)*idim)/2
	
	jdim <- d[j]
	jdimtri <- jdim^2
	
	# initialize the correlation number
	C <- 0
	
	for(k in 1:idimtri){
		# create the k'th basis matrix
		vk <- c(rep(0,(k-1)),1,rep(0,idimtri-k))
		ek <- matrix(0,idim,idim)
		ek[upper.tri(ek,diag=TRUE)] <- vk
		
		# count off-diagonal elements twice
		if(sum(diag(ek)) == 0){ek <- 2*ek}
		
		# include the n'th basis matrix into the tensor product
		Ek <- tIncl3(ek,i,d)
		
		# real and imaginary parts of the evolved k'th basis matrix
		Ekevre <- Ure %*% Ek %*% t(Ure) + Uim %*% Ek %*% t(Uim)
		Ekevim <- Uim %*% Ek %*% t(Ure) - Ure %*% Ek %*% t(Uim)
		
		for(l in 1:jdimtri){
			# create the l'th basis matrix
			vl <- c(rep(0,(l-1)),1,rep(0,jdimtri-l))
			
			# create basis matrices
			el <- matrix(vl,jdim,jdim)
			
			# include the l'th basis matrices into the tensor product
			El <- tIncl3(el,j,d)
			
			# commutator
			Creup <- Ekevre %*% El - El %*% Ekevre
			Cimup <- Ekevim %*% El - El %*% Ekevim
			
			# add the square root of trace of C squared to the corr. number
			C <- C + sqrt(sum(Creup^2) + sum(Cimup^2))
		}
	}
	
	# normalize, if wanted
	if(normalization){C <- C/(2*idim^2*jdim^2)}
	
	# return correlation number
	return(C)
}
