# function that implements the inclusion of a matrix M
# into an n-fold tensor product of matrix factors
# whose dimensions are encoded into the vector d

tIncl <- function(M,i,d){

	# packages
	library("Matrix")
	
	# number of tensor product factors
	nFac <- length(d)
	
	# check that i <=  nFac
	if((i > nFac) || (i<1)){
	print("tIncl error: position of factor out of bounds!")
	return()
	}
	
	# check that the matrix M is square
	if(ncol(M) != nrow(M)){
	print("tIncl error: matrix not square!")
	return()
	}
	
	# check that the matrix dimension agrees with
	# the dimension of the i'th factor
	if(ncol(M) != d[i]){
	print("tIncl error: matrix dimension != factor dimension!")
	return()
	}
	
	# compute the dimensions for the surrounding identity matrices
	if(i == 1){ lDim <- 1 } else { lDim <- prod(d[1:i-1], na.rm=TRUE) }
	if(i == nFac){ rDim <- 1 } else { rDim <- prod(d[i+1:nFac], na.rm=TRUE) }
	
	# form the tensor product
	lId <- as(diag(lDim), "dgCMatrix")
	rId <- as(diag(rDim), "dgCMatrix")
	
	if(is.complex(M)){
		MRe <- as(Re(M), "dgCMatrix")
		MIm <- as(Im(M), "dgCMatrix")
		
		PRe <- kronecker(kronecker(lId,MRe),rId)
		PIm <- kronecker(kronecker(lId,MIm),rId)
		
		return(as(PRe, "matrix") + 1i*as(PIm,"matrix"))
	} else {
		M <- as(M, "dgCMatrix")
		P <- kronecker(kronecker(lId,M),rId)
		return(P)
	}
}
