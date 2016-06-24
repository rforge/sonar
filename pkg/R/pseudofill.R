pseudoflood <- function(bild, alteFarbe, neueFarbe){

	.C("pseudoflood",
		Bild=as.double(as.vector(bild)),
		as.integer(nrow(bild)),
		as.integer(ncol(bild)),
		as.double(alteFarbe),
		as.double(neueFarbe),
		PACKAGE="sonar")$Bild
		
}

# flood.clust <- function(dataName="data.bin", mtx){
# 	saveMatrix(dataName,mtx)
# 	system.call <- paste("./FloodFill", dataName,"50000 1", sep=" ")
# 	system(system.call)
# 	mat1 <- loadMatrix(dataName)
# 	
# 	groups <- mat1[which(mat1!=0)]-1
# 	return(groups)
# }

# flood.clust.cpp <- function(bild){
# 
# 	# Bild transponieren, da in C erst nach Zeilen durchgegangen wird
# 	tbild <- t(bild)
# 	
# 	BildVec <- .C("floodfunc",
# 		Bild=as.double(as.vector(tbild)),
# 		as.integer(nrow(bild)),
# 		as.integer(ncol(bild)),
# 		PACKAGE="sonar")$Bild
# 	
# 	Bild <- matrix(BildVec, nrow=nrow(bild), ncol=ncol(bild), byrow=TRUE)
# 	
# 	groups <- Bild[which(Bild!=0)]-1
# 	return(groups)
# 
# }


# conv <- function(a, b){
# 	.C("convolve",
# 		as.double(a),
# 		as.integer(length(a)),
# 		as.double(b),
# 		as.integer(length(b)),
# 		ab = double(length(a) + length(b) - 1),
# 		PACKAGE="sonaR")$ab
# 	
# }