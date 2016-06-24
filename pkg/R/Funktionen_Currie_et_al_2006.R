################################################################################
# Fischprojekt
# Funktionen zur 3D-Glättung
# Datum der ersten Version: 28.03.12
# Letzte Änderung: 22.10.12
# Autor: Ludwig Bothmann
################################################################################

# # Pfad des Ordners festlegen
# # Uni:
# if(uni == TRUE) pfad.basis <- "/home/bothmannlu/Documents/"
# 
# # Daheim:
# if(uni == FALSE) pfad.basis <- "/Users/Ludwig/Documents/Uni/"

################################################################################
# Dieser Code ist wie folgt aufgebaut...
#
#	1.) Funktionen aus dem Paper von Currie et al. (2006)
#	2.) Eigene Funktionen zur 3D-Glättung
#
################################################################################

# library(Rniftilib)

require(splines)

# # Funktion plot.nifti2() von Stefanie Kalus laden, um nifti-Objekte anzusehen
# source("../../Diplomarbeit/R-Code/Stefanie/LookAtNifti_Funktionen_Stefanie.R")
# 
# # Farbig: plot.nifti4()
# source("../../Diplomarbeit/R-Code/Stefanie/LookAtNifti_Funktionen_Stefanie_Farbig.R")


##########################################################
# 1.) Funktionen aus dem Paper von Currie et al. (2006)
##########################################################


# Row tensor of a matrix X
Rten <- function(X){
	one <- matrix(1,1,ncol(X))
	kronecker(X,one)*kronecker(one,X)	
}

# H-transform of an array A by a matrix X
H <- function(X,A){
	
	d <- dim(A)
	M <- matrix(A, nrow=d[1])
	XM <- X %*% M
	array(XM, c(nrow(XM), d[-1]))
	
}

# Rotation of an array A
Rotate <- function(A){
	
	d <- 1:length(dim(A))
	d1 <- c(d[-1], d[1])
	aperm(A, d1)
}

# Rotated H-transform of an array A by a matrix X
#	=> entspricht rho?
RH <- function(X, A) {
	
	Rotate(H(X,A))
}

##########################################################
# 2.) Eigene Funktionen zur 3D-Glättung
##########################################################

#######################################################################
# Die Funktion kq.3D() soll die ML-Schätzer für eine 3-dimensionale 
# Glättung ausgeben. 
#
# Übergeben werden muss:
# 	- Y:		Datenarray der abhängigen Variable
# 	- df:		3D-Vektor, der angibt, wie viele Freiheitsgrade die marginalen 
#					B-spline-Basen haben sollen
# 	- degree:	3D-Vektor, der den Grad der Splines angibt
# 
#	
# Ausgegeben wird:
# 	- Theta: 	Array mit den ML-Schätzern
# 	- Y.hat:	Array mit dem geschätzten Effekt für jeden Voxel
# 	- Y.res:	Array mit den Residuen
#
########################################################################



#' @import splines
kq.3D <- function(Y, df, degree,
									res=FALSE){
	
	# Marginale Designmatrizen bestimmen
	X1 <- bs(1:dim(Y)[1], df=df[1], degree=degree[1], intercept=TRUE)
	X2 <- bs(1:dim(Y)[2], df=df[2], degree=degree[2], intercept=TRUE)
	X3 <- bs(1:dim(Y)[3], df=df[3], degree=degree[3], intercept=TRUE)
	
	# ML-Schätzer ausrechnen
	X1.tilde <- solve( t(X1) %*% X1 ) %*% t(X1)
	X2.tilde <- solve( t(X2) %*% X2 ) %*% t(X2)
	X3.tilde <- solve( t(X3) %*% X3 ) %*% t(X3)
	
	Theta <- RH( X3.tilde, RH(X2.tilde, RH(X1.tilde, Y)))

	# Geschätzter Effekt
	Y.hat <- RH(X3, RH(X2, RH(X1, Theta)))
	
	if(res){
		
		# Residuen
		Y.res <- Y - Y.hat
	
		back <- list(Theta = Theta, Y.hat = Y.hat, Y.res = Y.res)
		
	}else{
		
		back <- list(Theta = Theta, Y.hat = Y.hat)
		
	}
	
	return(back)
}

#######################################################################
# Die Funktion array2nifti() soll aus einem Array ein Nifti-File machen
#
# Übergeben werden muss:
# 	- A:		Array
#	
# Ausgegeben wird:
# 	- A.nifti:	Nifti-File 
#
########################################################################

# #' @import Rniftilib
# array2nifti <- function(A){
# 	
# 	A.nifti <- nifti.image.new()
# 	A.nifti$dim <- dim(A)
# 	A.nifti[] <- A*1
# 	
# 	return(A.nifti)
# 	
# }

#######################################################################
# Die Funktion plot.diff() soll die summierten quadrierten Differenzen 
#	in X3-Richtung plotten und ausgeben
#
# Übergeben werden muss:
# 	- A:		Array
#
# Optional kann übergeben werden:
#	- quad:		TRUE => Quadrierte Differenzen, sonst einfache
#	- cols:		Satz an Farben
#	- output:	TRUE => Ergebnise wird ausgegeben
#	
# Ausgegeben wird:
# 	- A.diff:	Matrix mit dem Ergebnis
#
########################################################################

plot.diff <- function(A, quad=TRUE, cols=rainbow(12), output=FALSE){
	
	# A.diff <- matrix(nrow=dim(A)[1], ncol=dim(A)[2])

	if(quad){	
	
		# Sys.time()
		# for(i in 1:dim(A)[1]){
			# for(j in 1:dim(A)[2]){
				# A.diff[i,j] <- sum(diff(A[i,j,])^2)
			# }
		# }
		# Sys.time()
		
		# Sys.time()
		# A.diff1 <- apply(A, 1:2, sum.diff2)
		# Sys.time()
		
		# sum((A.diff-A.diff1)^2)
		# => Geht gleich schnell und es kommt das gleiche raus
		# => Ich nehme das, was weniger Code ist
		
		A.diff <- apply(A, 1:2, sum.diff2)
	
	}

	if(!quad){	
		
		# Sys.time()
		# for(i in 1:dim(A)[1]){
			# for(j in 1:dim(A)[2]){
				# A.diff[i,j] <- sum(diff(A[i,j,]))
			# }
		# }
		
		# Sys.time()
		
		# Sys.time()
		# A.diff1 <- apply(A, 1:2, sum.diff)
		# Sys.time()
		
		# sum((A.diff-A.diff1)^2)
		# => Geht gleich schnell und es kommt das gleiche raus
		# => Ich nehme das, was weniger Code ist
		
		A.diff <- apply(A, 1:2, sum.diff)		
		
	}

	
	image(1:dim(A)[1], 1:dim(A)[2], A.diff, col=cols)
	
	if(output)	return(A.diff)
	
}



###############################################
# sum.diff() und sum.diff2() sind Hilfsfunktionen, 
#	die für plot.diff benötigt werden
###############################################

sum.diff <- function(a){
	
		# Summe der Differenzen
		sum(diff(a))
	}

sum.diff2 <- function(a){
		
		# Summe der quadrierten Differenzen
		sum(diff(a)^2)	
	}


