################################################################################
# Fischprojekt
# Funktionen zur Hotspot-Findung
# Datum der ersten Version: 11.04.12
# Letzte Aenderung: 06.06.13
# Autor: Ludwig Bothmann
################################################################################

require(pixmap)
require(cluster)

####################################################################
# Die Funktion hotspot() soll in einer Matrix ein Rechteck um den 
#	Fisch ziehen
#
# Uebergeben werden muss:
#	- A:	Die Matrix mit 0/1 Werten
#
# Ausgegeben wird:
#	- dims.A:	Matrix mit den Grenzen des Rechtecks	
#
####################################################################

hotspot <- function(A){
	
	# Indizes, wo die 1er stehen
	a <- which(A==1, arr.ind=TRUE)
	
	##########################################
	# i) Erst nach der zweiten Spalte ordnen
	##########################################
	
	a <- a[order(a[,2]),]

	# Solange der Indikator auf TRUE steht muss ich weiter machen
	
	ind <- TRUE
	
	while(ind){
		
		# Auf FALSE setzen, wenn nachher in beiden if-Abfragen FALSE rauskommt, bleibt es dabei
		ind <- FALSE
		
		if(diff(a[,2][1:2]) > 1){
					
			# Falls Abstand von erstem zum zweiten > 1 => erster raus
			a <- a[-1,]
			
			ind <- TRUE
		}
		
		# In die zweite Abfrage darf er nur, wenn er in der ersten nicht war.
		if(ind == FALSE){
		
			if(diff(a[,2][(dim(a)[1]-1):dim(a)[1]]) > 1){
		
				# Falls Abstand vom vorletzten zum letzten > 1 => letzter raus		
				a <- a[-dim(a)[1],]
		
				ind <- TRUE
			}
		}
	}



	##########################################
	# ii) Dann nach der ersten Spalte ordnen
	##########################################
	
	a <- a[order(a[,1]),]

	# Solange der Indikator auf TRUE steht muss ich weiter machen
	ind <- TRUE
	while(ind){
		
		# Auf FALSE setzen, wenn nachher in beiden if-Abfragen FALSE rauskommt, bleibt es dabei
		ind <- FALSE
		
		if(diff(a[,1][1:2]) > 1){
					
			# Falls Abstand von erstem zum zweiten > 1 => erster raus
			a <- a[-1,]
			
			ind <- TRUE
		}
		
		# In die zweite Abfrage darf er nur, wenn er in der ersten nicht war.
		if(ind == FALSE){
		
			if(diff(a[,1][(dim(a)[1]-1):dim(a)[1]]) > 1){
		
				# Falls Abstand vom vorletzten zum letzten > 1 => letzter raus		
				a <- a[-dim(a)[1],]
		
				ind <- TRUE
			}
		}
	}
	
	# Gibt es noch ein Problem?
	prob2 <- max(diff(a[,1]))>1 # => Falls TRUE gibt es wohl noch ein Problem / einen einzelnen Pixel innerhalb des Rechtecks um den Fisch
	
	a <- a[order(a[,2]),]
	prob1 <- max(diff(a[,2]))>1 # => Falls TRUE gibt es wohl noch ein Problem / einen einzelnen Pixel innerhalb des Rechtecks um den Fisch
	a <- a[order(a[,1]),]
		
		
	dim1 <- c(min(a[,1]),max(a[,1]))
	dim2 <- c(min(a[,2]),max(a[,2]))
	
	dims.A <- cbind(dim1, dim2)
	
	if(prob1 | prob2) print("Es koennte ein Problem geben.")
	
	
	return(dims.A)
	
	
}


####################################################################
# Die Funktion hotspot2() soll in einer Matrix ein Rechteck um den 
#	Fisch ziehen, spart sich aber das loeschen einzelner Werte
#
# Uebergeben werden muss:
#	- A:	Die Matrix mit 0/1 Werten
#
# Ausgegeben wird:
#	- dims.A:	Matrix mit den Grenzen des Rechtecks	
#
####################################################################

hotspot2 <- function(A){
	
	# Indizes, wo die 1er stehen
	a <- which(A!=0, arr.ind=TRUE)
			
	dim1 <- c(min(a[,1]),max(a[,1]))
	dim2 <- c(min(a[,2]),max(a[,2]))
	
	dims.A <- cbind(dim1, dim2)
	
	return(dims.A)
	
}

####################################################################
# Die Funktion hotspot.mult() soll in einer Matrix ein Rechteck um den 
#	Fisch ziehen. 
#	Variante von hotspot2() fuer mehrere Gruppen
#
# Uebergeben werden muss:
#	- A:	n x 3 - Matrix mit Koordinaten der Punkte + Gruppennamen
#	- i:	Gewuenschte Grupe
#
# Ausgegeben wird:
#	- dims.A:	Matrix mit den Grenzen des Rechtecks	
#
####################################################################

hotspot.mult <- function(A, i){
	
	# Teil der Matrix, der zur Gruppe i gehoert
	A.i <- A[which(A[,3]==i),]
	
	dim1 <- c(min(A.i[,1]),max(A.i[,1]))
	dim2 <- c(min(A.i[,2]),max(A.i[,2]))
	
	dims.A <- cbind(dim1, dim2)
	
	return(dims.A)
	
}


####################################################################
# Die Funktion clean() soll in einer Matrix Punktewolken loeschen, 
#	die kleiner sind als ein bestimmter Wert
#
# Uebergeben werden muss:
#	- A:	Die Matrix mit beliebigen Werten
#	- cut:	Die minimale Groesse einer Wolke, die nicht geloescht werden soll
#
# Ausgegeben wird:
#	- A:	Die bereinigte Matrix
#
####################################################################

#' @import cluster
clean <- function(A, cut){
	
	# falls nur Nullen in der Matrix stehen: direkt wieder ausgeben
	if(sum(A^2)==0) return(A)
	
	if(sum(A^2)!=0){
		y.clust <- which(A!=0, arr.ind=TRUE)
		
		if(dim(y.clust)[1]>1){
		
			########### Anmerkung 17.04.12: Hier koennte ich auch A!=0 machen, 
			#	dann muesste ich nicht die 0/1-Matrix uebergeben, sondern koennte 
			#	direkt die abgeschnittene Y.hat.cut uebergeben...
			#
			#	=> Dann auch in hotspot2()
			#
			#	=> Hab das umgesetzt, mal sehen ob es gescheit funktioniert
			###########
			
			# Mit Single-Linkage und euklidischer Distanz
			clust <- agnes(y.clust, 
								method="single",
								metric="euclidean")
			
			# Wenn der Abstand groesser ist als 1 => neue Gruppe
			groups <- cutree(as.hclust(clust), h=1)
			
			erg <- cbind(y.clust,groups)
		
			# erg nach Gruppen sortieren
			erg <- erg[order(erg[,3]),]
			
			# Wenn Punktewolke kleiner als "cut" Punkte, dann loeschen
			gr.del <- which(table(erg[,3])<cut)	
		
			# Diese Indizes muessen geloescht werden
			del <- erg[which(is.element(erg[,3],gr.del)),1:2]
			
			# Setze Werte auf 0
			
			if(is.vector(del)) 	A[del[1], del[2]] <- 0
			if(!is.vector(del)) A[del] <- 0
		
		}
		
		if(dim(y.clust)[1]==1) A <- 0
		# Noetig, da agnes() mindestens 2 Punkte zum Clustern braucht: Wenn nur einer da ist => loeschen.
		
			
		return(A)
	}
}

####################################################################
# Die Funktion clean.mult() soll in einer Matrix Punktewolken loeschen, 
#	die kleiner sind als ein bestimmter Wert.
#	Variante von clean() fuer multiple Hotspots. Der Output ist anders.
#
# Uebergeben werden muss:
#	- A:				Die Matrix mit beliebigen Werten
#	- cut:			Die minimale Groesse einer Wolke, die nicht geloescht werden soll
#	- maxdist:		Abstand, ab dem 2 Cluster getrennt werden
#	- fastclust:	TRUE: Meine Clusterfunktion, sonst agnes()
#	- speeditup:	TRUE: Noch etwas schneller, evtl. mehr Fehler
#	- t:				Zeitpunkt	
#
# Ausgegeben wird:
#	- A:	Die bereinigte Matrix
#
####################################################################

#' @import cluster
clean.mult <- function(A, 
							  cut, 
							  maxdist,
							  fastclust=FALSE,
							  speeditup=FALSE,
								floodclust=FALSE,
											 floodCpp=TRUE,
							  t,
											 pack=NULL){
	
	
	# falls nur Nullen in der Matrix stehen: direkt wieder ausgeben
# 	if(sum(A^2)==0) return(A)
	
	if(sum(A^2)!=0){
		
		if(fastclust & floodclust){
		
			cat("Achtung: fastclust und floodclust = TRUE => fastclust wird ausgefuehrt \n")		
		}
		
		if(floodCpp & any(c(fastclust, floodclust))){
		
			cat("Achtung: Mehrere Clustermethoden ausgewaehlt => flood.clust.cpp wird ausgefuehrt \n")		
		}
		
		y.clust <- which(A!=0, arr.ind=TRUE)
		
		# Wenn mehr als ein Pixel drin ist
		if(dim(y.clust)[1]>1){
			
			# Dateinamen fuer floodclust festlegen (noetig bei mehreren Kernen)
			if(is.null(pack)){				
				dataName <- "data.bin"
			}else{
				dataName <- paste("data_pack_",pack,".bin", sep="")
			}
			
			# Clustermethode auswaehlen:
			# Praeferiert ist floodCpp!
			if(floodCpp){
				
				if(t%%10 == 0) cat("Clustermethode bei t = ",t,": floodCpp \n")
				# 0/1 - Bild
				mtx <- (A!=0)*1
				groups <- flood.clust.cpp(bild=mtx)
				
			}else{
				if(!fastclust){
					
					if(!floodclust){
						
						if(t%%10 == 0) cat("Clustermethode bei t = ",t,": agnes \n")
						# Clusteranalyse mit Single-Linkage und euklidischer Distanz
						clust <- agnes(y.clust, 
											method="single",
											metric="euclidean")
						
						# Wenn der Abstand groesser ist als maxdist => neue Gruppe
						groups <- cutree(as.hclust(clust), h=maxdist)
						
					}else{
						
						# Braeuchte man hier eine andere Grenze? 
						#		Ja: ca. 390 Pixel, bei ssd-Festplatte sicher niedriger
						# 	Daheim ab 235 Pixeln schneller :-)
						#		Ausserdem waere es sinnvoll, alle drei Methoden
						#		zu vergleichen und dann je nach Pixelanzahl auszuwaehlen, 
						#		wenn floodclust & fastclust == TRUE ist
						if(nrow(y.clust)<235){
							
							if(t%%10 == 0) cat("Clustermethode bei t = ",t,": agnes \n")
							# Clusteranalyse mit Single-Linkage und euklidischer Distanz
							clust <- agnes(y.clust, 
														 method="single",
														 metric="euclidean")
							
							# Wenn der Abstand groesser ist als maxdist => neue Gruppe
							groups <- cutree(as.hclust(clust), h=maxdist)
							
						}else{
							
							if(t%%10 == 0) cat("Clustermethode bei t = ",t,": flood.clust \n")
							# 0/1 - Bild
							mtx <- (A!=0)*1
							groups <- flood.clust(dataName=dataName,
																		mtx=mtx)
						}
					}
					
				}else{
					
					if(nrow(y.clust)<250){
					
						if(t%%10 == 0) cat("Clustermethode bei t = ",t,": agnes \n")
						# Clusteranalyse mit Single-Linkage und euklidischer Distanz
						clust <- agnes(y.clust, 
													 method="single",
													 metric="euclidean")
						
						# Wenn der Abstand groesser ist als maxdist => neue Gruppe
						groups <- cutree(as.hclust(clust), h=maxdist)
						
					}else{
						
						if(t%%10 == 0) cat("Clustermethode bei t = ",t,": fastclust \n")
						# UPDATE 6.6.13:
						# Schneller als agnes() ab ca. 250 Pixeln (mit maxdist=5,
						#		speeditup=TRUE)
						groups <- fast.clust(y.clust, 
																 maxdist=maxdist,
																 speeditup=speeditup,
																 t=t)
	
					}		
				}
			}
			
			erg <- cbind(y.clust,groups)
			
			# erg nach Gruppen sortieren
			erg <- erg[order(erg[,3]),]
			
			# Wenn Punktewolke kleiner als "cut" Punkte, dann loeschen
			gr.del <- which(table(erg[,3])<cut)	
			
			# Diese Indizes muessen geloescht werden
			# Geaendert am 6.6.13, unten zweimal entsprechend
			del <- erg[which(is.element(erg[,3],as.numeric(names(gr.del)))),1:2]

			# ALT
# 			del <- erg[which(is.element(erg[,3],gr.del)),1:2]
			
			# Setze Werte auf 0
			
			if(is.vector(del)){
				# Falls del ein Vektor ist, also nur ein Punkt geloescht wird
				A[del[1], del[2]] <- 0
				clust.erg <- erg[-which(is.element(erg[,3],as.numeric(names(gr.del)))),]
			}else if(is.matrix(del) & dim(del)[1]!=0){
				# Falls del eine Matrix ist, also mehrere Punkte geloescht werden
				A[del] <- 0
				clust.erg <- erg[-which(is.element(erg[,3],as.numeric(names(gr.del)))),]
			}else{
				# Sonst, also wenn nichts geloescht werden muss
				clust.erg <- erg
			}
			
		}
		
		if(dim(y.clust)[1]==1) A <- 0
		# Noetig, da agnes() mindestens 2 Punkte zum Clustern braucht: Wenn nur einer da ist => loeschen.
		
		out <- list(Bild=A,         # Bereinigtes Bild
						Gruppen=clust.erg  # Matrix mit Koordinaten der Punkte
										# und zugehoeriger Gruppe (Cluster, Fisch)
		)
		
		return(out)
	}
}


####################################################################
# Die Funktion dreh() soll ein Bild um xx Grad gegen den Uhrzeigersinn 
#	drehen. Noetig, damit mit pixmap() die Bilder richtig gedreht sind. 
#	Bisher passiert in dieser Funkion real aber gar nichts, da ich 
#	die Ausrichtung der Bilder doch gut finde, so wie sie sind. 
#	Hier koennte ich aber wieder alles auskommentieren.
#
# Uebergeben werden muss:
#	- A:	Die Matrix mit beliebigen Werten
#
# Ausgegeben wird:
#	- A:	Die gedrehte Matrix
#
####################################################################


dreh <- function(A){
	# B <- 0
	# if(!is.null(dim(A))){
		# B <- t(A[,dim(A)[2]:1])
	# }
	# return(B)
	A
}

####################################################################
# Die Funktion pix.plot() soll einen pixmap-plot erstellen von einem Array
#
# Uebergeben werden muss:
#	- A:		Der Array
#	- cellres:	Seitenverhaeltnis der Pixel (1 = quadratisch)
#	- t.grid:	Zeitpunkte der geplotteten Frames
#
# Ausgegeben wird:
#	- nichts, es wird nur geplottet
#
####################################################################

#' @import pixmap
pixplot <- function(A, cellres=1, t.grid){

	if(!is.matrix(A)){
		n.time <- dim(A)[3]
		n <- ceiling(sqrt(n.time))
		
		par(mfrow=c(n,n), mar=c(1.7,1.8,1.7,0.2)+0.2)
			for(t in 1:n.time){
				x <- pixmapGrey(A[,,t], cellres=cellres)
				plot(x, main=paste("t = ",t.grid[t],sep=""))
				box()
			}
	}
	
	if(is.matrix(A)){
		n.time <- 1
		n <- 1
		
		par(mfrow=c(1,1), mar=c(1.7,1.8,1.7,0.2)+0.2)		
			x <- pixmapGrey(A, cellres=cellres)
			plot(x, main=paste("t = ",t.grid,sep=""))
			box()
			
	}
	


}




####################################################################
# Die Funktion fast.clust() soll die Cluster schneller finden als agnes()
#	und wird optional in clean.mult verwendet
#
# Uebergeben werden muss:
#	- y.clust:		Matrix mit Koordinaten der Pixel
#	- maxdist:		Wie weit duerfen Pixel entfernt sein um als Nachbarn
#							angesehen zu werden?
#	- speeditup:	Schnellere Variante
#	- t:				Nummer des Frames
#
# Ausgegeben wird:
#	- groups:		Vektor mit den Gruppenzugehoerigkeiten
#
####################################################################

######################################################################
# 6.6.13: Achtung, bei der Option speeditup passiert noch irgendwas
#	komisches, zumindest bei dem Plot Hotspots_alleframes.jpg
######################################################################

#' @import cluster
fast.clust <- function(y.clust,
							  maxdist=1,
							  speeditup=FALSE,
							  t){

	#################################
	# 1. Initialisierungen
	#################################
	
	# Distanzmatrix aller schwarzen Pixel
	Daisy <- as.matrix(daisy(y.clust))
	
	# T/F: kleiner als maxdist aber nicht 0?
	D <- Daisy<=maxdist & Daisy>0
	
	# T/F: zwischen maxdist und maxdist-1?
	D.k <- Daisy<=maxdist & Daisy>(maxdist-1)
	
	# Anzahl Pixel
	n <- nrow(y.clust)
	
	# Gruppenzugehoerigkeiten
	groups <- rep(0, n)
	
	# Wurden Nachbarn schon bestimmt?
	nb <- rep(0, n) 
	
	# Mit erstem Pixel anfangen
	i <- 1
	
	# Gruppe des 1. Pixels: 1
	i.cluster <- 1
	groups[1] <- i.cluster
	
	# Vektor der Pixel, die noch bearbeitet werden muessen
	j.vec <- NULL # NA?
	
	# Anzahl Iterationen
	n.iter <- 0
	
	# While-Schleife stoppen, wenn do.loop == FALSE
	do.loop <- TRUE
	
	#################################
	# 2. While-Schleife
	#################################

	while(do.loop){
	
		# Aussteigen, wenn mehr als n Iterationen vorbei
		n.iter <- n.iter + 1
		if(n.iter > n) {
			do.loop <- FALSE
			cat("Achtung: Mehr Iterationen als Elemente:",n.iter," \n")
		}
	
		# Nachbarn von i bestimmen
		j <- which(D[i,])
		
		# Nachbarn von i bestimmt
		nb[i] <- 1
		
		# Nachbarn in gleiche Gruppe stecken
		groups[j] <- i.cluster
		
		# Wenn es schnell gehen soll
		if(speeditup){
			
			# Nachbarn, die nah genug dran sind nicht weiter bearbeiten,
			#	aber nur von 0 auf 1 setzen, nicht andersum, sonst schau ich
			#	mir Pixel nochmal an, die ich schon hatte. 
			#	Sprich: Auf 1 setzen, 
			#	wenn es vorher schon so war oder, 
			#	wenn es auf 0 war
			#		UND D.k auch 0 ist, also wenn ich naeher dran bin als der
			#		aeussere Ring
			nb[j] <- nb[j] | (!D.k[i,j] & !nb[j])
		
		}
		
		# Welche Pixel sehe ich mir als naechstes an, die in diesem Cluster
		#	sind?
		# Nur Pixel, die ich noch nicht habe
		j.vec <- intersect(union(j.vec,j),which(nb==0))
		
		# Wenn es noch welche in dem Cluster gibt
		if(length(j.vec)>0){
		
			# Einen Punkt raussuchen und aus j.vec loeschen
			i <- j.vec[1]
			j.vec <- j.vec[-1]
		
		}else{
		
			# Wenn es ueberhaupt noch Punkte gibt, die ich anschauen muss
			if(length(which(nb==0))>0){
			
				# Neuen Gruppennamen uebergeben
				i.cluster <- i.cluster+1
				
				# Neuen Pixel aussuchen
				i <- which(nb==0)[1]
				
				# Hinweis CM 17.06.13: Diese Zeile fehlt, sonst gibt es
				#	Pixel mit Gruppe "0"
				# Diesem Pixel den neuen Gruppennamen geben
				groups[i] <- i.cluster
				
			}else{
				
				# Aufhoeren
				do.loop <- FALSE
				if(t%%10 == 0) cat("Clustern Frame",t,": Insgesamt",n.iter,"Iterationen \n")
# 				cat("maxdist=",maxdist, "\n")
# 				cat("fastclust=",fastclust, "\n")
# 				cat("speeditup=",speeditup, "\n")
			} # end else
		} #end else
		
	} # end while
	
	if(0 %in% groups) {
		
		# Pixel in Gruppe 0
		n.0 <- sum(groups %in% 0)
		
		cat("Achtung:",n.0,"Pixel mit Gruppe 0 \n")
	}
	
	# Zurueckgeben: Gruppenzugehoerigkeiten
	return(groups)
	
}

########################################################
# Funktionen fuer floodclust
########################################################

loadMatrix <- function(dataName = "data") {
	con <- file(dataName, 'rb')
	on.exit(close(con))
	
	d <- readBin(con, 'integer', 2)
	type <- readBin(con, 'integer', 0)
	t(structure(readBin(con, "double", prod(d)), dim=d))
}

saveMatrix <- function(dataName, mtx) {
	con <- file(dataName, 'wb')
	on.exit(close(con))
	
	writeBin(dim(mtx), con)
	writeBin(c(t(mtx)), con) 
	flush(con)
}


flood.clust <- function(dataName="data.bin", mtx){
	saveMatrix(dataName,mtx)
	system.call <- paste("./FloodFill", dataName,"50000 1", sep=" ")
	system(system.call)
	mat1 <- loadMatrix(dataName)
	
	groups <- mat1[which(mat1!=0)]-1
	return(groups)
}

# Funktion, die direkt C aufruft und nicht ueber die Festplatte muss

flood.clust.cpp <- function(bild){
	
	# Bild transponieren, da in C erst nach Zeilen durchgegangen wird
	tbild <- t(bild)
	
	BildVec <- .C("floodfunc",
								Bild=as.double(as.vector(tbild)),
								as.integer(nrow(bild)),
								as.integer(ncol(bild)),
								PACKAGE="sonar")$Bild
	
	Bild <- matrix(BildVec, nrow=nrow(bild), ncol=ncol(bild), byrow=TRUE)
	
	groups <- Bild[which(Bild!=0)]-1
	return(groups)
	
}

########################################################
# Funktionen flood.clust Ende
########################################################
