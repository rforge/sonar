################################################################################
# Fischprojekt
# Funktionen zur Glaettung und Hotspotfindung bei ddf-files
# Datum der ersten Version: 10.05.12
# Letzte Aenderung: 22.10.12
# Autor: Ludwig Bothmann
################################################################################

#############################
# Pfad des Ordners festlegen
#############################

# # Uni:
# if(uni == TRUE) pfad.basis <- "/home/bothmannlu/Documents/"
# 
# # Daheim:
# if(uni == FALSE) pfad.basis <- "/Users/Ludwig/Documents/Uni/"

#############################
# Funktionen laden
#############################

# # Funktionen ddf
# source("Funktionen_read_plot_ddf.R")
# 
# # Funktionen Hotspots
# source("Funktionen_Hotspots.R")
# 
# # Funktionen aus Currie et al. 2006
# source("Funktionen_Currie_et_al_2006.R")



####################################################################
# Die Funktion find.hotspots() soll ein lineares Array Modell 
#	schaetzen und dann die Hotspots finden.
#
# Uebergeben werden muss:
#	- Y:			Der Datensatz als 3D-Array wobei die 3. Dimension die Zeit ist
# 	- df:			3D-Vektor, der angibt, wie viele Freiheitsgrade die marginalen 
#					B-spline-Basen haben sollen
# 	- degree:		3D-Vektor, der den Grad der Splines angibt
#	- signal.neg: 	TRUE, falls das Signal des Fisches negativ ist. Default ist FALSE
#	- a.1:			Grenzwert, bei dem das gemittelte Signal abgeschnitten wird
#	- a.2:			Grenzwert, bei dem das absolute Signal abgeschnitten wird
#	- cut :			gibt die minimale Groesse an, die eine Punktewolke haben muss 
#						um nicht geloescht zu werden
#	- c:			Anzahl Pixel, die in jede Richtung zusaetzlich in das Hotspot-Rechteck
#						miteinbezogen werden
#
# Ausgegeben wird eine Liste die enthaelt:
#	- Teil_01:			Liste mit den binaerisierten Hotspots - nur Ausschnitt
#	- Ganz_01:			Liste mit den binaerisierten Hotspots
#	- Teil_Y.hat:		Liste mit den Hotspots - Y.hat - nur Ausschnitt
#	- Ganz_Y.hat:		Liste mit den Hotspots - Y.hat
#	- Teil_Y.roh:		Liste mit den Hotspots - Y.roh - nur Ausschnitt
#	- Dims:				Liste mit den Koordinaten der Hotspots
#	- Y.sum:				Summe ueber Ganz_Y.hat
#	- Y.hat:				Y.hat
#	- Y:					Die uebergebenen Daten
#	- df:					siehe oben
#	- Y.hat.zentr:		Y.hat.zentr
#	- Y.hat.cut:		Y.hat.cut
#
####################################################################


find.hotspots <- function(Y, df, degree=c(3,3,3), signal.neg=FALSE, 
								  a.1, a.2, 
								  cut, c=5){
	
	# Schaetzung
	kq1 <- kq.3D(Y, df = df, degree)
	
	Theta <- kq1[[1]]
	Y.hat <- kq1[[2]]
	# Y.res <- kq1[[3]]
	
	###################################################
	# MW ueber die Zeit fuer jeden Pixel abziehen
	###################################################
	
	
	Y.m <- apply(Y.hat, 1:2, mean)
	Y.M <- array(Y.m, dim=dim(Y))
	Y.hat.zentr  <- Y.hat - Y.M
	
	###############################################################################
	# Werte von Y.hat auf 0 setzen an Stellen, wo das 
	#	transformierte Y.hat gross ist.
	###############################################################################
	
	Y.hat.cut <- Y.hat
	
	if(signal.neg == TRUE)	Y.hat.cut[which(Y.hat.zentr > a.1)] <- 0	
	if(signal.neg == FALSE)	Y.hat.cut[which(Y.hat.zentr < a.1)] <- 0	
	
	##########################################
	# Werte von Y.hat auf 0 setzen an Stellen, wo das 
	#	Y.hat gross ist.
	##########################################
	
	if(signal.neg == TRUE) Y.hat.cut[which(Y.hat > a.2)] <- 0
	if(signal.neg == FALSE) Y.hat.cut[which(Y.hat < a.2)] <- 0

	############################################################
	# Rechtecke bestimmen, auch wenn noch Dreck im Bild ist
	############################################################
	
	
	# Liste initialisieren in die jeweils nur die Daten fuer die Hotpots 
	#	geschrieben werden fuer
	
	# ... 0/1 - Bilder
	Liste.y01 <- list()

	# ... 0/1 - Bilder, aber das ganze Bild
	Liste.y01.all <- list()
		
	# ... Y.hat
	Liste.y.hat <- list()
	
	# ... Y.hat, aber das ganze Bild
	Liste.y.hat.all <- list()
	
	# ... Y.roh
	Liste.y.roh <- list()
	
	# Dimensionen der Hotspots
	Liste.dims <- list()
	
	# Summenbild initialisieren
	Y.sum <- matrix(0,nrow=dim(Y)[1], ncol=dim(Y)[2])
	
	# In dieser for-Schleife werden fuer jeden Zeitpunkt die Hotspots berechnet
	
	
	for(t in 1:dim(Y)[3]){
		
#		print(t)
		
		# Y.hat.cut - Bild zu diesem Zeitpunkt
		A <- Y.hat.cut[,,t]
	
		# Bild reinigen (mit Clusteranalyse)
		B <- clean(A,cut)
	
		# 0/1 Bild des bereinigten Bildes
		C <- B
		C[which(C!=0)] <- 1
			
		# Summenbild aller bereinigten Bilder
		Y.sum <- Y.sum + B
		
		# Dimensionen des Rechtecks, falls mindestens ein Pixel schwarz
	
		# Falls im bereinigten Bild noch schwarze Pixel sind: Rechteck zeichnen
		if(sum(B)>0){
			
			# Dimensionen der Hotpots berechnen
			dims <- hotspot2(B)
		
			# => Ausschnitt aus den Bildern
				
			# 0/1 Bild
			Liste.y01[[t]] <- C[(dims[1,1]-c):(dims[2,1]+c),(dims[1,2]-c):(dims[2,2]+c)]
		
			# 0/1 Bild
			Liste.y01.all[[t]] <- C
			
			# Y.hat
			Liste.y.hat[[t]] <- Y.hat[(dims[1,1]-c):(dims[2,1]+c),(dims[1,2]-c):(dims[2,2]+c),t]
			
			# Y.hat - ganzes Bild
			Liste.y.hat.all[[t]] <- B
			
			# Y.roh
			Liste.y.roh[[t]] <- Y[(dims[1,1]-c):(dims[2,1]+c),(dims[1,2]-c):(dims[2,2]+c),t]
	
			# Dimensionen des Rechtecks
			Liste.dims[[t]] <- dims + matrix(c(-c,-c,+c,+c), byrow=TRUE, nrow=2)
	
		}
		
		
		# Falls im bereinigten Bild keine schwarzen Pixel mehr sind: Alles auf 0 setzen
		if(sum(B)==0){
			Liste.y01[[t]] <- 0
			Liste.y01.all[[t]] <- 0
			Liste.dims[[t]] <- 0
			Liste.y.hat[[t]] <- 0
			Liste.y.hat.all[[t]] <- 0
			Liste.y.roh[[t]] <- 0
		}
		
	}

	out <- list(Teil_01=Liste.y01, Ganz_01=Liste.y01.all, 
					Teil_Y.hat=Liste.y.hat, Ganz_Y.hat=Liste.y.hat.all, 
					Teil_Y.roh=Liste.y.roh, Dims=Liste.dims, Y.sum=Y.sum, 
					Y.hat=Y.hat, Y=Y, df=df,
					Y.hat.zentr=Y.hat.zentr,
					Y.hat.cut=Y.hat.cut)
	
	return(out)
}

####################################################################
# Die Funktion find.mult.hotspots() soll ein lineares Array Modell 
#	schaetzen und dann die Hotspots finden.
#	Variante von find.hotspots(), die mehrere Fische pro Bild zulaesst.
#
# Uebergeben werden muss:
#	- Y:			Der Datensatz als 3D-Array wobei die 3. Dimension die Zeit ist
# 	- df:			3D-Vektor, der angibt, wie viele Freiheitsgrade die marginalen 
#					B-spline-Basen haben sollen
# 	- degree:		3D-Vektor, der den Grad der Splines angibt
#	- signal.neg: 	TRUE, falls das Signal des Fisches negativ ist. Default ist FALSE
#	- a.1:			Grenzwert, bei dem das gemittelte Signal abgeschnitten wird
#	- a.2:			Grenzwert, bei dem das absolute Signal abgeschnitten wird
#	- cut :			gibt die minimale Groesse an, die eine Punktewolke haben muss 
#						um nicht geloescht zu werden
#	- c:			Anzahl Pixel, die in jede Richtung zusaetzlich in das Hotspot-Rechteck
#						miteinbezogen werden
#	- maxdist:	Abstand, ab dem 2 Cluster getrennt werden
#
# Ausgegeben wird eine Liste die enthaelt:
#	- Y.01:				Array mit den binaerisierten Hotspots
#	- Y.sum:				Summe ueber Ganz_Y.hat
#	- Y.hat:				Y.hat
#	- Y:					Die uebergebenen Daten
#	- df:					siehe oben
#	- Y.hat.zentr:		Y.hat.zentr
#	- Y.hat.cut:		Y.hat.cut
#	- n.groups:			Anzahl der Hotspots pro Zeitpunkt
#	- hots:				Liste fuer Hotspots, enthaelt zu jedem Zeitpunkt
#								und jeder Gruppe: y01, y01.all,
#			 					y.hat, y.roh, dims,
#								dims.c,
# 								track = NA,species
#
####################################################################


find.mult.hotspots <- function(Y, 
															 df, 
															 degree=c(3,3,3), 
															 signal.neg=FALSE,
															 a.1, 
															 a.2, 
															 cut, 
															 c=5,
															 maxdist=1,
															 fastclust=FALSE,
															 speeditup=FALSE,
															 floodclust=FALSE,
															 floodCpp=TRUE,
															 do.clust=TRUE,
															 pack=NULL,
															 plot.hots.mult=TRUE,
															 save.jpg=TRUE){
	
	cat(as.character(Sys.time())," Beginn Glaettung \n") 
	
	# Schaetzung
	kq1 <- kq.3D(Y, df = df, degree)
	
	# Theta <- kq1[[1]]
	Y.hat <- kq1[[2]]
	
	cat(as.character(Sys.time())," Ende Glaettung \n") 
	
	cat(as.character(Sys.time())," Beginn MW-Bereinigung \n") 
	
	###################################################
	# MW ueber die Zeit fuer jeden Pixel abziehen
	###################################################
	
	Y.m <- apply(Y.hat, 1:2, mean)
	Y.M <- array(Y.m, dim=dim(Y))
	Y.hat.zentr  <- Y.hat - Y.M
	
	###############################################################################
	# Werte von Y.hat auf 0 setzen an Stellen, wo das 
	#	transformierte Y.hat gross ist.
	###############################################################################
	
	Y.hat.cut <- Y.hat
	
	if(signal.neg == TRUE)	Y.hat.cut[which(Y.hat.zentr > a.1)] <- 0	
	if(signal.neg == FALSE)	Y.hat.cut[which(Y.hat.zentr < a.1)] <- 0	
	
	##########################################
	# Werte von Y.hat auf 0 setzen an Stellen, wo das 
	#	Y.hat gross ist.
	##########################################
	
	if(signal.neg == TRUE) Y.hat.cut[which(Y.hat > a.2)] <- 0
	if(signal.neg == FALSE) Y.hat.cut[which(Y.hat < a.2)] <- 0
	
	cat(as.character(Sys.time())," Ende MW-Bereinigung \n") 
	
	############################################################
	# Rechtecke bestimmen, auch wenn noch Dreck im Bild ist
	############################################################
	
	cat(as.character(Sys.time())," Beginn Initialisierung Listen \n") 
	
	# Liste initialisieren in die jeweils nur die Daten fuer die Hotpots 
	#	geschrieben werden fuer
	
	# ... 0/1 - Bilder
	# Liste.y01 <- list()
	
	# ... 0/1 - Bilder, aber das ganze Bild
	# Liste.y01.all <- list()
	# => besser als Array?
	
	Y.01 <- array(dim=dim(Y))
		
	# ... Y.hat
	# Liste.y.hat <- list()
	
	# ... Y.hat, aber das ganze Bild
	# Liste.y.hat.all <- list()
	# => besser als Array?
	# Ganz_Y.hat
	
	# ... Y.roh
	# Liste.y.roh <- list()
	# => besser als Array?
	
	# Ergebnisse fuer jeden Zeitpunkt
	hots <- list()
	
	# Liste in die nur die 0/1 - Hotspots gespeichert werden
	hots.y01 <- list()
	
	# Liste in die nur die Dimensionen der Hotspots gespeichert werden
	dims.l <- list()
	
	# Dimensionen der Hotspots
	#	Liste.dims <- list()
	
	# Summenbild initialisieren
	Y.sum <- matrix(0,nrow=dim(Y)[1], ncol=dim(Y)[2])
	
	# Anzahl Gruppen pro Bild
	n.groups <- vector(length=dim(Y)[3])
	
	cat(as.character(Sys.time())," Ende Initialisierung Listen \n") 
	
	# In dieser for-Schleife werden fuer jeden 
	#	Zeitpunkt die Hotspots berechnet
	
	# Diese Schleife ist denke ich das, 
	#	was am meisten Rechenzeit braucht
	if(do.clust){
			
		cat(as.character(Sys.time())," Beginn Clustern \n")
		
		for(t in 1:dim(Y)[3]){
			
	#		print(t)
			
			hots[[t]] <- list()
			hots.y01[[t]] <- list()
			dims.l[[t]] <- list()
			
			# Y.hat.cut - Bild zu diesem Zeitpunkt
			A <- Y.hat.cut[,,t]
			
			# Falls in A nur Nullen stehen oder weniger als cut Werte enthalten 
			#	sind
			#	=> Direkt abbrechen
			
			if(sum(A)!=0 & nrow(which(A!=0, arr.ind=TRUE)) >= cut ){
				
				# Bild reinigen (mit Clusteranalyse)
				A.clean <- clean.mult(A, 
															cut = cut, 
															maxdist = maxdist,
															fastclust=fastclust,
															speeditup=speeditup,
															floodclust=floodclust,
															floodCpp=floodCpp,
															t=t,
															pack=pack)
				
				# Clusterergebnis
				clust.erg <- A.clean$Gruppen
				# plot(clust.erg[,1:2], col=clust.erg[,3], pch=19)
				
				# Bereinigtes Bild
				B <- A.clean$Bild
				# image(B)
				
				# 0/1 Bild des bereinigten Bildes
				C <- B
				C[which(C!=0)] <- 1
				
				Y.01[,,t] <- C
				# image(C)
				
																						
				# Summenbild aller bereinigten Bilder
				Y.sum <- Y.sum + B
				
				# Anzahl Gruppen
        if(nrow(clust.erg)==0){
          n.groups[t] <- 0
        }else{        
				  n.groups[t] <- nrow(as.data.frame(table(clust.erg[,3])))
          # besser: length(unique(clust.erg[,3])) ? 
        }
        
        # LB 20.01.14: In naechste if-Klammer verschoben
        #   groups <- as.integer(rownames(table(clust.erg[,3])))

				
				# Falls im bereinigten Bild noch schwarze Pixel sind: 
				# Fuer jeden Hotspot die abzuspeichernden Kennzeichen berechnen
				
				if(sum(C)>0){	
					
				  groups <- as.integer(rownames(table(clust.erg[,3])))
				  # besser: unique(clust.erg[,3]) ? 
          
					# Gruppenbild des bereinigten Bildes (enthaelt als "Signal" den 
					#	Gruppennamen)		
					Group <- C
					for(j in 1:nrow(clust.erg)){
						Group[clust.erg[j,1], clust.erg[j,2]] <- clust.erg[j,3]
					}
					
					# image(Group)
				
					for(i in 1:n.groups[t]){
					
		# 				if(sum(clust.erg==1)!=0){
						# Was soll dieses if() ?
							
							j <- groups[i]
							
							# Dimensionen der Hotpots berechnen
							dims <- hotspot.mult(A=clust.erg, 
														i=j)
							
							# => Ausschnitt aus den Bildern
							
							# 0/1 Bild Ausschnitt
							y01 <- Group[(dims[1,1]-c):(dims[2,1]+c),
											  (dims[1,2]-c):(dims[2,2]+c)]
							
							y01[which(y01!=j)] <- 0
							y01[which(y01==j)] <- 1
							# image(y01)
							
							# 0/1 Bild komplett
							y01.all <- Group
							y01.all[which(y01.all!=j)] <- 0
							y01.all[which(y01.all==j)] <- 1
							# image(y01.all)
			
							
							# Y.hat Ausschnitt
							y.hat <- Y.hat[(dims[1,1]-c):(dims[2,1]+c),
												(dims[1,2]-c):(dims[2,2]+c),t]
							
							# Y.hat - ganzes Bild
							# y.hat.all <- ???
							
							# Y.roh Ausschnitt
							y.roh <- Y[(dims[1,1]-c):(dims[2,1]+c),
										  (dims[1,2]-c):(dims[2,2]+c),t]
							
							# Dimensionen des Rechtecks
							dims.c <- dims + matrix(c(-c,-c,+c,+c), byrow=TRUE, nrow=2)
							
							hots[[t]][[i]] <- list(y01 = y01,
														  y01.all = y01.all,
														  y.hat = y.hat,
														  y.roh = y.roh,
														  dims = dims,
														  dims.c = dims.c,
														  track = NA,
														  species = NA)
							
							# Nur 0/1 - Hotspot abspeichern
							hots.y01[[t]][[i]] <- y01
							
							# Dimensionen abspeichern
							dims.l[[t]][[i]] <- dims
							
		# 				}						
					}
				}
				
				# Falls im bereinigten Bild keine schwarzen Pixel mehr sind: 
				#	Alles auf 0 setzen
				
				if(sum(C)==0){			
					hots[[t]][[1]] <- list(y01 = 0,
												  y01.all = 0,
												  y.hat = 0,
												  y.roh = 0,
												  dims = 0,
												  dims.c = 0,
												  track = NA,
												  species = NA)
					
					hots.y01[[t]][[1]] <- 0
					dims.l[[t]][[1]] <- 0
				}
			}else{
				# Wenn in A  nur Nullen stehen oder weniger als cut Werte != 0 
				#	sind:
				#	Alles auf 0 setzen, wie, wenn in C nur Nullen stehen
				hots[[t]][[1]] <- list(y01 = 0,
											  y01.all = 0,
											  y.hat = 0,
											  y.roh = 0,
											  dims = 0,
											  dims.c = 0,
											  track = NA,
											  species = NA)	
				
				hots.y01[[t]][[1]] <- 0
				dims.l[[t]][[1]] <- 0
			}				
		}
		
		cat(as.character(Sys.time())," Ende Clustern \n")
		
		# Nur alles ausgeben, wenn auch geplottet werden soll, sonst reicht
		#	ein Teil und entlastet den Speicher
		
		if(plot.hots.mult){
			out <- list(Y.01=Y.01,
							Y.sum=Y.sum, 
							Y.hat=Y.hat, 
							Y=Y, 
							df=df,
							Y.hat.zentr=Y.hat.zentr,
							Y.hat.cut=Y.hat.cut,
							n.groups=n.groups,
							hots = hots,
							hots.y01 = hots.y01,
							dims.l=dims.l)
		}else{
			
			# Wenn die Filme gespeicher werden sollen, brauche ich 
			#	zusaetzlich Y 
			if(save.jpg){
				out <- list(Y.01=Y.01,
# 								Y.sum=Y.sum, 
# 								Y.hat=Y.hat, 
								Y=Y, 
								df=df,
# 								Y.hat.zentr=Y.hat.zentr,
# 								Y.hat.cut=Y.hat.cut,
								n.groups=n.groups,
								hots = hots,
								hots.y01 = hots.y01,
								dims.l=dims.l)
			}else{
				out <- list(Y.01=Y.01,
# 								Y.sum=Y.sum, 
# 								Y.hat=Y.hat, 
# 								Y=Y, 
								df=df,
# 								Y.hat.zentr=Y.hat.zentr,
# 								Y.hat.cut=Y.hat.cut,
								n.groups=n.groups,
								hots = hots,
								hots.y01 = hots.y01,
								dims.l=dims.l)
			}
		}
		
		return(out)
	}else{
		
		out <- list(Y.sum=Y.sum, 
								Y.hat=Y.hat, 
								Y=Y, 
								df=df,
								Y.hat.zentr=Y.hat.zentr,
								Y.hat.cut=Y.hat.cut)	
		return(out)
			
	}
}


####################################################################
# Die Funktion plot.hotspots() soll die Ergebnisplots nach Finden 
#	der hotspots erstellen.
#
# Uebergeben werden muss:
#	- hotspots:			Das Ergebnis der Funktion find.hotspot()
#	- pfad:				Der Dateipfad fuer die Plots
#	- t.plot:			Interessante Zeitpunkte
#	- mfrow.tplot:		Anzahl Zeilen und Spalten fuer den Plot - nur t.plot
#	- mfrow.all:		Anzahl Zeilen und Spalten fuer den Plot
#	- width, height, res, quality: Grafikparameter
#	- t.grid:			Vektor, der die Zeitpunkte der analysierten Frames angibt
#	- cellres:			Seitenverhaeltnis der Pixel
#
# Ausgegeben wird:
#	- nur die Plots
#
####################################################################

#' @import pixmap
plot.hotspots <- function(hotspots, pfad, t.plot, 
                          mfrow.tplot, mfrow.all, width=8, 
                          height=8, res=250, quality=100, 
                          t.grid, cellres=1,...){
	
	df <- hotspots$df
	
	##########
	# Y.hat
	##########
	Y.hat <- hotspots$Y.hat
	
	filename <- paste(pfad,"Y_hat_",df[1],"_",df[2],"_",df[3],".jpg", sep="") 	
	jpeg(filename,width=width,height=height, units="in", res=res, quality=quality) 
		pixplot(Y.hat, cellres=cellres, t.grid=t.grid)
	dev.off()
	
	##########
	# Y.roh
	##########
	Y <- hotspots$Y
	
	filename <- paste(pfad,"Y_roh.jpg", sep="") 
	jpeg(filename,width=width,height=height, units="in", res=res, quality=quality) 
		pixplot(Y, cellres=cellres, t.grid=t.grid)
	dev.off()

	##############################
	# Summe ueber die Zeitpunkte
	##############################
	Y.sum <- hotspots$Y.sum
	
	filename <- paste(pfad,"Y_hat_sum.jpg", sep="") 	
	jpeg(filename,width=width,height=height, units="in", res=res, quality=quality) 
		pixplot(Y.sum, cellres=cellres, t.grid="sum")		
	dev.off()
		
	####################
	# Hotspots als jpg
	####################
	
	Liste.y01 <- hotspots$Teil_01
	
	filename <- paste(pfad,"Hotspots.jpg", sep="") 	
	jpeg(filename,width=width,height=height, units="in", res=res, quality=quality) 
	
	par(mfrow=mfrow.tplot)
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				x <- pixmapIndexed(dreh(Liste.y01[[t]]), col=c("white","black"), cellres=1)
				plot(x, main=paste("t = ",t.grid[t], sep=""))
				box()
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
	dev.off()
	
	########################################
	# Hotspots + Y.hat + Y.roh als jpg
	########################################
	
	Liste.y.hat <- hotspots$Teil_Y.hat
	Liste.y.roh <- hotspots$Teil_Y.roh
	
	filename <- paste(pfad,"Hotspots01_hat_roh.jpg", sep="") 
	jpeg(filename,width=width,height=height, units="in", res=res, quality=quality) 
	
	par(mfrow=c(3,length(t.plot)), mar=c(1.7,1.8,1.7,0.2)+0.2)
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				x <- pixmapIndexed(dreh(Liste.y01[[t]]), col=c("white","black"), cellres=1)
				plot(x, main=paste("t = ",t.grid[t], sep=""))
				box()
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				x <- pixmapGrey(dreh(Liste.y.hat[[t]]), cellres=1)
				plot(x, main=paste("t = ",t.grid[t], sep=""))
				box()
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				x <- pixmapGrey(dreh(Liste.y.roh[[t]]), cellres=1)
				plot(x, main=paste("t = ",t.grid[t], sep=""))
				box()
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
	dev.off()
	
	########################################
	# Hotspots aller Zeitpunkte als jpg
	########################################
	
	filename <- paste(pfad,"Hotspots_alleframes.jpg", sep="") 	
	jpeg(filename,width=width,height=height, units="in", res=res, quality=quality) 
	
	par(mfrow=mfrow.all, mar=c(1.7,1.8,1.7,0.2)+0.2)
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				x <- pixmapIndexed(dreh(Liste.y01[[t]]), col=c("white","black"), cellres=1)
				plot(x, main=paste("t = ",t.grid[t], sep=""))
				box()
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}
		}
	dev.off()
	
	################################################################################
	# Hotspots aller Zeitpunkte als jpg - Y.hat - mit ganzem Ausschnitt
	################################################################################

	Liste.y.hat.all <- hotspots$Ganz_Y.hat

	filename <- paste(pfad,"Hotspots_alleframes_Y.hat.jpg", sep="") 
	jpeg(filename,width=width,height=height, units="in", res=res, quality=quality)
	
	par(mfrow=mfrow.all, mar=c(1.7,1.8,1.7,0.2)+0.2)
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				x <- pixmapGrey(dreh(Liste.y.hat.all[[t]]), cellres=cellres)
				plot(x, main=paste("t = ",t.grid[t], sep=""))
				box()
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
	dev.off()

}



####################################################################
# Die Funktion plot.hotspots.ddf() soll die Ergebnisplots nach Finden 
#	der hotspots erstellen. Allerdings beruecksichtigt sie im Gegensatz 
#	zu plot.hotspots() die Trichterform.
#
# Uebergeben werden muss:
#	- hotspots:			Das Ergebnis der Funktion find.hotspot()
#	- pfad:				Der Dateipfad fuer die Plots
#	- t.plot:			Interessante Zeitpunkte
#	- mfrow.tplot:		Anzahl Zeilen und Spalten fuer den Plot - nur t.plot
#	- mfrow.all:		Anzahl Zeilen und Spalten fuer den Plot
#	- bw:				1: Schwarz und weiss vertauscht
#	- xs, ys:			Koordinaten der Punkte
#	- t.grid:			Vektor, der die Zeitpunkte der analysierten Frames angibt
# 	- cex:				Groesse der Punkte
#	- width, height, res, quality: Grafikparameter
#
# Ausgegeben wird:
# 	- nur die Plots
#
####################################################################

plot.hotspots.ddf <- function(hotspots, pfad, t.plot, mfrow.tplot, mfrow.all, bw=1, xs, ys, t.grid, cex=1, width=8, height=8, res=250, quality=100){
	
	df <- hotspots$df
	dims <- hotspots$Dims
	
	##########
	# Y.hat
	##########
	Y.hat <- hotspots$Y.hat

	####################################################################	
	# ratio streckt bzw. staucht das Bild in die Breite. Das ist 
	#	noetig, damit die einzelnen Plots quadratisch sind und somit 
	#	die Seitenverhaeltnisse stimmen.
	####################################################################
	ratio <- mfrow.all[2]/mfrow.all[1]
	
	filename <- paste(pfad,"Y_hat_",df[1],"_",df[2],"_",df[3],".jpg", sep="") 	
	jpeg(filename,width=width*ratio,height=height, units="in", res=res, quality=quality) 
		par(mfrow=mfrow.all, mar=c(1.7,1.8,1.7,0.2)+0.2)
		for(t in t.plot) {
			plot.ddf(Y.hat[,,t], bw=bw, xs=xs, ys=ys, t=t.grid[t])		
		}
		
	dev.off()
	
	##########
	# Y.roh
	##########
	Y <- hotspots$Y

	ratio <- mfrow.all[2]/mfrow.all[1]
		
	filename <- paste(pfad,"Y_roh.jpg", sep="") 
	jpeg(filename,width=width*ratio,height=height, units="in", res=res, quality=quality) 
		par(mfrow=mfrow.all, mar=c(1.7,1.8,1.7,0.2)+0.2)
			for(t in t.plot) {
				plot.ddf(Y[,,t], bw=bw, xs=xs, ys=ys, t=t.grid[t])
			}
			
	dev.off()
	
	###############
	# Y.hat.zentr
	###############
	
	Y.hat.zentr <- hotspots$Y.hat.zentr
	
	filename <- paste(pfad,"Y_hat_zentr.jpg", sep="") 
	jpeg(filename,width=width*ratio,height=height, units="in", res=res, 
		  quality=quality) 
	par(mfrow=mfrow.all, mar=c(1.7,1.8,1.7,0.2)+0.2)
	for(t in t.plot) {
		plot.ddf(Y.hat.zentr[,,t], bw=bw, xs=xs, ys=ys, t=t.grid[t])
	}
	
	dev.off()

	###############
	# Y.hat.cut
	###############
	
	Y.hat.cut <- hotspots$Y.hat.cut
	
	filename <- paste(pfad,"Y_hat_cut.jpg", sep="") 
	jpeg(filename,width=width*ratio,height=height, units="in", res=res, 
		  quality=quality) 
	par(mfrow=mfrow.all, mar=c(1.7,1.8,1.7,0.2)+0.2)
	for(t in t.plot) {
		plot.ddf(Y.hat.cut[,,t], bw=bw, xs=xs, ys=ys, t=t.grid[t])
	}
	
	dev.off()
	
	
	##############################
	# Summe ueber die Zeitpunkte
	##############################
	Y.sum <- hotspots$Y.sum
	
	filename <- paste(pfad,"Y_hat_sum.jpg", sep="") 	
	jpeg(filename,width=width,height=height, units="in", res=res, quality=quality) 
		plot.ddf(Y.sum, bw=bw, xs=xs, ys=ys, t="sum")		
	dev.off()
	
	####################
	# Hotspots als jpg
	####################
	
	Liste.y01 <- hotspots$Teil_01

	ratio <- mfrow.tplot[2]/mfrow.tplot[1]
	
	filename <- paste(pfad,"Hotspots.jpg", sep="") 	
	jpeg(filename,width=width*ratio,height=height, units="in", res=res, quality=quality) 
	
	par(mfrow=mfrow.tplot, mar=c(1.7,1.8,1.7,0.2)+0.2)
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				
				# x und y Koordinaten der zu plottenden Punkte
				xs.plot <- xs[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				ys.plot <- ys[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				plot.ddf(Liste.y01[[t]], bw=bw, xs=xs.plot, ys=ys.plot, bin=TRUE, cex=cex, t=t.grid[t])
				
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
	dev.off()
	
	########################################
	# Hotspots + Y.hat + Y.roh als jpg
	########################################
	
	Liste.y.hat <- hotspots$Teil_Y.hat
	Liste.y.roh <- hotspots$Teil_Y.roh
	
	ratio <- length(t.plot)/3
	
	filename <- paste(pfad,"Hotspots01_hat_roh.jpg", sep="") 
	jpeg(filename,width=width*ratio,height=height, units="in", res=res, quality=quality) 
	
	par(mfrow=c(3,length(t.plot)), mar=c(1.7,1.8,1.7,0.2)+0.2)
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				
				# x und y Koordinaten der zu plottenden Punkte
				xs.plot <- xs[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				ys.plot <- ys[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				plot.ddf(Liste.y01[[t]], bw=bw, xs=xs.plot, ys=ys.plot, bin=TRUE, cex=cex, t=t.grid[t])

			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
		
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				
				# x und y Koordinaten der zu plottenden Punkte
				xs.plot <- xs[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				ys.plot <- ys[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				plot.ddf(Liste.y.hat[[t]], bw=bw, xs=xs.plot, ys=ys.plot, bin=FALSE, cex=cex, t=t.grid[t])
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
		
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				
				# x und y Koordinaten der zu plottenden Punkte
				xs.plot <- xs[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				ys.plot <- ys[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				plot.ddf(Liste.y.roh[[t]], bw=bw, xs=xs.plot, ys=ys.plot, bin=FALSE, cex=cex, t=t.grid[t])
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
	dev.off()
	
	########################################
	# Hotspots aller Zeitpunkte als jpg
	########################################
	
	ratio <- mfrow.all[2]/mfrow.all[1]

	filename <- paste(pfad,"Hotspots_alleframes.jpg", sep="") 	
	jpeg(filename,width=width*ratio,height=height, units="in", res=res, quality=quality) 
	
	par(mfrow=mfrow.all, mar=c(1.7,1.8,1.7,0.2)+0.2)
		for(t in 1:t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				
				# x und y Koordinaten der zu plottenden Punkte
				xs.plot <- xs[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				ys.plot <- ys[dims[[t]][1,1]:dims[[t]][2,1],dims[[t]][1,2]:dims[[t]][2,2]]
				plot.ddf(Liste.y01[[t]], bw=bw, xs=xs.plot, ys=ys.plot, bin=TRUE, cex=cex, t=t.grid[t])
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}
		}
	dev.off()
	
	################################################################################
	# Hotspots aller Zeitpunkte als jpg - Y.hat - mit ganzem Ausschnitt
	################################################################################

	Liste.y.hat.all <- hotspots$Ganz_Y.hat

	ratio <- mfrow.all[2]/mfrow.all[1]
	
	filename <- paste(pfad,"Hotspots_alleframes_Y.hat.jpg", sep="") 
	jpeg(filename,width=width*ratio,height=height, units="in", res=res, quality=quality)
	
	par(mfrow=mfrow.all, mar=c(1.7,1.8,1.7,0.2)+0.2)
		for(t in t.plot) {
			if(is.matrix(Liste.y01[[t]])) {
				
				plot.ddf(Liste.y.hat.all[[t]], bw=bw, xs=xs, ys=ys, bin=FALSE, t=t.grid[t])
			}		
			if(!is.matrix(Liste.y01[[t]])) {
				plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}	
		}
	dev.off()

}


####################################################################
# Die Funktion plot.hotspots.mult() soll die Ergebnisplots nach Finden 
#	der hotspots erstellen. 
#	Variante von plot.hotspots.ddf(), die mehrere Hotspots pro Bild
#		zulaesst und dementsprechend als Input das Ergebnis von
#		find.mult.hotspots() benoetigt
#
# Uebergeben werden muss:
#	- hotspots:			Das Ergebnis der Funktion find.mult.hotspot()
#	- pfad:				Der Dateipfad fuer die Plots
#	- t.plot:			Interessante Zeitpunkte
#	- mfrow.tplot:		Anzahl Zeilen und Spalten fuer den Plot - nur t.plot
#	- mfrow.all:		Anzahl Zeilen und Spalten fuer den Plot
#	- bw:					1: Schwarz und weiss vertauscht
#	- xs, ys:			Koordinaten der Punkte
#	- t.grid:			Vektor, der die Zeitpunkte der analysierten Frames angibt
# 	- cex:				Groesse der Punkte
#	- width, height, res, quality: Grafikparameter
#
# Ausgegeben wird:
# 	- nur die Plots
#
####################################################################

plot.hotspots.mult <- function(hotspots, 
															 pfad, 
															 t.plot, 
															 mfrow.tplot, 
															 mfrow.all, 
															 bw=1, 
															 xs, 
															 ys, 
															 t.grid, 
															 cex=1, 
															 width=NULL, 
															 height=NULL, 
															 res=250, # 250
															 quality=100, #100,
															 plot.hots=TRUE
){

	######################################################################
	# Update 6.6.13: 1:n.frames durch t.plot ersetzt, damit nicht immer
	#	alles geplottet wird
	######################################################################
	
	
	# Groesse des Plots festlegen
	if(is.null(width)){
		width.tplot  <- 2*mfrow.tplot[2] 
		width.all  <- 2*mfrow.all[2]
	}else{
		width.tplot  <- width 
		width.all  <- width
	}
	
	if(is.null(height)){
		height.tplot  <- 2*mfrow.tplot[1] 
		height.all  <- 2*mfrow.all[1]
	}else{
		height.tplot  <- height 
		height.all  <- height
	}
	
	
	
	df <- hotspots$df
	# dims <- hotspots$Dims
	# geht nicht mehr...

	
# 	n.frames <- dim(Y)[3]
	# Update 15.04.13 nach Hinweis Carolin:
	n.frames <- dim(hotspots$Y)[3]
	
	
	##########
	# Y.hat
	##########
	Y.hat <- hotspots$Y.hat
	
	####################################################################	
	# ratio streckt bzw. staucht das Bild in die Breite. Das ist 
	#	noetig, damit die einzelnen Plots quadratisch sind und somit 
	#	die Seitenverhaeltnisse stimmen.
	####################################################################
	# ratio <- mfrow.all[2]/mfrow.all[1]
	
	filename <- paste(pfad,"Y_hat_",df[1],"_",df[2],"_",df[3],".jpg", sep="") 	
	jpeg(filename,width=width.all,height=height.all, units="in", 
		  res=res, quality=quality) 
		par(mfrow=mfrow.all, mar=c(1.7,1.8,2.7,0.2)+0.2)
			for(t in t.plot) {
				plot.ddf(Y.hat[,,t], bw=bw, xs=xs, ys=ys, t=t.grid[t], 
								 cex.main=3)		
			}
	
	dev.off()
	
	##########
	# Y.roh
	##########
	Y <- hotspots$Y
	
	# ratio <- mfrow.all[2]/mfrow.all[1]
	
	filename <- paste(pfad,"Y_roh.jpg", sep="") 
	jpeg(filename,width=width.all,height=height.all, units="in", res=res, quality=quality) 
		par(mfrow=mfrow.all, mar=c(1.7,1.8,2.7,0.2)+0.2)
			for(t in t.plot) {
				plot.ddf(Y[,,t], bw=bw, xs=xs, ys=ys, t=t.grid[t], 
								 cex.main=3)		
			}
	
	dev.off()
	
	###############
	# Y.hat.zentr
	###############
	
	Y.hat.zentr <- hotspots$Y.hat.zentr
	
	filename <- paste(pfad,"Y_hat_zentr.jpg", sep="") 
	jpeg(filename,width=width.all,height=height.all, units="in", res=res, 
		  quality=quality) 
		par(mfrow=mfrow.all, mar=c(1.7,1.8,2.7,0.2)+0.2)
			for(t in t.plot) {
				plot.ddf(Y.hat.zentr[,,t], bw=bw, xs=xs, ys=ys, t=t.grid[t], 
								 cex.main=3)		
			}
	
	dev.off()
	
	###############
	# Y.hat.cut
	###############
	
	Y.hat.cut <- hotspots$Y.hat.cut
	
	filename <- paste(pfad,"Y_hat_cut.jpg", sep="") 
	jpeg(filename,width=width.all,height=height.all, units="in", res=res, 
		  quality=quality) 
		par(mfrow=mfrow.all, mar=c(1.7,1.8,2.7,0.2)+0.2)
			for(t in t.plot) {
				plot.ddf(Y.hat.cut[,,t], bw=bw, xs=xs, ys=ys, t=t.grid[t], 
								 cex.main=3)		
			}
	
	dev.off()
	
	if(plot.hots){
		
		##############################
		# Summe ueber die Zeitpunkte
		##############################
		Y.sum <- hotspots$Y.sum
		
		filename <- paste(pfad,"Y_hat_sum.jpg", sep="") 	
		jpeg(filename,width=8,height=8, units="in", res=res, quality=quality) 
			plot.ddf(Y.sum, bw=bw, xs=xs, ys=ys, t="sum", 
							 cex.main=3)		
		dev.off()
		

			
		####################
		# Hotspots als jpg
		####################
		
		hots <- hotspots$hots
		
		# Anzahl Zeilen und Spalten fuer den Plot bestimmen, haengt von
		#	Anzahl Zeitpunkte und Hotspots pro Zeitpunkt ab
		
		n.groups <- hotspots$n.groups
		ng.plot <- n.groups
		ng.plot[which(ng.plot==0)] <- 1
	# 	nRow <- ceiling(sqrt(sum(ng.plot[t.plot])))
	# 	
	# 	if(nRow==0){
	# 		nRow <- 1
	# 	}
		
		mfrow.tplot <- n2mfrow(sum(ng.plot[t.plot]))
			
		ratio <- mfrow.tplot[2]/mfrow.tplot[1]
		
		filename <- paste(pfad,"Hotspots.jpg", sep="") 	
		jpeg(filename,width=2*mfrow.tplot[2],height=2*mfrow.tplot[1], units="in", 
			  res=res, quality=quality) 
		
		par(mfrow=mfrow.tplot, mar=c(1.7,1.8,2.7,0.2)+0.2)
		for(t in t.plot) {
			if(is.matrix(hots[[t]][[1]]$y01)) {
				
				# Schleife ueber alle Hotspots in dem Bild
				for(i in 1:n.groups[t]){
					
					dims <- hots[[t]][[i]]$dims
					
					# x und y Koordinaten der zu plottenden Punkte
					xs.plot <- xs[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
					ys.plot <- ys[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
					plot.ddf(hots[[t]][[i]]$y01, bw=bw, xs=xs.plot, ys=ys.plot, 
								bin=TRUE, cex=cex, t=t.grid[t],
									 cex.main=3)
				}		
			}		
			if(!is.matrix(hots[[t]][[1]]$y01)) {
				plot(0,0,col="white", 
					  main=paste("t = ",t.grid[t],": No Fish",sep=""),
					  axes=FALSE, xlab="", ylab="", cex.main=1.5)
				box()
			}	
		}
		dev.off()
		
		########################################
		# Hotspots + Y.hat + Y.roh als jpg
		########################################
		
	# 	# Liste.y.hat <- hotspots$Teil_Y.hat
	# 	# Liste.y.roh <- hotspots$Teil_Y.roh
	# 	
	# 	ratio <- sum(n.groups[t.plot])/3
	# 	
	# 	filename <- paste(pfad,"Hotspots01_hat_roh.jpg", sep="") 
	# 	jpeg(filename,width=8*ratio,height=8, units="in", 
	# 		  res=res, quality=quality) 
	# 	
	# 	par(mfrow=c(3,sum(n.groups[t.plot])), mar=c(1.7,1.8,1.7,0.2)+0.2)
	# 	
	# 	# Y.01
	# 	for(t in t.plot) {
	# 		if(is.matrix(hots[[t]][[1]]$y01)) {
	# 			
	# 			# Schleife ueber alle Hotspots in dem Bild
	# 			for(i in 1:n.groups[t]){
	# 				
	# 				dims <- hots[[t]][[i]]$dims
	# 				
	# 				# x und y Koordinaten der zu plottenden Punkte
	# 				xs.plot <- xs[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
	# 				ys.plot <- ys[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]				
	# 				plot.ddf(hots[[t]][[i]]$y01, bw=bw, xs=xs.plot, 
	# 							ys=ys.plot, bin=TRUE, cex=cex, 
	# 							t=t.grid[t])
	# 			}
	# 		}		
	# 		if(!is.matrix(hots[[t]][[1]]$y01)) {
	# 			plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
	# 			box()
	# 		}	
	# 	}
	# 	
	# 	# Y.hat
	# 	for(t in t.plot) {
	# 		if(is.matrix(hots[[t]][[1]]$y01)) {
	# 			
	# 			# Schleife ueber alle Hotspots in dem Bild
	# 			for(i in 1:n.groups[t]){
	# 				
	# 				dims <- hots[[t]][[i]]$dims
	# 				
	# 				# x und y Koordinaten der zu plottenden Punkte
	# 				xs.plot <- xs[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
	# 				ys.plot <- ys[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]				
	# 				plot.ddf(hots[[t]][[i]]$y.hat, bw=bw, xs=xs.plot, 
	# 							ys=ys.plot, bin=FALSE, cex=cex, 
	# 							t=t.grid[t])
	# 			}
	# 		}		
	# 		if(!is.matrix(hots[[t]][[1]]$y01)) {
	# 			plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
	# 			box()
	# 		}	
	# 	}
	# 	
	# 	# Y.roh
	# 	for(t in t.plot) {
	# 		if(is.matrix(hots[[t]][[1]]$y01)) {
	# 			
	# 			# Schleife ueber alle Hotspots in dem Bild
	# 			for(i in 1:n.groups[t]){
	# 				
	# 				dims <- hots[[t]][[i]]$dims
	# 				
	# 				# x und y Koordinaten der zu plottenden Punkte
	# 				xs.plot <- xs[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
	# 				ys.plot <- ys[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]				
	# 				plot.ddf(hots[[t]][[i]]$y.roh, bw=bw, xs=xs.plot, 
	# 							ys=ys.plot, bin=FALSE, cex=cex, 
	# 							t=t.grid[t])
	# 			}
	# 		}		
	# 		if(!is.matrix(hots[[t]][[1]]$y01)) {
	# 			plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
	# 			box()
	# 		}	
	# 	}
	# 	dev.off()
		
	# 	########################################
	# 	# Hotspots aller Zeitpunkte als jpg
	# 	########################################
	# 	
	# 	nRow.all <- ceiling(sqrt(sum(n.groups)))
	# 	
	# 	if(nRow.all==0){
	# 		nRow.all <- 1
	# 	}
	# 	
	# 	mfrow.all.g <- c(nRow.all, nRow.all)
	# 	
	# 	# ratio <- mfrow.all.g[2]/mfrow.all.g[1]
	# 	
	# 	filename <- paste(pfad,"Hotspots_alleframes.jpg", sep="") 	
	# 	jpeg(filename,width=2*nRow.all,height=2*nRow.all, units="in", 
	# 		  res=res, quality=quality) 
	# 	
	# 	# Update 17.06.13: mfrow cleverer uebergeben
	# 	par(mfrow=n2mfrow(sum(n.groups) + length(which(n.groups==0))), mar=c(1.7,1.8,1.7,0.2)+0.2)
	# # 	par(mfrow=mfrow.all.g, mar=c(1.7,1.8,1.7,0.2)+0.2)
	# 	
	# 	# Y.01
	# 	for(t in 1:n.frames) {
	# 		if(is.matrix(hots[[t]][[1]]$y01)) {
	# 			
	# 			# Schleife ueber alle Hotspots in dem Bild
	# 			for(i in 1:n.groups[t]){
	# 				
	# 				dims <- hots[[t]][[i]]$dims
	# 				
	# 				# x und y Koordinaten der zu plottenden Punkte
	# 				xs.plot <- xs[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
	# 				ys.plot <- ys[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]				
	# 				plot.ddf(hots[[t]][[i]]$y01, bw=bw, xs=xs.plot, 
	# 							ys=ys.plot, bin=TRUE, cex=cex, 
	# 							t=t.grid[t])
	# 			}
	# 		}		
	# 		if(!is.matrix(hots[[t]][[1]]$y01)) {
	# 			plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
	# 			box()
	# 		}
	# 	}
	# 	dev.off()
		
	# 	################################################################################
	# 	# Hotspots aller Zeitpunkte als jpg - Y.hat - mit ganzem Ausschnitt
	# 	################################################################################
	# 	
	# 	# ratio <- mfrow.all[2]/mfrow.all[1]
	# 	
	# 	filename <- paste(pfad,"Hotspots_alleframes_Y_01.jpg", sep="") 
	# 	jpeg(filename,width=width.all,height=height.all, units="in", 
	# 		  res=res, quality=quality)
	# 	
	# 	par(mfrow=mfrow.all, mar=c(1.7,1.8,1.7,0.2)+0.2)
	# 	for(t in t.plot) {
	# 		if(is.matrix(hots[[t]][[1]]$y01)) {
	# 				
	# 			plot.ddf(hotspots$Y.01[,,t], bw=bw, xs=xs, 
	# 							ys=ys, bin=FALSE, cex=cex, 
	# 							t=t.grid[t])
	# 
	# 		}		
	# 		if(!is.matrix(hots[[t]][[1]]$y01)) {
	# 			plot(0,0,col="white", main=paste("t = ",t.grid[t],": No Fish",sep=""),axes=FALSE, xlab="", ylab="")
	# 			box()
	# 		}	
	# 	}
	# 	dev.off()
	}
	
}
































