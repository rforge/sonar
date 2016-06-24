###########################################################################
# Fischprojekt
# Funktionen zum Tracken der Hotspots, also zur Zuordnung der Hotspots
#		zu den richtigen Obkjekten
# Datum der ersten Version: 21.08.12
# Letzte Änderung: 22.10.12
# Autor: Ludwig Bothmann
###########################################################################

# Für Funktion cart()

# source("Funktionen_Merkmale_bestimmen.R")

####################################################################
# Die Funktion track.func() soll die Hotspots den richtigen Objekten 
#		zuordnen, so dass diese über die Zeit verfolgbar sind.
#
# Da es hier sein kann, dass ein Hotspot mehreren Objekten zugeordnet
#		wird und das unter Umständen sehr viele werden, schreibe ich 
#		noch eine Funktion split.tracks(), welche die Doppel- oder
#		Mehrfach-Zuordnungen wieder auflösen bzw. reduzieren woll.
#
# Übergeben werden muss:
#	- hots.mult:	Ergebnis der Funktion find.mult.hotspots()
#	- t.plot:		Zeitpunkte, die von Interesse sind
#	- max.dist.cs:	Maximaler Abstand zweier Schwerpunkte, damit die 
#							entsprechenden Hotspots dem gleichen Objekt
#							zugeordnet werden. Pro Zehntelsekunde.
#							Default = 0.05 m
#	- rep.time:		Zeit zwischen Frames in Zehntelsekunden, im Normalfall
#							entspricht das der Anzahl Frames dazwischen
#	- cart.coord:	Kartesische Koordinaten der Pixel
#
# Ausgegeben wird: Liste mit
#	- tracker:		Liste mit so viel Elementen wie Hotspots, Eintrag 
#							dem (oder den) Objekt(en), dem der Hotspot zugeordnet 
#							wurde
#	- anz.obj:		Anzahl verschiedener Objekte
#	- matcher:		Matcher-Matrix
#	- t.plot:		siehe oben
#	
#
####################################################################


track.func <- function(hots.mult,
							  t.plot,
							  max.dist.cs = 0.05, # gewählt nach Beispiel 0.25 in 5 Frames
							  rep.time,
							  cart.coord
							  )
{
	
	########################################################################
	# Maximale Distanz zweier Schwerpunkte zu aufeinanderfolgenden 
	#	Zeitpunkten, so dass sie zu einem Objekt gezählt werden.
	#	Berechnet sich als Konstante*Zeitabstand
	#	Konstante eventuell nach Geschwindigkeit der Fische wählen
	#	=> s = v*t
	#	v in Strecke/Zehntelsekunde (in Metern)
	#	t in Zehntelsekunden (also # Frames zwischen 2 Zeitpunkte,
	#		da die normal in Zehntelsekunden aufgezeichnet werden)
	########################################################################
	max.dist <- max.dist.cs * rep.time
	
	
	# Anzahl der Hotspots pro Zeitpunkt
	n.groups <- hots.mult$n.groups
	
	# Matrix initialisieren, die die Hotspots einzelnen Objekten / 
	#	Individuen zuordnet
	matcher <- cbind(1:sum(n.groups[t.plot]),
						  rep(t.plot, n.groups[t.plot]),
						  0,
						  0,
						  0)
	
	colnames(matcher) <- c("id", "t", "Gruppe in t","s.x","s.y")
	
	# 3. Spalte: Gruppe innerhalb des Zeitpunkts t
	# => Kann direkt befüllt werden, es geht einfach immer bei 1 los
	# Achtung: Das müsste schon so stimmen, bei der Clusterung ist es zwar 
	#	nicht immer 1:x, aber ich speichere es dann so ab. Dort:
	#	hots[[t]][[i]] statt hots[[t]][[j]].
	
	matcher[1,3] <- 1
	
	# Falls es überhaupt mehr als eine ID gibt
	
	if(nrow(matcher)>1){
		for(i in 2:nrow(matcher)){
			
			if(matcher[i,2]==matcher[(i-1),2]){
				matcher[i,3] <- matcher[(i-1),3]+1
			}else{
				matcher[i,3] <- 1
			}
		}
	}
		
	
	###################
	# matcher fertig
	###################
	
	# Liste "tracker" erstellen, die für jeden Hotspot später die Objekte
	#	enthält, zu denen er passen könnte
	
	tracker <- vector("list",nrow(matcher))
	
	#############################################
	# Schwerpunkt jedes Hotspots berechnen
	#############################################
	
	# Element hots rausziehen
	hots <- hots.mult$hots
	
	for(i in 1:nrow(matcher)){
		
		# Zeitpunkt des Hotspots
		t <- matcher[i,2]
		
		# Gruppe / ID innerhalb des Zeitpunkts
		j <- matcher[i,3] 
		
		# 0/1 - Bild des Hotspots - Ausschnitt
		y <- hots[[t]][[j]]$y01
		
		if(is.matrix(y)){
		
			# bisher wie in Funktionen_Merkmale_bestimmen, dort wird erst der
			#	Schwerpunkt auf den Indizes bestimmt und dann umgerechnet,
			#	sollte ich das ändern?
		
			# Dimensionen des Hotspots
			dims <- hots[[t]][[j]]$dims
				
			matcher[i,4:5] <- schwerpunkt(y, dims, cart.coord)
						
		}			
	}
	
	#############################################
	# Schwerpunkt berechnen - fertig
	#############################################
	
	##################################################
	# Jetzt kommt das eigentliche Tracking
	#	Alle Zeitpunkte durchgehen und immer sehen,
	#	ob der Hotspot zu einem oder mehreren des
	#	letzten Zeitpunkts passt.
	##################################################
	
	t.max <- length(t.plot)
	
	for(i in 1:t.max){
		
		# wahrer Zeitpunkt innerhalb der analysierten
		t <- t.plot[i]
		
		# Beim ersten Zeitpunkt bekommen alle eine neue Nummer
		if(i==1){
			
			# ID's der Hotspots zum ersten Zeitpunkt
			ids <- which(matcher[,2]==t)
			
			for(j in ids){
				tracker[[j]] <- j
			}		
			
			# Anzahl verschiedener Objekte
			anz.obj <- length(ids)
		
		}else{
			
			# Ab dem zweiten Zeitpunkt müssen die Schwerpunkte immer mit denen
			#	des vorherigen Zeitpunkts verglichen werden.
			
			# Zeitpunkt davor
			t_1 <- t.plot[i-1]
			
			# Indizes der Hotspots zum Zeitpunkt t und t-1
			ind.t <- which(matcher[,2]==t)
			ind.t_1 <- which(matcher[,2]==t_1)
			
			# Falls zum Zeitpunkt t kein Hotspot erkannt wurde:
			#	direkt zum nächsten Zeitpunkt springen
			
			if(length(ind.t)!=0){
				
				# Falls zum Zeitpunkt t-1 kein Hotspot erkannt wurde:
				#	neue Nummern vergeben
				if(length(ind.t_1)==0){
				
					for(j in 1:length(ind.t)){
						
						# id des j-ten Hotspots zum Zeitpunkt t
						id <- ind.t[j]
						
						# Zähler um eins erhöhen
						anz.obj <- anz.obj + 1
						
						# und übergeben
						tracker[[id]] <-  anz.obj
					}
					
					
				}else{
				
					# Distanzmatrix:
					dist.matrix <- as.matrix(as.matrix(dist(matcher[c(ind.t_1,
															ind.t),
														 4:5])
													)[-(1:length(ind.t_1)),1:length(ind.t_1)])
					
					# TRUE, wenn Abstand kleiner gleich max.dist
					track.matrix <- dist.matrix<=max.dist
					
					# Alle Hotspots des Zeitpunkts t durchgehen und entweder
					#	Objektnummern des vorherigen Zeitpunktes nehmen oder
					#	neue Objektnummern vergeben.
					
					for(j in 1:length(ind.t)){
						
						# id des j-ten Hotspots zum Zeitpunkt t
						id <- ind.t[j]
						
						
						# Falls mindestens einmal TRUE vorkommt wird Nummer von dem 
						#	oder den Hotspot(s) des vorherigen Zeitpunktes genommen
						if(sum(track.matrix[j,]) >= 1){
							
							# Hotspots bei t-1 zu denen der Hotspot passt
							merge.obj <- which(track.matrix[j,])
							
							# Ausschnitt von matcher, nur Zeilen der gesuchten ID's
							matcher.t_1 <- matcher[which(matcher[,2]==t_1),]
							
							# id's dieser Hotspots (exakte Zuordnung der Hotspots,
							#		gleiche id, wie Element in tracker)
							if(is.vector(matcher.t_1)){ # Falls matcher.t_1 Vektor)
								id.merge <- matcher.t_1[merge.obj]
							}else{
								id.merge <- matcher.t_1[merge.obj,1]
							}
							
							# Vektor, in den die Objektnummern reingeschrieben werden
							track.id <- NULL
							
							# merge.obj durchgehen und alle Objektnummern rausziehen von
							#	den entsprechenden Hotspots des vorherigen Zeitpunktes
							for(k in 1:length(merge.obj)){
							
								track.id <- union(track.id,tracker[[id.merge[k]]])
								
							}
							
							# Falls track.id nicht NULL ist 
							#	=> übergeben, sonst bleibt das Element
							#	tracker[[id]] leer. Dürfte aber eigentlich
							#	eh nicht vorkommen...
							
							if(!is.null(track.id)){
								tracker[[id]] <-  track.id
							}	
								
						}else{	# Sonst wird eine neue Nummer vergeben
							
							# Zähler um eins erhöhen
							anz.obj <- anz.obj + 1
							
							# und übergeben
							tracker[[id]] <-  anz.obj
							
						}				
					}
				}	
			}
		}	
	}	
	
	#######################
	# Tracking fertig
	#######################
	
	# Zurückgegeben wird tracker
	
	out <- list(tracker=tracker, 
					anz.obj=anz.obj,
					matcher=matcher,
					t.plot=t.plot)
	
	return(out)
	
}
	

####################################################################
# Die Funktion schwerpunkt() soll den Schwerpunkt eines Hotspots bestimmen
#
# Übergeben werden muss:
#	- y:				Das 0/1 - Bild
#	- dims: 			Die Dimension des 0/1 Bildes, also wo y im ganzen Bild liegt
#	- cart.coord:	Kartesische Koordinaten jedes Pixels
#
# Ausgegeben wird:
#	- schwerp
#
####################################################################

schwerpunkt <- function(y, 
								dims,
								cart.coord){
	
	# Anzahl Pixel
	anz <- sum(y)

	# Breite und Länge des Rechtecks in Pixeln
	l.dim <- apply(dims,2,diff)[1]+1
	b.dim <- apply(dims,2,diff)[2]+1
	
	# Marginale Schwerpunkte absolut als Indizes:
	
	s.y <- sum(apply(y, 1, sum) * 1:l.dim) / anz + dims[1,1] - 1
	s.x <- sum(apply(y, 2, sum) * 1:b.dim) / anz + dims[1,2] - 1
	
	# Schwerpunkt kartesisch			
	schwerp <- cart(c(s.y, s.x),cart.coord)
	
	return(schwerp)
	
	}
	

####################################################################
# Die Funktion plot.track() soll ein Track-Objekt plotten:
#		Großes Schachbrett mit #Zeilen = #Objekte, #Spalten = #Zeitpunkte
#
# Übergeben werden muss:
#	- track.ddf:	Ergebnis von track.func()
#	- hots.mult:	Ergebnis von find.mult.hotspots()
#	- pfad:			Dateipfad wohin geplottet werden soll
#	- filename:		Name der Datei
#	- cex:			Plotparameter, default=1
#	- t.grid:		Gitter der Frame-Zeitpunkte
#	- cart.coord:	Kartesische Koordinaten der Pixel
#
# Ausgegeben wird:
#	- nichts, nur plot
#
####################################################################

plot.track <- function(track.ddf,
							  hots.mult,
							  pfad,
							  filename=NULL,
							  cex=1,
							  t.grid,
							  cart.coord){

	xs <- cart.coord[,,1]
	ys <- cart.coord[,,2]
	
	anz.obj <- track.ddf$anz.obj
	tracker <- track.ddf$tracker
	matcher <- track.ddf$matcher
	t.plot <-  track.ddf$t.plot
	
	hots <- hots.mult$hots
	
	# Anzahl Zeilen und Spalten bestimmen (aber maximal 200 Spalten)
	nRow <- anz.obj
	nCol <- min(length(t.plot),100)
	
	# Grafikparameter
	height <- nRow * 1
	width <-  nCol * 1
	
	# Dateiname
	if(is.null(filename)){
		filename  <- paste(pfad,"Tracking.jpg", sep="")
	}else{
		filename  <- paste(pfad,filename, sep="")	
	}
	
	jpeg(filename,
		  width=width,
		  height=height, 
		  units="in", 
		  res=200, 
		  quality=100)
	
	# Plot initialisieren
	par(mfrow=c(nRow, nCol), mar=c(1.7,1.8,1.7,0.2)+0.2)
	
	for(i in 1:nRow){	
		for(j in 1:nCol){
			
			# Zeitpunkt t in t.plot
			t <- t.plot[j]
			
			# ID's (Hotspots) die diesem Objekt zugeordnet sind
			# ids <- which(is.element(tracker,i))
			
			ids <- NULL
			
			# n.objects enthält die Anzahl Objekte zu diesem Hotspot
			# Achtung, ist was anderes als Anzahl Hotpots zu einem Objekt
			
			n.objects <- NULL
			
			for(k in 1:length(tracker)){
				if(length(intersect(tracker[[k]],i))>0){
					ids <- union(ids, k)	
					n.objects[k] <- length(tracker[[k]])
				}
			}
			
			# davon die zum Zeitpunkt t 
			if(length(ids)==1){
				ids.t <- ids[which(matcher[ids,][2]==t)]
			}else{			
				ids.t <- ids[which(matcher[ids,][,2]==t)]		
			}
			
			# Wenn das Objekt i zum Zeitpunkt t nicht da ist: Leerer Plot 
			if(length(ids.t)==0){
				plot(0,0,col="white", 
					  main=paste("t = ",t.grid[t],": No Fish",sep=""),
					  axes=FALSE, xlab="", ylab="")
				box()			
			}
			# Wenn dem Objekt i zum Zeitpunkt t ein Hotspots zugeordnet wird:
			#		Plotten
			if(length(ids.t)==1){
				
				# ID des Hotspots innerhalb des Zeitpunkts t
				i.t <- matcher[ids.t,3]
			
				dims <- hots[[t]][[i.t]]$dims
				
				# x und y Koordinaten der zu plottenden Punkte
				xs.plot <- xs[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
				ys.plot <- ys[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
				plot.ddf(hots[[t]][[i.t]]$y01, bw=0, xs=xs.plot, ys=ys.plot, 
							bin=TRUE, cex=cex, t=t.grid[t], col=i)	
				
				# ID des Hotspots
				k <- matcher[ids.t,1]
				if(n.objects[k]>1) {
				
					legend("topright", paste("Achtung!", n.objects[k], "Objekte"))				
					
				}
				
			}
			
			# Wenn dem Objekt i zum Zeitpunkt t mehrere Hotspots zugeordnet 
			#	werden: Ersten Hotspot plotten und Anmerkung, dass und wie 
			#	viele es eigentlich wären.
			if(length(ids.t)>1){
				
				# ID des Hotspots innerhalb des Zeitpunkts t
				i.t.mult <- matcher[ids.t,3]
				
				# Anzahl Hotspots
				n.hots <- length(i.t.mult)
				
				# Kleinster Index
				i.t <- min(i.t.mult)
				
				dims <- hots[[t]][[i.t]]$dims
				
				# x und y Koordinaten der zu plottenden Punkte
				xs.plot <- xs[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
				ys.plot <- ys[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
				plot.ddf(hots[[t]][[i.t]]$y01, bw=0, xs=xs.plot, ys=ys.plot, 
							bin=TRUE, cex=cex, t=t.grid[t], col=i)	
				legend("topright", paste("Achtung!", n.hots, "Hotspots"))				
			}
		}
	}	
	dev.off()
	
}


####################################################################
# Die Funktion track.eval() soll das Tracking evaluieren.
#		Für jedes Objekt wird bestimmt, ob es doppelte Hotspots oder
#		doppelte Objekte gibt. Zusätzlich wird die Anzahl Hotspots für jedes
#		Objekt bestimmt.
#
# Übergeben werden muss:
#	- track.ddf:		Ergebnis von track.func()
#
# Ausgegeben wird:
#	- obj.unique:		Matrix, die für jedes Objekt einen Indikator enthält,
#								ob es funktioniert hat (s. oben), die Anzahl 
#								Hotspots zu dem Objekt und die Nummer des 
#								Problems, falls es ein Problem gab
#	- print-Meldung:	Wie viele sind unproblematisch?
#
####################################################################

track.eval <- function(track.ddf){

	# data.frame, der für jedes Objekt die Info enthält, ob es "unique" ist
	#	zusätzlich die Anzahl der Hotspots zu dem Objekt
	#	zusätzlich das aufgetrene Problem (mehrere Objekte für ein Hotspot
	#		oder mehrere Hotspots für ein Objekt)
	anz.obj <- track.ddf$anz.obj
	objects <- 1:anz.obj
	obj.unique <- data.frame(object=objects,unique=NA,n.hots=0,problem=0)
	
	# tracker aus track.ddf
	tracker <- track.ddf$tracker
	
	# matcher aus track.ddf
	matcher <- track.ddf$matcher
	
	# object.ids initialisieren: Liste, die alle ids und die Zeitpunkte der
	#	ids enthält
	object.ids <- vector("list",length=anz.obj)
	
	# Objekte durchgehen
	for(i in objects){
		
		# ID's (Hotspots) die diesem Objekt zugeordnet sind
		# ids <- which(is.element(tracker,i))
		
		ids <- NULL
		
		# n.objects enthält die Anzahl Objekte zu diesem Hotspot
		# Achtung, ist was anderes als Anzahl Hotpots zu einem Objekt und
		#	Zeitpunkt
		
		n.objects <- rep(0,length(tracker))
		
		for(k in 1:length(tracker)){
			if(length(intersect(tracker[[k]],i))>0){
				ids <- union(ids, k)	
				n.objects[k] <- length(tracker[[k]])
			}
		}
		
		# Wenn Länge eines tracker[[ids]] größer als 1, also einem Hotspot
		#	außer dem aktuellen Objekt noch ein anderes zugeordnet ist
		#	=> Problematisch
		#	Sonst: Weitersehen
		if(max(n.objects)>1){
		
			# unique = FALSE
			obj.unique[i,2] <- FALSE
			
			# problem = mehrere Objekte
			obj.unique[i,4] <- 1
			
		}else{
			
			# Zeitpunkte der Hotspots bestimmen
			t.vec <- matcher[ids,2]
			
			# Falls eine Zeitpunkt öfter vorkommt, also mehrere Hotspots
			#	eines Zeitpunktes einem Objekt zugeordnet werden
			#	=> Problematisch
			# Sonst: Gut
			if(max(table(t.vec))>1){
				
				# unique = FALSE
				obj.unique[i,2] <- FALSE
				
				# problem = mehrere Hotspots
				obj.unique[i,4] <- 2
			}else{
				
				# Alles in Ordnung, scheint unproblematisch zu sein
				# unique = TRUE
				obj.unique[i,2] <- TRUE
				
			}			
		}
		
		# Anzahl der Hotspots
		obj.unique[i,3] <- length(ids)
		
		# IDs abspeichern und Zeitpunkte dazu
		object.ids[[i]] <- matrix(matcher[ids,1:2], ncol=2)
		
		colnames(object.ids[[i]]) <- c("id","t")
		
	}
	
	n.unique <- apply(obj.unique,MARGIN=2,FUN=sum)[2]
		
	cat(paste("Von insgesamt",anz.obj,"Objekten scheinen",n.unique,"unproblematisch zu sein \n"))
# 	cat("\n Mögliche Probleme: \n 1: Mehrere Objekte zu einem Hotspot \n 2: Mehrere Hotspots zu einem Objekt \n")
# 	cat("Hinweis: in 1 kann 2 enthalten sein \n")
	
	out <- list(track.ev=obj.unique,
					object.ids=object.ids)
	
	return(out)
}
	
	
####################################################################
# Die Funktion track.split() stellt den 2. (3.?) Teil des Trackings dar.
#		Hier sollen Tracks getrennt werden, wenn zu einem Zeitpunkt einem
#		Objekt mehr als ein Hotspot zugeordnet wurde (und die Tracks
#		gelöscht werden, die zu kurz sind)
#
# Übergeben werden muss:
#	- track.ddf:		Ergebnis von track.func()
#	- min.length:		Minimale Länge eines Tracks. Default=1
#
# Ausgegeben wird:
#	- track.ddf.up:	Update des Inputs
#
####################################################################

track.split <- function(track.ddf, min.length=1){	
	
	# Anzahl Objekte
	anz.obj <- track.ddf$anz.obj
	objects <- 1:anz.obj
	
	# tracker aus track.ddf
	tracker <- track.ddf$tracker
	
	# matcher aus track.ddf
	matcher <- track.ddf$matcher
	
	
	# Objekte durchgehen
	for(i in objects){
	
		# ID's (Hotspots) die diesem Objekt zugeordnet sind
		
		ids <- NULL
		
		# n.objects enthält die Anzahl Objekte zu diesem Hotspot
		# Achtung, ist was anderes als Anzahl Hotpots zu einem Objekt und
		#	Zeitpunkt
		
		n.objects <- rep(0,length(tracker))
		
		for(k in 1:length(tracker)){
			if(length(intersect(tracker[[k]],i))>0){
				ids <- union(ids, k)	
				n.objects[k] <- length(tracker[[k]])
			}
		}
	
		# Nur weiter, wenn das Objekt immer das einzige auf den Hotspots ist
		if(sum(n.objects>1)==0){
		
			# Zeitpunkte der Hotspots bestimmen
			t.vec <- matcher[ids,2]
			
			# Falls ein Zeitpunkt öfter vorkommt, also mehrere Hotspots
			#	eines Zeitpunktes einem Objekt zugeordnet werden
			#	=> Lösen
			# Sonst: direkt nächstes Objekt
			# Auch weiter, wenn es mal 3 oder mehr sind, dafür habe ich 
			#	noch keine Lösung implementiert
			if(max(table(t.vec))==2){
				
				# Erster Zeitpunkt wo es zwei Hotspots sind
				t.start <- min(which(table(t.vec)==2))
			
				# ACHTUNG: Wenn es zwischendurch nur einer ist, was dann?
				
				# Wenn vorher nur ein Zeitpunkt ist: Konstante Prognose
				if(t.start==2){
				
					# ID des vorherigen Zeitpunkts und die beiden des aktuellen
					id.t <- ids[(t.start-1):(t.start+1)]
					
					# Der mit dem größeren Abstand
					id.max <- id.t[which.max(as.matrix(dist(matcher[id.t,4:5]))[1,])]
					
					# bekommt eine neue Nummer
					tracker[[id.max]] <- anz.obj + 1
					
					# Nummer merken
					obj.new <- anz.obj + 1 
					
					# Anzahl Objekte um 1 raufsetzen
					anz.obj <- anz.obj + 1
					
				}
	
				# Sonst: Lineare Prognose
				if(t.start>2){
					
					# ID's der vorherigen beiden
					id.t_1 <- ids[(t.start-2):(t.start-1)]
					
					# ID's der aktuellen
					id.t <- ids[(t.start):(t.start+1)]
					
					# Differenz von t-2 auf t-1
					diff.t_1 <- diff(matcher[id.t_1,4:5])
					
					# Prognose zum Zeitpunkt t
					prog.t <- diff.t_1 + matcher[id.t_1[2],4:5]
					
					# Abstand der beiden Schwerpunkte zur Prognose
					abst <- as.matrix(dist(rbind(prog.t,matcher[id.t,4:5])))[1,2:3]
					
					# ID mit größerem Abstand
					id.max <- id.t[which.max(abst)]
					
					# bekommt eine neue Nummer
					tracker[[id.max]] <- anz.obj+1
					
					# Nummer merken
					obj.new <- anz.obj+1
					
					# Anzahl Objekte um 1 raufsetzen
					anz.obj <- anz.obj+1
					
				}
				
				# Jetzt alle anderen Zeitpunkt durchgehen und immer lineare
				#	Prognose der früheren Objektnummer machen, der Hotspot
				#	der weiter weg ist, bekommt die andere Nummer
				
				# Erster Zeitpunkt, wo es zwei Hotspots sind
				t.min <- t.vec[t.start]
				
				# Letzter Zeitpunkt, wo es zwei Hotspots sind
				t.max <- unique(t.vec)[max(which(table(t.vec)==2))]
				
				# Nur sinnvoll, wenn mehr als ein Zeitpunkt mit 2 Hotspots
				if(t.max>t.min){
					
					for(t in (t.min+1):t.max){
						
						# ID's (Hotspots) die dem alten Objekt zugeordnet sind
						
						ids <- NULL
						
						for(k in 1:length(tracker)){
							if(length(intersect(tracker[[k]],i))>0){
								ids <- union(ids, k)	
							}
						}
						
						# ID's der vorherigen beiden Zeitpunkte
						id.t_1 <- matcher[ids,][which(matcher[ids,2] %in% c(t-2, t-1)),1]
						
						# ID's der aktuellen
						id.t <- matcher[ids,][which(matcher[ids,2]==t),1]
						
						# Nur weiter, falls zum Zeitpunkt t 2 Hotspots zum Objekt 
						#	i gehören.
						if(length(id.t)==2){
							
							# Differenz von t-2 auf t-1
							diff.t_1 <- diff(matcher[id.t_1,4:5])
							
							# Prognose zum Zeitpunkt t
							prog.t <- diff.t_1 + matcher[id.t_1[2],4:5]
							
							# Abstand der beiden Schwerpunkte zur Prognose
							abst <- as.matrix(dist(rbind(prog.t,matcher[id.t,4:5])))[1,2:3]
							
							# ID mit größerem Abstand...
							id.max <- id.t[which.max(abst)]
							
							# ...bekommt eine neue Nummer
							tracker[[id.max]] <- obj.new
						}
						
# 						# Falls es sogar mehr sind: Warnung ausgeben
						if(length(id.t)>2){
							cat(paste("Achtung:",length(id.t),"Hotspots bei Objekt",i,"zum Zeitpunkt",t))
						}
						# => Brauche ich nicht, da ich eh nur weiter mache, wenn
						#		es maximal 2 sind...
						
					}
				}			
			}
			if(max(table(t.vec))>2){
				cat("Mehr als zwei Hotspots bei einem Objekt \n")
			}
		}				
	}
	
# 	# Lösche alle Tracks, die kürzer sind als min.length
# 	
# 	objects <- 1:anz.obj
# 	obj.length <- data.frame(object=objects,length=0)
# 	
# 	for(i in objects){
# 		
# 		ids <- NULL
# 		
# 		for(k in 1:length(tracker)){
# 			if(length(intersect(tracker[[k]],i))>0){
# 				ids <- union(ids, k)	
# 			}
# 		}
# 		
# 		# Länge des Tracks bestimmen
# 		obj.length[i,2] <- length(ids)
# 		
# 		# Wenn kürzer als min.length
# 		if(obj.length[i,2]<min.length){
# 			
# 			# Setze alle Hotspots mit diesen ids auf Objekt=0 
# 			#	(und lösche die Zeile aus dem matcher)
# 			
# 			for(q in ids){
# 				tracker[[q]] <- 0
# # 				matcher <- matcher[-q,]
# 			}
# 			
# 			# Verringere anz.obj um 1
# 			anz.obj <- anz.obj-1
# 		}
# 	}

	
	# Zurückgegeben wird Liste wie vorher + obj.length
	
	out <- list(tracker=tracker, 
					anz.obj=anz.obj,
					matcher=matcher,
					t.plot=track.ddf$t.plot
# 					,obj.length=obj.length
					)
	
	return(out)
}
	

####################################################################
# Die Funktion track.split.over() stellt den 3.(2.) Teil des Trackings dar.
#		Hier sollen weitere Probleme abgefangen werden, die mit sich 
#		überlappenden oder zu nah aneinanderliegenden Hotspots entstehen.
#
# Übergeben werden muss:
#	- track.ddf:		Ergebnis von track.func()
#
# Ausgegeben wird:
#	- track.ddf.up:	Update des Inputs
#
####################################################################

track.split.over <- function(track.ddf){
	
	# Alles relevante aus track.ddf rausziehen
	# Anzahl Objekte
	anz.obj <- track.ddf$anz.obj
# 	objects <- 1:anz.obj
	
	# tracker aus track.ddf
	tracker <- track.ddf$tracker
	
	# matcher aus track.ddf
	matcher <- track.ddf$matcher
	
	# Schleife über alle Objekte
	#	While-Schleife, da Objekte dazukommen können, und ich die dann auch
	#	gleich analysieren will. Solange i noch kleiner gleich der Anzahl
	#	Objekte ist (die sich in der Schleife ändern kann) wird das ganze
	#	Programm durchgezogen
	#
	
	i <- 1
	
	while(i <= anz.obj){
		
		# n.objects bestimmen, also den Vektor, der für alle Hotspots, bei 
		#	denen i mitspielt die Anzahl der Objekte enthält, die auf 
		#	den Hotspot zeigen
		
		n.objects <- rep(0,length(tracker))
		
		# zusätzlich die ID's derjenigen Hotspots, bei denen i vorkommt
		ids <- NULL
		
		
		for(k in 1:length(tracker)){
			if(length(intersect(tracker[[k]],i))>0){
				ids <- union(ids, k)	
				n.objects[k] <- length(tracker[[k]])
			}
		}
		
		# Falls i Teil einer Doppelbelegung ist (nicht dreifach!) wird 
		#	weitergemacht, sonst abgebrochen		
		if(max(n.objects)==2){
		
# 			# Zeitpunkte der Hotspots bestimmen
			t.vec <- matcher[ids,2]
			
			# Bestimme 1. Zeitpunkt der Doppelbelegung und letzten Zeitpunkt
			#	überhaupt in dem i vorkommt
			
			# 1. ID mit Doppelbelegung
			id.min <- min(which(n.objects==2))
			
			t.min <- matcher[id.min,2]			
			t.max <- max(t.vec)
			
			# Flags für Überschneidung und neue Nummer auf 0 setzen
			flag.overlap <- FALSE
			flag.newnumber <- FALSE
			
			# Andere Objektnummer
			j <- setdiff(tracker[[id.min]],i)
			
			# Statt einer for(t in t.min:t.max)-Schleife hier wieder
			#	eine while-Schleife, damit abgebrochen werden kann,
			#	wenn eine neue Nummer vergeben wurde
			t <- t.min
			while(is.element(t,t.min:t.max)){
			
				# ID's zum Zeitpunkt t
				id.t <- matcher[ids,][which(matcher[ids,2]==t),1]
				
				# Wenn das mehr als zwei sind ist irgendwas komisch
				if(length(id.t)>2){
					cat("Mehr als zwei ID's")
				}else{
					
					
					# Wenn nur ein Hotspot und  darin Doppelbelegung
					#	=> Hotspot-Objekt auf 0 setzen
					#	=> Flag für Überschneidung = TRUE
					if(length(id.t)==1){
						if(length(tracker[[id.t]])==2){
							# Die zweite Bedingung könnte ich mir wahrscheinlich
							#	sparen, eigentlich müssten es ab t.min ja immer
							#	mindestens Doppelbelegungen sein.
						
							flag.overlap <- TRUE
							tracker[[id.t]] <- 0
						}
					}
					
					# Wenn 2 Hotspots
					if(length(id.t)==2){
						
						# Anzahl Objekte zu den beiden Hotspots			
						l.1 <- length(tracker[[id.t[1]]])
						l.2 <- length(tracker[[id.t[2]]])
						
						# Wenn eine ID eine Doppelbelegung und 
						#	die andere ID eine Einzelbelegung hat
						if((l.1==2 & l.2==1)| (l.1==1 & l.2==2)){
							
							# Die ID mit der Einzelbelegung und Doppelbelegung
							if(l.1==1){
								id.einzeln <- id.t[1]
								id.doppelt <- id.t[2]
							}else{
								id.einzeln <- id.t[2]
								id.doppelt <- id.t[1]
							}
							
							# Das Objekt der Einzelbelegung
							obj.einzeln <- tracker[[id.einzeln]]
							
							# Die andere ID bekommt das andere Objekt zugewiesen
							if(obj.einzeln==i){
								tracker[[id.doppelt]] <- j
							}else{
								tracker[[id.doppelt]] <- i
							}							
						}
						
						# Wenn beide ID's Doppelbelegung
						if(l.1==2 & l.2==2){
							
							# Wenn es noch keine Überschneidung gab
							if(!flag.overlap){
								
								# Wenn nur ein Zeitpunkt vorher ist:
								#	Konstante Prognose
								if(min(which(t.vec==t.min))==2){
									
									# ID des letzten einzelnen Auftretens von i
									id.t_1 <- max(which(is.element(tracker,i)))
									
									# Der mit dem größeren Abstand
									id.max <- id.t[which.max(as.matrix(dist(matcher[c(id.t_1,id.t),4:5]))[1,2:3])]
																	
									# bekommt die andere Nummer
									tracker[[id.max]] <- j
									
									# Der mit dem kleineren Abstand 
									id.mini <- id.t[which.min(as.matrix(dist(matcher[c(id.t_1,id.t),4:5]))[1,2:3])]
									
									# bekommt i
									tracker[[id.mini]] <- i									
									
								# Sonst: Lineare Prognose
								}else{
									
									# ID des letzten einzelnen Auftretens von i
									id.t_1 <- which(is.element(tracker,i))[-(1:(length(which(is.element(tracker,i)))-2))]
									
									# Differenz von t-2 auf t-1
									diff.t_1 <- diff(matcher[id.t_1,4:5])
									
									# Prognose zum Zeitpunkt t
									prog.t <- diff.t_1 + matcher[id.t_1[2],4:5]
									
									# Abstand der beiden Schwerpunkte zur Prognose
									abst <- as.matrix(dist(rbind(prog.t,matcher[id.t,4:5])))[1,2:3]
									
									# ID mit größerem Abstand
									id.max <- id.t[which.max(abst)]
									
									# bekommt Nummer j
									tracker[[id.max]] <- j
									
									# ID mit kleinerem Abstand
									id.mini <- id.t[which.min(abst)]
									
									# bekommt Nummer i
									tracker[[id.mini]] <- i
								}								
							}
							
							# Wenn es bereits eine Überschneidung gab
							if(flag.overlap){
								
								# Wenn noch keine neue Nummern vergeben wurden:
								#	dies tun, sonst:
								# Prognose konstant oder linear
								#
								# Update, neue Idee:
								#	Ab hier kümmere ich mich nicht um die weitere
								#	Verfolgung sondern setze einfach alle danach 
								#	auf (k,l). Durch die while-Schleife komme
								#	ich bei diesen Objekten ja noch vorbei und 
								#	springe dann in den Fall 
								#	"Noch keine Überschneidung, aber 2 Hotspots
								#	mit Doppelbelegung"
								# => Mal probieren ob das geht, da würde ich mir
								#	einiges an Code sparen
								if(!flag.newnumber){
									
									# Flag auf TRUE setzen
									flag.newnumber <- TRUE
									
									# Neue Nummern vergeben									
									k <- anz.obj+1
									l <- anz.obj+2
									
									tracker[[id.t[1]]] <- k
									tracker[[id.t[2]]] <- l
									
									anz.obj <- l
									
									# Alle Objekte der ID's danach auf (k,l) setzen
									ids.follow <- ids[-(1:which(ids==id.t[2]))]
									
									for(p in ids.follow){
										tracker[[p]] <- c(k,l)
									}
									
									# und Schleife abbrechen
									t <- t.max
									
									
								} # Sonst: Gar nichts
							}
						}
						
						# Wenn beide ID's Einzelbelegung: Darf nicht vorkommen
						if(l.1==1 & l.2==1){
							cat("Zwei Einzelbelegungen!")
						}
					}
				}	
				
				# t um eins raufsetzen
				# Wenn neue Nummern vergeben wurden, wurde t auf t.max gesetzt
				#	=> t ist dann t.max + 1 was zum Abbruch der Schleife führt
				t <- t+1
			}			
		}		
		
		# Falls i Teil einer Dreifach- oder Mehrbelegung ist: Meldung 
		if(max(n.objects)>2){
			cat(paste("Objekt",i,"ist Teil einer",max(n.objects),"er- Belegung \n"))
		} 
		
		# Falls i immer nur Einfachbelegungen hat geht es einfach weiter 
		#	mit i+1
		
		
		# Setze i um eins rauf
		i <- i+1
	}
	
	# Zurückgegeben wird Liste wie üblich
	
	out <- list(tracker=tracker, 
					anz.obj=anz.obj,
					matcher=matcher,
					t.plot=track.ddf$t.plot)
	
	return(out)
	
}


####################################################################
# Die Funktion track.eval.total() soll den Tracking-Prozess für den
#	ganzen Film evaluieren, also eine Tabelle erstellen wie viele 
#	problematisch sind oder nicht und wie viele Tracks welcher Länge es 
#	gibt.
#
# Übergeben werden muss:
#	- track.evals:		track.evals aus dem Preprocessing
#	- min.length:		Minimale Länge interessierender Tracks
#	- save.plot:		TRUE: Table abspeichern
#	- pfad.mult:		Pfad für den Plot
#	- fisch:				Fischart (Für den Titel des Plots)
#	- a.2:				Für Dateinamen
#
# Ausgegeben wird:
#	- plot
#	- track.summary:	Übersichtstabelle, wie viele (un)problematisch sind
#
####################################################################

track.eval.total <- function(track.evals,
									  min.length,
									  save.plot=FALSE,
									  pfad.mult,
									  fisch="Fischart",
									  a.2){
	
	tracks <- 0
	
	for(i in 1:length(track.evals)){	
		tracks <- tracks + length(which(track.evals[[i]][,3]>=min.length))
	}
	
	# Anzahl Tracks die länger sind als min.length
	#tracks
	
	
	###########
	# Anzahl problematische / unproblematische (von denen, die größer sind als
	#		min.length)
	###########
	
	no.prob <- 0
	
	for(i in 1:length(track.evals)){	
		no.prob <- no.prob + 
			sum(track.evals[[i]][which(track.evals[[i]][,3]>=min.length),2])
	}
	
	#no.prob
	#tracks
	
	###########
	# Histogramm der Pfadlängen
	###########
	
	track.length <- NULL
	
	for(i in 1:length(track.evals)){	
		track.length <- c(track.length,track.evals[[i]][,3])
	}
	
	#track.length
	
	if(!is.null(track.length)){
# 		if(!save.plot){
# 			hist(track.length)
# 			plot(table(track.length), ylim=c(0,50))
# 		}
		
		# table.forelle <- table(track.length)
		
		if(save.plot){
			filename <- paste(pfad.mult,"Track_Table_a2_",a.2,".jpg",sep="")
			jpeg(filename,
				  width=13,
				  height=8, 
				  units="in", 
				  res=200, 
				  quality=100)
			plot(table(track.length), main=fisch)
			dev.off()
		}
	}
	
	###########
	# Tracks mit Problem 1: Mehrere Objekte zu einem Hotspot 
	#	Achtung: Hier könnte auch 2 aufgetreten sein!
	###########
	
	prob.1 <- 0
	
	for(i in 1:length(track.evals)){	
		prob.1 <- prob.1 + 
			length(which(track.evals[[i]][which(track.evals[[i]][,3]>=min.length),4]==1))
	}
	
#	prob.1
	
	
	###########
	# Tracks mit Problem 2: Mehrere Hotspots zu einem Objekt 
	###########
	
	prob.2 <- 0
	
	for(i in 1:length(track.evals)){	
		prob.2 <- prob.2 + 
			length(which(track.evals[[i]][which(track.evals[[i]][,3]>=min.length),4]==2))
	}
	
#	prob.2
	
	tracks.summary <- data.frame(tracks=tracks,
										  unproblematisch=no.prob,
										  problem.1=prob.1,
										  problem.2=prob.2)
#	tracks.summary

	
	return(tracks.summary)
}
	
	