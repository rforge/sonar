################################################################################
# Fischprojekt
# Funktionen zur Bestimmung der Merkmale eines Fischs
# Datum der ersten Version: 19.04.12
# Letzte Änderung: 22.10.12
# Autor: Ludwig Bothmann
################################################################################

####################################################################
# Die Funktion merkmale.func() soll die Merkmale aller Objekte 
#	bestimmen, für die Hotspots gefunden wurden. Das Tracking muss
#	bereits erfolgt sein.
#
#	Für die Vektoren wird eine neue Idee im Vergleich zu früher verfolgt:
#		Berechnung basierend auf Polarkoordinaten, nicht basierend
#		auf einer Drehung der Index-Matrix
#
# Übergeben werden muss:
#	- hotspots:			Ergebnis-Liste aus find.mult.hotspots()
#	- track.ddf:		Ergebnis-Liste aus track.func() oder track.split()
#								oder track.split.over()
#	- track.eval.out:	Ergebnis-Liste aus track.eval()
#	- n.angle:			Anzahl der Vektoren pro Viertel
#  - cart.coord:		3D-Array mit kartesischen Koordinaten der Pixel
#	- t.plot:			Die Indizes der interessanten Zeitpunkte
#								=> Brauche ich das noch?
# - only.watch:		Trapezmerkmale durch watch-hands-Merkmale ersetzen
#
# Ausgegeben wird Liste mit:
#	- merkmale:			Liste mit Matrix der Merkmale für jedes Objekt
#	- v:					Vektor mit Geschwindigkeiten der Objekte
#	- richtung:			Vektor: Kommt das Objekt von rechts => TRUE
#
####################################################################


merkmale.func <- function(hotspots, 
													track.ddf,
													track.eval.out,
													n.angle,
													cart.coord,
													t.plot,
													only.watch){
	
	# Alle relevanten Objekte aus dem Input rausziehen
	
	hots <- hotspots$hots
	Y.01 <- hotspots$Y.01
	
	tracker <- track.ddf$tracker
	anz.obj <- track.ddf$anz.obj
	matcher <- track.ddf$matcher
	
	object.ids <- track.eval.out$object.ids
	track.ev <- track.eval.out$track.ev
	
	xs <- cart.coord[,,1]
	ys <- cart.coord[,,2]
	
	##################################
	# Fläche eines Pixels bestimmen	
	##################################
	
	fl.pix <- vector(length=dim(Y.01)[1])
	# Vektor, da die Fläche nur vom Abstand von der Linse, nicht aber
	#	vom Winkel abhängt.
	
	for(i in 2:(dim(Y.01)[1]-1)){
		b <- abstand(c(i,1), c(i,3),cart.coord)		
		l <- abstand(c(i-1,2), c(i+1,2),cart.coord)		
		fl.pix[i] <- (b/2) * (l/2)
	}
	
	# Ersten und letzten separat, weil da die Nachbarn fehlen
	fl.pix[1] <- fl.pix[2]
	fl.pix[dim(Y.01)[1]] <- fl.pix[(dim(Y.01)[1]-1)]
	
	################################################################################
	# Merkmale bestimmen
	################################################################################
	
	# Outputliste initialisieren
	merkmale.list <- vector("list", length=anz.obj)
	
	# Vektoren für Geschwindigkeit und Richtung initialisieren
	v.vec <- vector(length=anz.obj) 
	richt.vec <- vector(length=anz.obj)
	
	# 26.11.12: Mean und sd von Richtung und Geschwindigkeit
	dir.vel.list <- vector("list", length=anz.obj)
	
	# Schleife über alle Objekte
	for(i in 1:anz.obj){

		# Nur weiter, wenn es nicht schon als problematisch erkannt wurde,
		#	Sonst bleibt NULL
		if(track.ev[i,2]==TRUE){
		
			# Merkmal-Matrix initialisieren
			merkmale <- matrix(nrow=(11+n.angle*8+1), 
									 ncol=nrow(object.ids[[i]]))
			
			# Merkmal-Vektor initialisieren
			merk <- vector(length=nrow(merkmale))
			
			# Zeilennamen festlegen
			rownames(merkmale) <- c("Breite (Dim.x)", "Laenge (Dim.y)", "Flaeche", 
											"Seitenverhaeltnis", "Volumen", "Mitte Dim.x", 
											"Mitte Dim.y", 
											"Marginaler Schwerpunkt Dim.x Ind", 
											"Marginaler Schwerpunkt Dim.y Ind", 
											"Marginaler Schwerpunkt Dim.x kart", 
											"Marginaler Schwerpunkt Dim.y kart",
											paste(90/n.angle * 0:((4*n.angle)-1), "Grad", sep=" "),
											paste(90/n.angle * 0:((4*n.angle)-1), "Grad wahr", sep=" "),
											"Hauptachse in Grad zur Horiz.")
			
			##################################################################	
			# Hotspots durchgehen und Merkmale - ohne Vektoren - bestimmen
			##################################################################
			
			# Schleife über Hotspots zu Objekt i
			ids <- object.ids[[i]][,1]
			
			for(k in 1:length(ids)){
				
				# ID des Hotspots
				id <- ids[k]
				
				# Zeitpunkt des Hotspots
				t <- matcher[id,2]
				
				# Gruppe innerhalb des Zeitpunkts
				j <- matcher[id,3]
				
				# 0/1 - Bild des Hotspots - Ausschnitt
				y <- hots[[t]][[j]]$y01
				
				# Nur weiter, falls y eine Matrix ist also irgendwas drin 
				#	steht
				if(is.matrix(y)){
					
					# 0/1 - Bild des Hotspots - Ganzes Bild
					y01.all <- hots[[t]][[j]]$y01.all
					
					# Dimensionen des Rechtsecks
					dims <- hots[[t]][[j]]$dims
					
					# 1. Breite des Rechtecks			
					p1 <- c(dims[1,1],dims[1,2])
					p2 <- c(dims[1,1],dims[2,2])
					p3 <- c(dims[2,1],dims[1,2])					
					breite <- abstand(p1, p2,cart.coord)
					merk[1] <- breite
					
					# 2. Laenge des Rechtecks	
					laenge <- abstand(p1, p3,cart.coord)
					merk[2] <- laenge
					
					# 3. Fläche des Rechtecks
					merk[3] <- merk[1] * merk[2]
					
					# 4. Seitenverhältnis
					merk[4] <- merk[2] / merk[1]
					
					# 5. Anzahl aktivierter Pixel / Volumen				
					# Einfach die Summe reicht nicht: 
					#	Fläche eines Pixels nutzen
					anz <- sum(y01.all)
					vol <- sum(y01.all * fl.pix)				
					merk[5] <- vol
					
					# 6. & 7. Koordinaten des Mittelpunkts des Rechtecks (kartesich)
					merk[6] <- cart(c(mean(dims[,1]), mean(dims[,2])),cart.coord)[1]
					merk[7] <- cart(c(mean(dims[,1]), mean(dims[,2])),cart.coord)[2]
					
					# 8. & 9. Marginale Schwerpunkte absolut als Indizes
					indizes <- which(y01.all==1, arr.ind=TRUE)
					merk[9:8] <- apply(X=indizes,MARGIN=2,FUN=mean)

					# 10. & 11. Marginale Schwerpunkte absolut kartesisch								
					
					# x- und y-Koordinaten
					x.s <- xs[which(y01.all==1)]
					y.s <- ys[which(y01.all==1)]
					
					# Schwerpunkt
					s.x <- mean(x.s)
					s.y <- mean(y.s)
					merk[10:11] <- c(s.x, s.y)
					
					###########################
					# ENDE: Alles abspeichern
					###########################
					
					if(!only.watch){
						merkmale[,k] <- merk
					}else{
						merkmale[5:11,k] <- merk[5:11]
					}				
					
				}
			}
			
			##################################################################
			# Geschwindigkeit und Schwimmrichtung (von rechts oder links)
			##################################################################
			
			# Ersten und letzten Zeitpunkt bestimmen für den Merkmale 
			#	berechnet wurden
			t.min <- min(which(!is.na(merkmale[5,])))
			t.max <- max(which(!is.na(merkmale[5,])))
			
			v <- sqrt((merkmale[10,t.min] - merkmale[10,t.max])^2 
						 + (merkmale[11,t.min] - merkmale[11,t.max])^2) / (t.max-t.min+1)
			names(v) <- "Geschwindigkeit (Meter /Frame)"
			
			richtung <- (merkmale[10,t.min]>merkmale[10,t.max])
			names(richtung) <- "von rechts"
			
			##################################################################
			# 26.11.12: Ich möchte die Bewegung noch besser abbilden, 
			#		deshalb wird jetzt der durchschnittliche Bewegungswinkel plus
			#		Standardabweichung sowie die durchschnittliche zurückgelegte
			#		Strecke plus Standardabweichung berechnet
			##################################################################
			
			# Nur, wenn mindestens 3 IDs, sonst kann die Standardabweichung
			#		nicht berechnet werden. (Da sonst nur 1 Differenz)
			
			if(length(ids)>2){
				
				# Differenzen kartesisch
				diff.cart <- diff(t(merkmale[(10:11),]))
				
				# In Polarkoordinaten umrechnen
# 			diff.pol.alt <- cart2pol(x=diff.cart[,1], y=diff.cart[,2])
				diff.pol <- cart2pol90(x=diff.cart[,1], y=diff.cart[,2])
				
				# Mean und sd des Winkels
				alpha.mean <- mean(diff.pol[,1])
				alpha.sd <- sd(diff.pol[,1])
				
				# ACHTUNG: Wenn ein Fisch mal knapp über 0, mal knapp unter 0 ist,
				#		also fast 2*pi, ergibt sich im Mittel pi, was natürlich total 
				#		daneben ist... Kann man das lösen?
				
				# Ich könnte es umgehen, indem ich sage, dass 0 (und 2 pi) oben
				#		sein soll, dann habe ich das Problem nur quer zum Strom, 
				#		sonst könnte es aber sogar problematisch sein zu unterscheiden,
				#		ob das Objekt mit oder gegen den Strom schwimmt.
				# => cart2pol90 schreiben, die das macht.
				
				# UPDATE: Habe es implementiert, sollte funktionieren
				
				# Mean und sd des Radius
				rad.mean <- mean(diff.pol[,2])
				rad.sd <- sd(diff.pol[,2])
				
				# Hier dürfte es das Problem nicht geben
				
				# Zusammenfassen zu einem Vektor
				dir.vel.vec <- c(alpha.mean, alpha.sd, rad.mean, rad.sd)
				
			}else{
				
				# Sonst leerer Vektor
				dir.vel.vec <- rep(0,4)
				
			}
			
			names(dir.vel.vec) <- c("alpha.mean", "alpha.sd", 
															"rad.mean", "rad.sd")
				
			##########################################################################
			# Vektoren bestimmen
			##########################################################################	
			
			# Die Zeitpunkte bestimmen, in denen das Objekt da ist
			k.not.na <- which(!is.na(merkmale[5,]))
			
			# Schleife über alle Hotspots zu Objekt i
			
			for(k in k.not.na){

				# ID des Hotspots
				id <- ids[k]
				
				# Zeitpunkt des Hotspots
				t <- matcher[id,2]
				
				# Gruppe innerhalb des Zeitpunkts
				j <- matcher[id,3]
				
				#####################################
				# Vektorlängen bestimmen
				#####################################		
				
				# 0/1 - Bild des Hotspots - Ganzes Bild
				y01.all <- hots[[t]][[j]]$y01.all
								
				# x- und y-Koordinaten
				x.s <- xs[which(y01.all==1)]
				y.s <- ys[which(y01.all==1)]
				xy.coord <- cbind(x.s, y.s)
				
				# Schwerpunkt
				s.x <- merkmale[10,k]
				s.y <- merkmale[11,k]
				anz <- sum(y01.all)
				s.xy <- cbind(rep(s.x,anz),rep(s.y,anz))
				
				# Zentriertes Bild
				y.zent <- xy.coord - s.xy
				
				
				#####################################
				# Winkel zur Horizontalen bestimmen
				#####################################
				
				# Steigung, der Regressionskoeffizient
				beta <- lm(y.zent[,2]~y.zent[,1])$coef[2]
				
				# Steigung in Bogenmaß
				alpha.start <- atan(beta)
				
				# 12. etc. Uhrzeigerlängen
				
				# Vektorlängen berechnen
				clock <- clock.vecs.func(y.zent, 
												n.angle, 
												alpha.start=alpha.start, 
												von.rechts=richtung)

				# Falls die Trapezvariablen nicht benutzt werden sollen,
				#		brauche ich jetzt noch Länge, Breite, Seitenverhältnis
				#		und Seitenprodukt
				if(only.watch){
					
					if(is.na(clock$vecs["0"]) | is.na(clock$vecs["180"]) | is.na(clock$vecs["90"]) | is.na(clock$vecs["270"])){
						stop("Fehler in merkmale.func(): Laenge und Breite kann bei diesem Wert von n.angle nicht bestimmt werden")
					}
					
					merkmale[1,k] <- clock$vecs["90"] + clock$vecs["270"]
					merkmale[2,k] <- clock$vecs["0"] + clock$vecs["180"]
					
					if(merkmale[1,k]>0 & merkmale[2,k]>0){
						merkmale[3,k] <- merkmale[2,k] * merkmale[1,k]
						merkmale[4,k] <- merkmale[2,k] / merkmale[1,k]
					}else{
						merkmale[3,k] <- NA
						merkmale[4,k] <- NA
# 						cat("Achtung mit Seitenverhaeltnis und Flaeche (merkmale.func()) \n")
					}
					# Ist das ok so? Das Problem ist, dass manchmal 0 auftreten kann,
					#		dann ist das Verhältnis unendlich
				}
				
				# Normieren, so dass X0.Grad = 1, X0.Grad bleibt aber auf altem Wert
				#	als Referenz.
				# Nur, wenn X0.Grad nicht 0 ist
				if(clock$vecs[1]>0){
					clock$vecs[-1] <- clock$vecs[-1] / clock$vecs[1]
				}
				
				# In merkmale abspeichern
				merkmale[12:(11+4*n.angle),k] <- clock$vecs
				merkmale[(12+4*n.angle):(11+8*n.angle),k] <- clock$angles
				merkmale[(11+8*n.angle+1),k] <- alpha.start*180/pi		
			}
			
			merkmale.list[[i]] <- merkmale
			v.vec[i] <- v
			richt.vec[i] <- richtung
			dir.vel.list[[i]] <- dir.vel.vec		
			
		}
	}

	
	
	out <- list(merkmale.list=merkmale.list, 
							v.vec=v.vec, 
							richt.vec=richt.vec,
							dir.vel.list=dir.vel.list)
	
	return(out)
	
}


####################################################################
# Die Funktion cart() berechnet die kartesischen Koordinaten 
#	eines oder mehrerer Punkte.
#
# Übergeben werden muss:
#	- vec:		Die Indizes des Punkts (Vektor) oder der Punkte (nx2-Matrix) 
#
# Ausgegeben wird:
#	- coords:	Die kartesischen Koordinaten im gleichen Format
#
# ACHTUNG: cart.coord muss bisher global definiert sein.
# => UPDATE 18.10.12: Jetzt wird es übergeben
#
####################################################################



cart <- function(vec, cart.coord){
	
	# Nur ein Punkt
	if(is.null(dim(vec))){
		coords <- cart.coord[round(vec[1]), round(vec[2]),]
		return(coords)
	}
	
	# Mehrere Punkte
	if(!is.null(dim(vec))){
		if(dim(vec)[2]!=2) print("Matrix muss zwei Spalten haben!")
		
		# Koordinaten bestimmen
		coords <- matrix(nrow=dim(vec)[1], ncol=dim(vec)[2])
		for(i in 1:dim(vec)[1]){
			coords[i,] <- cart.coord[round(vec[i,1]), round(vec[i,2]),]
		}
		
		return(coords)
	}	
}

####################################################################
# Die Funktion abstand() berechnet den Abstand von zwei Punkten, 
#	übergeben werden die Indizes, nicht die kartesischen Koordinaten.
#
# Übergeben werden muss:
#	- p1, p2:		Die beiden Punkte
#	- cart.coord:	Kartesische Koordinaten aller Pixel
#
# Ausgegeben wird:
#	- abst:		Der Abstand
#
####################################################################


abstand <- function(p1, p2, cart.coord){
	
	punkt1 <- cart(p1, cart.coord)
	punkt2 <- cart(p2, cart.coord)	
	
	# Abstand berechnen
	abst <- sqrt(sum(apply(cbind(punkt1, punkt2), 1, diff)^2))
	
	return(abst)
	
}


####################################################################
# Die Funktion clock.vecs.func() soll in einem Bild die Länge 
#	der "Uhrzeiger" bestimmen, allerdings für die Analyse basierend auf den
#	ddf-Files.
#
# Hier wird die Idee mit den Polarkoordinaten umgesetzt
#
# Übergeben werden muss:
#	- y.zent:		Matrix mit den Koordinaten der Pixel (um 0 zentriert)
#	- n.angle:		Gewünschte Anzahl an Uhrzeigern pro Viertel
#	- alpha.start: Winkel der Hauptachse zur Horizontalen
#	- von.rechts:	TRUE: Objekt kommt von rechts, schaut also nach links
#
# Ausgegeben wird List mit:
#	- vecs:		Vektor mit den Längen der einzelnen Vektoren
#	- angles:	Wahre Winkel IN GRAD NICHT IM BOGENMAß!!!
#
####################################################################

clock.vecs.func <- function(y.zent, 
									 n.angle, 
									 alpha.start, 
									 von.rechts){
	
	# Vektor der Ergebnisse
	vecs <- vector(length=4*n.angle)
	
	# Winkeleinheit bestimmen
	alpha.unit <- 90/n.angle
	
	names(vecs) <- paste(alpha.unit * 0:((4*n.angle)-1), sep=" ")
	
	# Vektor der wahren Winkel
	angles <- vector(length=4*n.angle)
	names(angles) <- names(vecs) 

	# Winkel für die ein Uhrzeiger bestimmt werden soll (im Bogenmaß)
	ang <- alpha.unit*pi/180 * 0:((4*n.angle)-1)
	
	# Umrechnung der kartesischen Koordinaten in Polarkoordinaten
	x <- y.zent[,1]
	y <- y.zent[,2]
	pol.coord <- cart2pol(x,y)
	
	# alpha.start soll in (0,2*pi) sein
	if(alpha.start<0) {alpha.start <- alpha.start + 2*pi}
	
	# Schleife über die Winkel
	for(w in 1:length(ang)){
		
		# Winkel in dessen Nähe gesucht werden soll:
		if(von.rechts==TRUE){
			
			# Kommt von rechts
			winkel <- pi + alpha.start + ang[w]
			
			# UPDATE 17.07.13: +ang[w] statt -ang[w], da ich ja immer gegen
			#		den Uhrzeiger laufen muss. Ich sehe den Fisch ja von oben,
			#		nicht von der Seite
			# => Sollte so funktionieren
						
		}else{
			
			# Kommt von links
			winkel <- alpha.start + ang[w]
			
		}
		
		# winkel soll in [0,2*pi) sein
		if(winkel>=4*pi){ 
			winkel <- winkel - 4*pi
		}else{
			if(winkel>=2*pi){ 
				winkel <- winkel - 2*pi
			}else{
				if(winkel<0){ 
					winkel <- winkel + 2*pi
				}
			}
		}
		
		# => mehr als 4pi oder weniger als -2pi dürfte nicht vorkommen
		# => Falls doch: abbrechen
		if(winkel < 0 | winkel > 2*pi) stop("Winkel in clock.vecs.func() nicht in (0, 2*pi)")
		
		# Umgebung von 5 Grad
		region <- winkel + c(-1,1)*5*pi/180
		
		# Pixel in dieser Region bestimmen, also wenn alpha in region .
		#	Problematisch ist es, wenn die Region über [0,2*pi) rausgeht
		
		if(region[1]>=0 & region[2]<2*pi){

			# => Alles ok
			in.reg <- which(pol.coord[,1]>=region[1] & pol.coord[,1]<=region[2])
			
		}else if(region[1]<0){
 
			# untere Grenze unter 0
			# => von 0 bis obere Grenze ok, oder
			# => von untere Grenze + 2*pi bis 2*pi ok
			in.reg <- which(pol.coord[,1]>=(region[1]+2*pi) | pol.coord[,1]<=region[2])
			
		}else if(region[2]>=2*pi){
			
			# obere Grenze 2*pi oder mehr
			# => von untere Grenze bis 2*pi ok, oder
			# => von 0 bis obere Grenze - 2*pi
			in.reg <- which(pol.coord[,1]>=region[1] | pol.coord[,1]<=(region[2]-2*pi))
			
		}
		
		
		# Falls da nichts sein sollte: auf 10 Grad gehen
		if(length(in.reg)==0){
			
			# Umgebung von 10 Grad
			region <- winkel + c(-1,1)*10*pi/180
			
			# Pixel in dieser Region
			if(region[1]>=0 & region[2]<2*pi){
				in.reg <- which(pol.coord[,1]>=region[1] & pol.coord[,1]<=region[2])				
			}else if(region[1]<0){
				in.reg <- which(pol.coord[,1]>=(region[1]+2*pi) | pol.coord[,1]<=region[2])				
			}else if(region[2]>=2*pi){
				in.reg <- which(pol.coord[,1]>=region[1] | pol.coord[,1]<=(region[2]-2*pi))				
			}			
		}
		
		# Falls da immer noch nichts sein sollte: auf 15 Grad gehen
		if(length(in.reg)==0){
			
			# Umgebung von 15 Grad
			region <- winkel + c(-1,1)*15*pi/180
			
			# Pixel in dieser Region
			if(region[1]>=0 & region[2]<2*pi){
				in.reg <- which(pol.coord[,1]>=region[1] & pol.coord[,1]<=region[2])				
			}else if(region[1]<0){
				in.reg <- which(pol.coord[,1]>=(region[1]+2*pi) | pol.coord[,1]<=region[2])				
			}else if(region[2]>=2*pi){
				in.reg <- which(pol.coord[,1]>=region[1] | pol.coord[,1]<=(region[2]-2*pi))				
			}			
		}
		
		# Falls da immer noch nichts sein sollte: r auf 0 setzen
		if(length(in.reg)==0){
			
			r.max <- 0
			alpha.max <- winkel*180/pi
			
		}else{
		
			# Welcher hat den größten Radius?
			ind.max <- which.max(pol.coord[in.reg,2])
			
			# Wie groß ist der Radius?
			r.max <- pol.coord[in.reg,2][ind.max]
			
			# Welcher Winkel ist das in Wirklichkeit?
			alpha.max <- pol.coord[in.reg,1][ind.max]
		}
		
		# Abspeichern und zum nächsten Vektor
		vecs[w] <- r.max
		angles[w] <- alpha.max*180/pi
		
	}
	
	
	out <- list(vecs=vecs, angles=angles)
	
	return(out)
}



####################################################################
# Die Funktion cart2pol() rechnet kartesische Koordinaten in 
#	Polarkoordinaten um, und zwar so dass die Winkel in [0,2*pi) liegen
#
# Übergeben werden muss:
#	- x:			Kartesische x-Koordinaten		
#	- y:			Kartesische y-Koordinaten
#
# Ausgegeben wird:
#	- alpha:		Winkel
#	- r:			Radius
#
####################################################################


cart2pol <- function(x,y){
	
	r <- sqrt(x^2+y^2)
	alpha <- atan2(y,x)
	alpha[which(alpha<0)] <- alpha[which(alpha<0)] + 2*pi
	out <- cbind(alpha,r)
	return(out)
	
}

####################################################################
# Die Funktion cart2pol90() rechnet kartesische Koordinaten in 
#	Polarkoordinaten um, und zwar so dass die Winkel in [0,2*pi) liegen.
#		Bei dieser Variante ist 0 (und 2 pi) oben
#
# Übergeben werden muss:
#	- x:			Kartesische x-Koordinaten		
#	- y:			Kartesische y-Koordinaten
#
# Ausgegeben wird:
#	- alpha:		Winkel
#	- r:			Radius
#
####################################################################


cart2pol90 <- function(x,y){
	
	r <- sqrt(x^2+y^2)
	alpha <- atan2(y,x)
	alpha[which(alpha<0)] <- alpha[which(alpha<0)] + 2*pi
	
	# 90 Grad abziehen
	alpha90 <- alpha - pi/2
	
	# Wieder auf [0, 2pi) ändern
	alpha90[which(alpha90<0)] <- alpha90[which(alpha90<0)] + 2*pi
	
	out <- cbind(alpha90,r)
	return(out)
	
}


####################################################################
# Die Funktion plot.vecs.func() soll die Vektoren des Fischs plotten.
#
# Übergeben werden muss:
#	- merkmale.out:	Ergebnis-Liste aus merkmale.func()
#	- track.ddf:		Ergebnis-Liste aus track.func() oder track.split()
#								oder track.split.over()
#	- track.eval.out:	Ergebnis-Liste aus track.eval()
#	- n.angle:			Anzahl der Vektoren pro Viertel
#	- objects:			Vektor der zu plottenden Objekte, alle falls NULL
#	- t.grid:			Vektor, der die Zeitpunkte der analysierten Frames 
#								angibt
#	- norm:				TRUE: Schnauze ist rechts und X0.Grad = 1
#							FALSE: So, wie es im Film ist
#	- do.save:			Speichern?
#	- pfad:				Dateipfad für den Plot
#	- filename:			Optionaler Dateiname innerhalb des Pfades
#
# Ausgegeben wird:
#	- nur plot
#
####################################################################

plot.vecs.func <- function(merkmale.out,
									track.ddf,
									track.eval.out,
									n.angle, 
									objects=NULL,
									t.grid, 
									norm=FALSE,
									do.save=FALSE,
									pfad=NULL,
									filename=NULL){
	
	
	# Anzahl Objekte
	anz.obj <- track.ddf$anz.obj
	
	# Hotspots und Zeitpunkte zu jedem Objekt
	object.ids <- track.eval.out$object.ids
	
	# Merkmale
	merkmale.list <- merkmale.out$merkmale.list
		
	x.s <- vector(length=4*n.angle)
	y.s <- vector(length=4*n.angle)
	
	# Wenn keine zu plottenden Objekte übergeben werden: Alle plotten, 
	if(is.null(objects)) {
		objects <- 1:anz.obj
	}
	
	# Anzahl Zeilen und Spalten für den Plot bestimmen, hängt von
	#	Anzahl Zeitpunkte und Anzahl Objekte ab
	nRow <- length(objects)
	nCol <- min(length(t.grid),100) # Maximal 100 Spalten
	
	
	
	# Dateinamen festlegen
	if(do.save==TRUE){
	
		if(is.null(filename)){
			if(norm)	filename <- "Vectors_Norm"
			if(!norm)	filename <- "Vectors_True"
		}
		
		file.pfad <- paste(pfad, filename,".jpg", sep="")

		# Grafikparameter
		height <- nRow * 1
		width <-  nCol * 1
		
		jpeg(file.pfad,
			  width=width,
			  height=height, 
			  units="in", 
			  res=200, 
			  quality=100)
				
	}
	
	par(mfrow=c(nRow, nCol), mar=c(1.7,1.8,1.7,0.2)+0.2)
	
	# Alle Objekte durchgehen
	for(i in objects){

		# Merkmalsmatrix rausziehen
		merkmale <- merkmale.list[[i]]
		
		# Nur weiter, wenn Merkmale bestimmt wurden, also wenn es nicht
		#	wegen in track.eval() gefundener Probleme gar nicht erst versucht
		#	wurde. Dann ist merkmale =NULL.
		#
		# Sonst: Leere Plots
		# Mögliche Veränderung: Sonst weniger Zeilen und diese Objekte gar 
		#	nicht plotten. Zum besseren Vergleich lasse ich es aber erstmal
		#	so.
		if(!is.null(merkmale)){
		
			# Alle Zeitpunkte durchgehen
			for(t in 1:nCol){ #length(t.grid)){
	
				# Falls zu dem Zeitpunkt das Objekt nicht vorkommt: leerer Plot
				if(!is.element(t, object.ids[[i]][,2])) {
					plot(0,0,col="white", main=paste("t = ",t.grid[t],sep=""),axes=FALSE, xlab="", ylab="")
					box()
				}else{
					
					# Richtige Spalte rausfinden
					j <- which(object.ids[[i]][,2]==t)
					
					# Sonst: Stern plotten
					
					# Nicht normieren, also alles "echt"
					if(norm==FALSE){
						
						# Nur, wenn X0.Grad nicht 0 ist
						if(merkmale[12,j]>0){
						
							# ...wieder auf echte Längen umrechnen
							merkmale[13:(11+4*n.angle),j] <- 
								merkmale[13:(11+4*n.angle),j] * merkmale[12,j] 
						}
						# ... sonst einfach so lassen: Wurden dann in der Merkmals
						#	funktion ja nicht verändert
						
						# Wahre Winkel
						alphas <- merkmale[(12+4*n.angle):(11+8*n.angle),j]
					}
					
					# Normieren, also "Schnauze rechts" und deren Länge = 1
					if(norm==TRUE){
						
						# Nur, wenn X0.Grad nicht 0 ist
						if(merkmale[12,j]>0){
						
							# X0.Grad auf 1 setzen
							merkmale[12,j] <- 1
						
							# Winkel (nicht die wahren, sondern die, für die sie sich 
							#	ausgeben) im Gradmaß
							alphas <- c(90/n.angle * 0:((4*n.angle)-1))
						}
					}
						
					# Falls norm==TRUE und X0.Grad=0: leerer Plot
					if(norm==TRUE & merkmale[12,j]==0){
						
						plot(0,0,col="white", main=paste("t = ",t.grid[t],sep=""),axes=FALSE, xlab="", ylab="")
						box()
					}else{
						
						# Sonst normal weiter
						
						# Kartesische Koordinate berechnen
						for(p in 0:((4*n.angle)-1)){
							x.s[p+1] <- merkmale[12+p,j] * cos(alphas[p+1]*pi/180)
							y.s[p+1] <- merkmale[12+p,j] * sin(alphas[p+1]*pi/180)
						}
						
						# Bestimmung von xlim und ylim für den Plot, damit
						#	die Einheiten auf der x- und y-Achse immer gleich sind.
						xlim <- c(min(x.s), max(x.s))
						x.l <- diff(xlim)
						ylim <- c(min(y.s), max(y.s))
						y.l <- diff(ylim)
						if(x.l>y.l)  ylim <- mean(ylim) + c(-1,1) * x.l/2
						if(x.l<=y.l) xlim <- mean(xlim) + c(-1,1) * y.l/2
						
						# Stern plotten					
						plot(0,0, xlim=xlim, ylim=ylim, xlab="xs", ylab="ys",main=paste("t = ",t.grid[t],sep=""))
						points(x.s, y.s, pch=16)
						for(p in 0:((4*n.angle)-1)){
							lines(c(x.s[p+1],0), c(y.s[p+1],0))
						}				
					}				
				}
			}	
		}else{
			
			# Alle Zeitpunkte durchgehen
			for(t in 1:length(t.grid)){
				
				# Leerer Plot
				plot(0,0,col="white", main=paste("t = ",t.grid[t],sep=""),axes=FALSE, xlab="", ylab="")
				box()
			}
		}
	}
	
	if(do.save==TRUE){
		dev.off()	
	}	
}

