###########################################################################
# Fischprojekt
# Funktionen zur Vorbereitung der Daten auf die Klassifikation
# Datum der ersten Version: 10.09.12
# Letzte Aenderung: 5.11.12
# Autor: Ludwig Bothmann
###########################################################################

# require(ape)

####################################################################
# Die Funktion collect.hot() erstellt den Datensatz auf Hotspotebene
#
# Uebergeben werden muss:
#	- sum.prepro:						Ergebnis aus summary.func()
#	- merkmale.out.list:		Aus Preprocessing
#	- fisch:								ja / nein
# - art:									Aal, Forelle, Treibgut
#	- n.angle:							Anzahl Uhrzeiger pro Viertel
#	- fisch.short:						Kurzname der Fischart bzw Treibgut,
#												damit die IDs eindeutig sind
#
# Ausgegeben wird:
#	- data.hot:							Der Datensatz als data.frame
#
####################################################################


collect.hot <- function(sum.prepro,
								merkmale.out.list,
								fisch,
								art,
								n.angle,
								fisch.short){
	
	n.hots <- sum.prepro[,5]
	ID <- rep(sum.prepro[,1], times=n.hots)
	p.iter <- rep(sum.prepro[,2], times=n.hots)
	obj.iter <- rep(sum.prepro[,3], times=n.hots)
	
	
	# 1. Spalte: IDs
	data.hot <- data.frame(ID = ID)
	
	# 2.Spalte: Nummer des Hotspots innerhalb ID
	data.hot$t[1] <- 1
	
	# Falls es ueberhaupt mehr als eine ID gibt
	
	if(nrow(data.hot)>1){
		for(i in 2:nrow(data.hot)){			
			if(data.hot$ID[i]==data.hot$ID[(i-1)]){
				data.hot$t[i] <- data.hot$t[(i-1)]+1
			}else{
				data.hot$t[i] <- 1
			}
		}
	}
	
	# 3. + 4. Spalte 
	data.hot$fisch <- factor(fisch)
	data.hot$art <- factor(art)
	
	# 5. Spalte und folgende: Merkmale
	data.hot[,5:(4+5+4*n.angle)] <- 0
	
	do.loop <- TRUE
	
	for(i in 1:nrow(data.hot)){

# 		print(i)
		
		p <- p.iter[i]
		obj <- obj.iter[i]
		t <- data.hot[i,2]		
		data.hot[i,5:(4+5+4*n.angle)] <- merkmale.out.list[[p]][[1]][[obj]][c(1:5,12:(11+4*n.angle)),t]
		
		# ACHTUNG: SO IST 0 GRAD DIE WAHRE LAeNGE UND ALLES ANDERE IST
		#		IM VERHAeLTNIS DAZU ANGEGEBEN, ERSCHEINT MIR ABER SINNVOLL
		#		ZU SEIN
	}
	
	# Variablennamen uebergeben
	names(data.hot)[5:(4+5+4*n.angle)] <- make.names(as.vector(rownames(merkmale.out.list[[p]][[1]][[obj]][c(1:5,12:(11+4*n.angle)),])))
	names(data.hot)[5:6] <- c("Breite","Laenge")
	
	# ID aendern in "fisch_ID" damit es nachher eindeutig ist, da
	#	in verschiedenen Filmen die selben IDs vergeben werden koennen
	data.hot$ID <- paste(fisch.short,data.hot$ID, sep="_")
	
	return(data.hot)
	
}


####################################################################
# Die Funktion data.hotspot() erstellt den Datensatz auf Hotspotebene.
#		Hierfuer wird nicht der Input von summary.func() benoetigt. Es kann also
#		mit dieser Funktion erst der Datensatz auf Hotspotebene und danach der 
#		auf Objektebene berechnet werden. 
#
# Uebergeben werden muss:
# - track.ddf.list:
# - track.eval.out.list:
# - t.grid.list:
#	- merkmale.out.list:		Aus Preprocessing
#	- fisch:								ja / nein
# - art:									Aal, Forelle, Treibgut
#	- n.angle:							Anzahl Uhrzeiger pro Viertel
#	- fisch.short:						Kurzname der Fischart bzw Treibgut,
#												damit die IDs eindeutig sind
#
# Ausgegeben wird:
#	- data.hot:							Der Datensatz als data.frame.
#
# Anmerkung: Hier werden die wahren Laengen der Uhrzeiger ausgegeben,
#		wenn man etwas anderes moechte, muss man das nachher noch aendern.
#
####################################################################


data.hotspot <- function(track.ddf.list,
												 track.eval.out.list,
												 merkmale.out.list,
												 t.grid.list,
												 fisch,
												 art,
												 n.angle,
												 fisch.short){
	
	# N = # Uhrzeiger
	N <- n.angle*4
	
	# Anzahl Objekte insgesamt:
	anz.obj.total <- 0
	
	for(p in 1:length(t.grid.list)){
		
		if(!is.null(track.ddf.list[[p]]$anz.obj)){
			anz.obj.total <- anz.obj.total + track.ddf.list[[p]]$anz.obj
		}
	}
	
	########
	# Track.ev.long sind alle track.evs untereinander, das heisst:
	#  	- Objektnummer innerhalb Paket
	#		- unique = TRUE / FALSE
	#		- Anzahl Hotspots
	#		- Problem?
	
	track.ev.long <- NULL
	for(p in 1:length(t.grid.list)){
		
		if(!is.null(track.ddf.list[[p]]$anz.obj)){
			track.ev.long <- rbind(track.ev.long,
														 cbind( rep(p,track.ddf.list[[p]]$anz.obj), 
														 			 track.eval.out.list[[p]]$track.ev))
			
		}
	}
	
	# header enthaelt zusaetzlich die global durchnummerierten Objekte
	header <- cbind(1:anz.obj.total, track.ev.long)
	
	colnames(header)[1:6] <- c("ObjectID", "p", "objectInP", "noprob", 
																"n.hots", "problem")
	
	# Hilfsvariablen
	n.hots <- header[,5]
	ID <- rep(header[,1], times=n.hots)
	p.iter <- rep(header[,2], times=n.hots)
	obj.iter <- rep(header[,3], times=n.hots)
	
	
	#########
	# Ab hier wird der Datensatz nach und nach aufgebaut
	#########
	
	
	# 1. Spalte: IDs
	data.hot <- data.frame(ID = ID)
	
	# 2.Spalte: Nummer des Hotspots innerhalb ID
	data.hot$t[1] <- 1
	
	# Falls es ueberhaupt mehr als eine ID gibt
	
	if(nrow(data.hot)>1){
		for(i in 2:nrow(data.hot)){			
			if(data.hot$ID[i]==data.hot$ID[(i-1)]){
				data.hot$t[i] <- data.hot$t[(i-1)]+1
			}else{
				data.hot$t[i] <- 1
			}
		}
	}
	
	# 3. + 4. Spalte 
	data.hot$fisch <- factor(fisch)
	data.hot$art <- factor(art)
	
	# 5. Spalte und folgende: Merkmale
	data.hot[,5:(4+5+4*n.angle)] <- 0
	
	do.if <- TRUE
	
	# Alle Hotspots durchgehen
	for(i in 1:nrow(data.hot)){
		
# 		print(i)
		
		p <- p.iter[i]
		obj <- obj.iter[i]
		t <- data.hot[i,2]		
		merks <- merkmale.out.list[[p]][[1]][[obj]][c(1:5,12:(11+4*n.angle)),t]
		
		# Falls 0.Grad-Zeiger > 0: Uhrzeiger zurueckrechnen auf echte Werte
		if(!is.null(merks)){
			if(merks[6]>0){
				merks[-(1:6)] <- merks[-(1:6)]*merks[6]
			}
		
		data.hot[i,5:(4+5+4*n.angle)] <- merks
		}
		
		# Variablennamen einmal merken
		if(do.if==TRUE & length(make.names(as.vector(rownames(merkmale.out.list[[p]][[1]][[obj]][c(1:5,12:(11+4*n.angle)),]))))>0){
			var.names <- make.names(as.vector(rownames(merkmale.out.list[[p]][[1]][[obj]][c(1:5,12:(11+4*n.angle)),])))
			do.if <- FALSE
		}
		# ACHTUNG: So haben alle wieder ihre wahre Laenge
	}
	
	# Variablennamen uebergeben
	names(data.hot)[5:(4+5+4*n.angle)] <- var.names
	names(data.hot)[5:6] <- c("Breite","Laenge")
	
	# ID aendern in "fisch_ID" damit es nachher eindeutig ist, da
	#	in verschiedenen Filmen die selben IDs vergeben werden koennen
	data.hot$ID <- paste(fisch.short,data.hot$ID, sep="_")
	
	return(data.hot)
	
}

####################################################################
# Die Funktion data.hot.st.zentr() standardisiert die Laengen der Uhrzeiger
#		auf Summe = 1 und gibt die Summe zusaetzlich mit aus.
#		Ausserdem werden die Uhrzeiger danach zentriert.
#
# Uebergeben werden muss:
# - data.hot.roh:		Der unstandardisierte Datensatz als data.frame.
# - center:					TRUE => Zentrieren, sonst nur auf Summe 1 
#											standardisieren
#
# Ausgegeben wird:
#	- data.hot.zentr:	Der standardisierte Datensatz als data.frame.
#
####################################################################


data.hot.st.zentr <- function(data.hot.roh,
															center=TRUE){
	
	# Alle Spalten ausser den ersten 9 beinhalten die Uhrzeiger:
	watch.hands <- data.hot.roh[,-(1:9)]

	# Summe berechnen
	Summe <- apply(watch.hands, MARGIN=1, sum)
	
	# Durch Summe teilen
	watch.hands.st <- watch.hands / Summe
	
	# Summe als Variable hinzufuegen
	data.hot.st <- data.hot.roh
	data.hot.st[,-(1:9)] <- watch.hands.st
	data.hot.st$Summe <- Summe
	
	# Zentrieren (nicht die ersten 4 Spalten)
	data.hot.zentr <- data.hot.st
	
	if(center){
		data.hot.zentr[,-(1:4)] <- t(t(data.hot.zentr[,-(1:4)]) - 
			apply(data.hot.zentr[,-(1:4)], MARGIN=2, mean, na.rm=TRUE))
	}
	
	# Test, ob das funktioniert hat:
	#	round(apply(data.hot.zentr[,-(1:4)], MARGIN=2, mean),10)
	# => Ja
	
# 	# X.330Grad entfernen, da perfekte Abhaengigkeit zu den anderen Uhrzeigern
# 	ncols <- ncol(data.hot.zentr)
# 	data.hot.zentr <- data.hot.zentr[,-(ncols-1)]
	
	# Datensatz ausgeben
	return(data.hot.zentr)
}


####################################################################
# Die Funktion obj.base() soll die Designmatrix fuer das Baselinemodell
#		berechnen
#
# Uebergeben werden muss:
# - data.hot:		Datensatz auf Hotspotebene
#	- long:				optional: FALSE = nur Uhrzeiger, TRUE = Alle Kovs
# - only.watch:	TRUE: Es werden Laenge und Breite basierend auf
#									watch hands benutzt => Max statt Mean
# - FUN:				Funktion, mit der die Elemente zusammengefasst werden 
#									sollen, default = mean
# - FUN.ow:			Funktion, mit der die Elemente der ersten beiden
#									Spalten zusammengefasst werden sollen, wenn only.watch
#									=TRUE und long=TRUE, default=max
#
#
# Ausgegeben wird:
#	- baseline:		Designmatrix fuer Baselinemodell
#
####################################################################


obj.base <- function(data.hot, 
										 long=FALSE,
										 only.watch=FALSE,
										 FUN=mean,
										 FUN.ow=max){

	# IDs der Objekte
	ID <- unique(data.hot$ID)
	
	# Anzahl Objekte
	n.ID <- length(ID)
	
	# Nummern der Spalten bestimmen, deren Mittelwerte ich berechnen will
	if(long){
		cols <- c(1:ncol(data.hot))[-(1:4)]
	}else{
		cols <- c(1:ncol(data.hot))[-(1:9)]
	}

	# Designmatrix initialisieren
	baseline <- matrix(0, ncol=length(cols), nrow=n.ID)
	colnames(baseline) <- colnames(data.hot)[cols]
	rownames(baseline) <- ID
	
	# Alle Kovariablen durchgehen
	for(i in 1:length(cols)){
		
		# Spaltennummer
		j <- cols[i]
		
		if(only.watch & long){
			if(i %in% 1:2){
				
				# Max statt mean
				baseline.i <- tapply(data.hot[,j], 
														 INDEX = data.hot$ID,
														 FUN=FUN.ow, 
# 														 FUN=mean, # mit mean ist es dasselbe 
														 #					wie mit only.watch=FALSE
														 na.rm=TRUE)
			}else{
				
				# Mittelwert pro ID bestimmen
				baseline.i <- tapply(data.hot[,j], 
														 INDEX = data.hot$ID,
														 FUN=FUN, 
														 na.rm=TRUE)
			}
		}else{
			
			# Mittelwert pro ID bestimmen
			baseline.i <- tapply(data.hot[,j], 
													 INDEX = data.hot$ID,
													 FUN=FUN, 
													 na.rm=TRUE)
		}
		
		# Abspeichern
		baseline[,i] <- baseline.i
	}
	
	# X.330 bzw. letzten rausschmeissen. Das entspricht der vorletzten 
	#		Spalte, weil in der letzten Spalte die Summe steht.
	baseline <- baseline[,-(length(cols)-1)]
	
	# In data.frame umwandeln
	data.base <- data.frame(baseline)
	
	# Erst in Matrix umwandeln 
	#	(nur noetig, falls nur ein Objekt vorhanden ist
	if(is.vector(baseline)){
		baseline.m <- matrix(baseline, ncol=length(baseline))
		colnames(baseline.m) <- names(baseline)
		rownames(baseline.m) <- ID
		
		# In data.frame umwandeln
		data.base <- data.frame(baseline.m)
	}
		
	return(data.base)
}
	
####################################################################
# Die Funktion obj.sdcorr() soll die Designmatrix fuer das Modell
#		mit Standardabweichungen und Korrelationen berechnen
#
# Uebergeben werden muss:
# - data.hot:		Datensatz auf Hotspotebene
#	- long:				optional: FALSE = nur Uhrzeiger, TRUE = Alle Kovs
#								ACHTUNG: TRUE bisher nicht implementiert
#
# Ausgegeben wird:
#	- sdcorr:		Designmatrix fuer SdCorr-Modell
#
####################################################################


obj.sdcorr <- function(data.hot, 
											 long=FALSE){

	# IDs der Objekte
	ID <- unique(data.hot$ID)
	
	# Anzahl Objekte
	n.ID <- length(ID)
	
	if(long){
		cat("Achtung: Option long=TRUE noch nicht implementiert \n => long=FALSE durchgefuehrt \n")
	}
	
	# Nummern der Spalten bestimmen, deren Sd und Korr ich berechnen will
	cols <- c(1:ncol(data.hot))[-c(1:9,ncol(data.hot))]
	
	# Anzahl Kovariablen
	N <- ncol(data.hot) - 10
	
	# Anzahl Kovariablen pro Viertel
	n.angle <- N/4
	
	# Anzahl Spalten in der Designmatrix
	anz.col <- N 	+						# fuer die SDs
						N*(N-1)/2 	# fuer die Korrelationen
	
	# Designmatrix initialisieren
	sdcorr <- matrix(0, ncol=anz.col, nrow=n.ID)
	
	# Zeilen- und Spaltennamen
	rownames(sdcorr) <- ID
	colnames(sdcorr) <- 1:anz.col
	colnames(sdcorr)[1:N] <- paste("Stdev_",90/n.angle * 0:(N-1), "_Grad", sep="")
	namen <- paste(rep(1:N,N),"_",rep(1:N,each=N),sep="")
	S <- matrix(namen, ncol=N)
	w <- which(upper.tri(S)==TRUE, arr.ind=TRUE)
	colnames(sdcorr)[(N+1):anz.col] <- S[w]
	
	# Alle IDs durchgehen
	for(i in 1:length(ID)){
	
		# ID des aktuellen Objekts
		ID.i <- ID[i]
	
		# Welche Hotspots gehoeren zu diesem Objekt?
		hots <- which(data.hot$ID==ID.i)
		
		# Standardabweichungen berechnen
		stdevs <- apply(data.hot[hots, cols], 2, sd, na.rm=TRUE)
		
		# Korrelationen berechnen
		corrs <- cor(data.hot[hots, cols])
		
		# Abspeichern
		sdcorr[i,1:N] <- stdevs
		sdcorr[i,-(1:N)] <- corrs[w]
			
	}
	
	# In data.frame umwandeln
	data.sdcorr <- data.frame(sdcorr)

	return(data.sdcorr)	
}
	

####################################################################
# Die Funktion obj.svd() soll die Designmatrix fuer das SVD-Modell
# 	berechnen
#
# Uebergeben werden muss:
# - data.hot:		Datensatz auf Hotspotebene
#	- long:				optional: FALSE = nur Uhrzeiger, TRUE = Alle Kovs
#								ACHTUNG: TRUE bisher nicht implementiert
# - nu, nv:			Parameter fuer SVD, default = 3
#
# Ausgegeben wird:
#	- svd.design:		Designmatrix fuer SVD-Modell
#
####################################################################


obj.svd <- function(data.hot, 
										long=FALSE,
										nu=3,
										nv=3){
	
	# IDs der Objekte
	ID <- unique(data.hot$ID)
	
	# IDs fuer alle Hotspots
	ID.hot <- data.hot$ID
	
	# Anzahl Objekte
	n.ID <- length(ID)
	
	if(long){
		cat("Achtung: Option long=TRUE noch nicht implementiert \n => long=FALSE durchgefuehrt \n")
	}
	
	# Nummern der Spalten bestimmen, in denen die Uhrzeiger stehen
# 	cols <- c(1:ncol(data.hot))[-c(1:9,ncol(data.hot))]
	
	# UPDATE 6.5.13: X.330 raus und dafuer Summe rein
	cols <- c(1:ncol(data.hot))[-c(1:9,(ncol(data.hot)-1))]
	
	# Anzahl Uhrzeiger
	N <- length(cols)
	
	# Designmatrix initialisieren
	svd.design <- matrix(0, ncol=N*3, nrow=n.ID)
	
	# Zeilen- und Spaltennamen
	rownames(svd.design) <- ID
	colnames(svd.design) <- c(paste("DA",(0:(N-2))*360/N,sep="."),"DA.Sum",
														paste("DF",(0:(N-2))*360/N,sep="."),"DF.Sum",
														paste("DT",(0:(N-2))*360/N,sep="."),"DT.Sum")
	
	# Alt:
# 	colnames(svd.design) <- c(paste("DA",(0:(N-1))*360/N,sep="."),
# 														paste("DF",(0:(N-1))*360/N,sep="."),
# 														paste("DT",(0:(N-1))*360/N,sep="."),)
	
	# Matrix mit Originaluhrzeigern
	X <- data.hot[,cols]
	
	########################
	# 1. SVD durchfuehren
	########################
	
	# Nur Zeilen mit Aalen, Forellen, Treibgut
	X.aal 		<- data.hot[which(data.hot[,4]=="Aal"),cols]
	X.forelle <- data.hot[which(data.hot[,4]=="Forelle"),cols]
	X.treib 	<- data.hot[which(data.hot[,4]=="Treibgut"),cols]
	
	# Singulaerwertzerlegung
	svd.aal 		<-  svd(X.aal,     nu=nu, nv=nv)
	svd.forelle <-  svd(X.forelle, nu=nu, nv=nv)
	svd.treib 	<-  svd(X.treib,   nu=nu, nv=nv)
	
	# Spalteneffekte
	V.aal 		<- svd.aal$v
	V.forelle <- svd.forelle$v
	V.treib 	<- svd.treib$v
	
	################################
	# 2. Lineare Modelle schaetzen
	################################
	
	# Brauche ich nicht, kann ich gleich in Prognose machen
	
	########################
	# 3. X approximieren
	########################
	
	X.hat.aal 		<- as.matrix(X) %*% V.aal 		%*% t(V.aal)
	X.hat.forelle <- as.matrix(X) %*% V.forelle %*% t(V.forelle)
	X.hat.treib 	<- as.matrix(X) %*% V.treib 	%*% t(V.treib)
	
	########################
	# 4. Scores berechnen
	########################
	
	# Differenzenmatrizen
	Diff.aal 			<- X - X.hat.aal
	Diff.forelle 	<- X - X.hat.forelle
	Diff.treib 		<- X - X.hat.treib
	
	# Pro Objekt Differenzen quadrieren und aufsummieren (und durch Anzahl
	#		Hotspots teilen). 
	# Fuer jeden Uhrzeiger separat.
	
	# Uhrzeiger durchgehen
	for(k in 1:length(cols)){
	
		# Objekte durchgehen
		for(i in 1:length(ID)){
		
			# Differenzen quadrieren und aufsummieren fuer Aal, Forelle, Treibgut			
			d.aal 		<- mean(Diff.aal[which(ID.hot==ID[i]),k]^2)
			d.forelle <- mean(Diff.forelle[which(ID.hot==ID[i]),k]^2)
			d.treib 	<- mean(Diff.treib[which(ID.hot==ID[i]),k]^2)
			
			# Abspeichern
			svd.design[i,k] 					<- d.aal
			svd.design[i,(k + N)] 		<- d.forelle
			svd.design[i,(k + N*2)] 	<- d.treib
			
		}
	
	}
	
	#################################
	# 5. Resultierende Designmatrix
	#################################
	
	# Schon in 4. passiert
	
	# In data.frame umwandeln
	data.svd <- data.frame(svd.design)
	
	return(data.svd)
}


####################################################################
# Die Funktion collect.hot.zentr() erstellt den Datensatz auf Hotspotebene.
#	Dies ist eine Version von collect.hot(), wo zusaetzlich die Zentroide
#	abgespeichert werden.
#
# Uebergeben werden muss:
#	- sum.prepro:						Ergebnis aus summary.func()
#	- merkmale.out.list:		Aus Preprocessing
#	- fisch:								ja / nein
# - art:									Aal, Forelle, Treibgut
#	- n.angle:							Anzahl Uhrzeiger pro Viertel
#	- fisch.short:						Kurzname der Fischart bzw Treibgut,
#												damit die IDs eindeutig sind
#
# Ausgegeben wird:
#	- data.hot:							Der Datensatz als data.frame
#
####################################################################


collect.hot.zentr <- function(sum.prepro,
								merkmale.out.list,
								fisch,
								art,
								n.angle,
								fisch.short){
	
	n.hots <- sum.prepro[,5]
	ID <- rep(sum.prepro[,1], times=n.hots)
	p.iter <- rep(sum.prepro[,2], times=n.hots)
	obj.iter <- rep(sum.prepro[,3], times=n.hots)
	
	
	# 1. Spalte: IDs
	data.hot <- data.frame(ID = ID)
	
	# 2.Spalte: Nummer des Hotspots innerhalb ID
	data.hot$t[1] <- 1
	
	# Falls es ueberhaupt mehr als eine ID gibt
	
	if(nrow(data.hot)>1){
		for(i in 2:nrow(data.hot)){			
			if(data.hot$ID[i]==data.hot$ID[(i-1)]){
				data.hot$t[i] <- data.hot$t[(i-1)]+1
			}else{
				data.hot$t[i] <- 1
			}
		}
	}
	
	# 3. + 4. Spalte 
	data.hot$fisch <- factor(fisch)
	data.hot$art <- factor(art)
	
	# 5. Spalte und folgende: Merkmale & Zentroide
	data.hot[,5:(4+5+4*n.angle + 2)] <- 0
	
	do.loop <- TRUE
	
	for(i in 1:nrow(data.hot)){
		
		# 		print(i)
		
		p <- p.iter[i]
		obj <- obj.iter[i]
		t <- data.hot[i,2]		
		data.hot[i,5:(4+5+4*n.angle)] <- merkmale.out.list[[p]][[1]][[obj]][c(1:5,12:(11+4*n.angle)),t]
		
		# ACHTUNG: SO IST 0 GRAD DIE WAHRE LAeNGE UND ALLES ANDERE IST
		#		IM VERHAeLTNIS DAZU ANGEGEBEN, ERSCHEINT MIR ABER SINNVOLL
		#		ZU SEIN
		
		data.hot[i,(4+5+4*n.angle+1):(4+5+4*n.angle+2)] <-  merkmale.out.list[[p]][[1]][[obj]][(10:11),t]
	}
	
	# Variablennamen uebergeben
	names(data.hot)[5:(4+5+4*n.angle)] <- make.names(as.vector(rownames(merkmale.out.list[[p]][[1]][[obj]][c(1:5,12:(11+4*n.angle)),])))
	names(data.hot)[(4+5+4*n.angle+1):(4+5+4*n.angle+2)] <- c("zentr.x","zentr.y")
	names(data.hot)[5:6] <- c("Breite","Laenge")
	
	# ID aendern in "fisch_ID" damit es nachher eindeutig ist, da
	#	in verschiedenen Filmen die selben IDs vergeben werden koennen
	data.hot$ID <- paste(fisch.short,data.hot$ID, sep="_")
	
	return(data.hot)
	
}

####################################################################
# Die Funktion collect.hot.diff() erstellt einen Datensatz auf 
#	Hotspotebene. Hier werden die Differenzen zu den gleichen Vektoren der
#	jeweils naechsten Hotspots berechnet. Pro Objekt ist also eine
#	Zeile weniger im Enddatensatz als bei collect.hot
#
# Uebergeben werden muss:
#	- sum.prepro:						Ergebnis aus summary.func()
#	- merkmale.out.list:		Aus Preprocessing
#	- fisch:								ja / nein
# - art:									Aal, Forelle, Treibgut
#	- n.angle:							Anzahl Uhrzeiger pro Viertel
#	- fisch.short:						Kurzname der Fischart bzw Treibgut,
#												damit die IDs eindeutig sind
#
# Ausgegeben wird:
#	- data.hot:							Der Datensatz als data.frame
#
####################################################################


collect.hot.diff <- function(track.ddf.list,
									  track.eval.out.list,
									  merkmale.out.list,
									  t.grid.list,
									  n.angle,
									  sum.prepro,
									  art,
									  fisch,
									  fisch.short,
									  quad=FALSE){
	
	n.hots <- sum.prepro[,5]
	ID <- rep(sum.prepro[,1], times=n.hots-1)
	p.iter <- rep(sum.prepro[,2], times=n.hots-1)
	obj.iter <- rep(sum.prepro[,3], times=n.hots-1)
	ID.uni <- sum.prepro[,1]
	
	# 1. Spalte: IDs
	data.hot <- data.frame(ID = ID)
	
	# 2.Spalte: Nummer des Hotspots innerhalb ID
	data.hot$t[1] <- 1
	
	# Falls es ueberhaupt mehr als eine ID gibt
	
	if(nrow(data.hot)>1){
		for(i in 2:nrow(data.hot)){			
			if(data.hot$ID[i]==data.hot$ID[(i-1)]){
				data.hot$t[i] <- data.hot$t[(i-1)]+1
			}else{
				data.hot$t[i] <- 1
			}
		}
	}
	
	# 3. + 4. Spalte 
	data.hot$fisch <- factor(fisch)
	data.hot$art <- factor(art)
	
	# 5. Spalte und folgende: Merkmale
	data.hot[,5:(4+4*n.angle)] <- 0
	
	
	# N = # Uhrzeiger
	N <- n.angle*4
	
	# Anzahl Objekte insgesamt:
	anz.obj.total <- length(ID.uni)
	
	# Track.ev.long sind alle track.evs untereinander
	track.ev.long <- NULL
	for(p in 1:length(t.grid.list)){
		
		if(!is.null(track.ddf.list[[p]]$anz.obj)){
			track.ev.long <- rbind(track.ev.long,
										  cbind( rep(p,track.ddf.list[[p]]$anz.obj), 
										  		 track.eval.out.list[[p]]$track.ev))
			
		}
	}
	
	track.ev.long <- track.ev.long[ID.uni,]
	
	# klass.out matcht alles moegliche zueinander
	klass.out <- cbind(ID.uni, track.ev.long)
	
	colnames(klass.out)[1:6] <- c("ObjectID", "p", "objectInP", "noprob", 
											"n.hots", "problem")
	
# 	# Merkmale dazu
# 	# 	klass.out[7:(N*2+6] <- 0
# 	anz.col <- (N*2+6+(N*(N-1)/2) + 4)
# 	klass.out[7:anz.col] <- 0
# 	colnames(klass.out)[7:(N+6)] <- paste("Mean_",90/n.angle * 0:((4*n.angle)-1), "_Grad", sep="")
# 	colnames(klass.out)[(N+7):(N*2+6)] <- paste("Stdev_",90/n.angle * 0:((4*n.angle)-1), "_Grad", sep="")
# 	namen <- paste(rep(1:N,N),"_",rep(1:N,each=N),sep="")
# 	S <- matrix(namen, ncol=N)
# 	w <- which(upper.tri(S)==TRUE, arr.ind=TRUE)
# 	colnames(klass.out)[(N*2+7):(anz.col-4)] <- S[w]
# 	
# 	colnames(klass.out)[(anz.col-3):anz.col] <- c("alpha.mean", "alpha.sd", 
# 																 "rad.mean", "rad.sd")
# 	
# 	vel <- vector(length=anz.obj.total)
	
	hot.diff<- NULL
	
	for(i in 1:anz.obj.total){
		
		# 		print(i)
		
		# Nur wenn kein Problem und mehr als ein Hotspots
		if(klass.out[i,4]==TRUE & klass.out[i,5]>1){
			
			# Zeitraum
			p <- klass.out[i,2]
			
			# Objekt innerhalb Zeitraum
			j <- klass.out[i,3]
			
			merks <- matrix(merkmale.out.list[[p]][[1]][[j]][12:(11+4*n.angle),], 
								 nrow=N)
			
# 			# Fuer alle Hotspots dieses Objekts...
# 			for(k in 1:ncol(merks)){
# 				# ... Uhrzeiger zurueckrechnen auf echte Werte
# 				if(merks[1,k]>0){
# 					merks[-1,k] <- merks[-1,k]*merks[1,k]
# 				}
# 			}
			
			# Differenz berechnen
			if(quad==FALSE){
				
				merks.diff <- apply(merks,1,diff)
				
			}else{
				
				# Quadrierte Differenzen
				merks.diff <- apply(merks,1,diff)^2
				
			}
			
		}
		
		hot.diff <- rbind(hot.diff, merks.diff)
	}
	
	# ID aendern in "fisch_ID" damit es nachher eindeutig ist, da
	#	in verschiedenen Filmen die selben IDs vergeben werden koennen
# 	ID <- paste(fisch.short,ID, sep="_")
	
	hot.diff.out <- data.frame(ID=ID, hot.diff)
	
	hot.diff.out$ID <- paste(fisch.short,hot.diff.out$ID, sep="_")
		
	colnames(hot.diff.out)[2:(N+1)] <- paste("Diff_",90/n.angle * 0:((4*n.angle)-1), "_Grad", sep="")
	
	return(hot.diff.out)
	
}


####################################################################
# Die Funktion prepare.and.count() soll die preprocesseten Daten
#		soweit vorbereiten, dass danach direkt gezaehlt werden kann
#		und zaehlen
#
# Uebergeben werden muss:
#	- 
#
# Ausgegeben wird:
#	- out.list:			Liste mit Objekten, die ich zum erstellen des
#										Outputs benoetige
#
####################################################################

#' Prepare data for counting and count
#' 
#' This function prepares the output of the previous steps (features...) for
#' the classification, classifies the objects and counts number of objects
#' per class
#' 
#' @param track.ddf.list Tracking, output of \code{\link{preprocess.data}}
#' @param track.eval.out.list Tracking evaluation, output of \code{\link{preprocess.data}}
#' @param merkmale.out.list Features, output of \code{\link{preprocess.data}}
#' @param t.grid.list Grids of time points, output of \code{\link{preprocess.data}}
#' @param n.angle Number of watch hands per quarter, i.e., total number of watch hands is \code{4 x n.angle}
#' @param min.length Minimal length of tracks. Tracks with less hotspots are deleted.
#' @param center \code{TRUE}: Center variables on hotspot level, default is \code{FALSE}
#' @param long \code{TRUE}: Average of all hotspot variables are used, default is \code{TRUE}
#' @param only.watch \code{TRUE}: Use watch hands instead of trapezoid variables - 
#'  strongly recommended, default is \code{TRUE}
#' @param file.regel File name of classification rule for fish species
#' @param klass.regel Type of classification rule for fish species, i.e., lda, qda...
#' @param sdcorr \code{TRUE}: Standard deviations and correlations of variables 
#'  on hotspot level are used, default is \code{FALSE}
#' @param FUN Function for gathering hotspot information on object level, default is \code{mean}
#' @param FUN.ow Function for gathering hotspot information on object level 
#'  when not using the trapezoid variables, default is \code{mean}
#' @param lda.obj Object of learned classification rule if lda is used, is 
#'  loaded from \code{file.regel}
#' @param qda.obj Object of learned classification rule if qda is used, is 
#'  loaded from \code{file.regel}
#' @param svm.obj Object of learned classification rule if svm is used, is 
#'  loaded from \code{file.regel}
#' @return Result of classification
#' @export
prepare.and.count <- function(track.ddf.list,
															track.eval.out.list,
															merkmale.out.list,
															t.grid.list,
															n.angle,
															min.length,
															center,
															long,
															only.watch,
															file.regel="Regel_2.RData",
															klass.regel,
															sdcorr,
															FUN,
															FUN.ow,
															lda.obj=NULL,
															qda.obj=NULL,
															svm.obj=NULL){

	
	
	##################################################
	# a) Datensaetze auf Hotspotebene
	#
	# => Uhrzeiger haben hier noch ihre wahren Laengen
	# => Es sind alle Objekte enthalten
	#
	# Funktion, die das macht: data.hotspot()
	##################################################
		
	# Datensatz auf Hotspotebene erstellen
	# MW: Welche der "fixen" Variabeln sind notwendig?
	# LB: Die Funktion braucht irgendwas, es ist aber egal, was 
	data.hot.roh1 <- data.hotspot(track.ddf.list=track.ddf.list,
																track.eval.out.list=track.eval.out.list,
																merkmale.out.list = merkmale.out.list,
																t.grid.list=t.grid.list,
																fisch="ja",
																art="Aal",
																n.angle=n.angle,
																fisch.short="a1")
	
	# Summary.func rechnen um die Motionvariablen zu bekommen
	sum.prepro.roh1 <- summary.func(track.ddf.list=track.ddf.list,
																	track.eval.out.list=track.eval.out.list,
																	merkmale.out.list=merkmale.out.list,
																	t.grid.list=t.grid.list,
																	n.angle=n.angle)
	
	# MW: Warnungen sollten abgefangen werden.
	
	# Nur die Motionvariablen abspeichern
	motion.roh1 <- sum.prepro.roh1[,-(1:(ncol(sum.prepro.roh1)-5))]
	rownames(motion.roh1) <- paste("a1",sum.prepro.roh1$ObjectID, sep="_")
	
	
	#########################################
	# b) Objekte selektieren
	#########################################
	
	# IDs auswaehlen: Mindestens min.length Hotspots und unproblematisch
	sum.prepro.roh <- sum.prepro.roh1[which(sum.prepro.roh1$n.hots>=min.length 
																					& sum.prepro.roh1$noprob),]
	
	# NAs loeschen
	nas <- unique(which(is.na(sum.prepro.roh),arr.ind=TRUE)[,1])
	
	if(length(nas)>0){
		sum.prepro.roh <- sum.prepro.roh[-nas,]
	}else{
		sum.prepro.roh <- sum.prepro.roh
	}
	
	# Welche Objekte sollen ausgewaehlt werden?
	which.select <- sum.prepro.roh$ObjectID
	
	# Nur weiter, wenn mindestens eins ausgewaehlt wird
	if(length(which.select)>0){
		
		# Wie heissen die IDs von den Objekten?
		ID.select <- paste("a1_",which.select,sep="")
		
		# Nur die Hotspots dieser Objekte auswaehlen
		data.hot.roh.select <- data.hot.roh1[which(is.element(data.hot.roh1$ID,
																													ID.select)),]
		
		# Nur die Motionvariablen dieser Objekte
		motion.roh.select <- motion.roh1[which.select,]
		
		# 	# Die wahre Klasse dieser Objekte
		# 	klasse.roh1 <-  data.frame(ID=ID.select, art="Aal")
		# 	# MW: Gibt es im Echtzeitsystem nicht
		
		
		#########################################
		# c) Datensaetze auf Hotspotebene standardisieren 
		#				und zentrieren wie in Text 6, Kap.2 
		#				beschrieben
		#########################################
		
		data.hot.1 <- data.hot.st.zentr(data.hot.roh = data.hot.roh.select,
																		center=center)	
		
		
		
		########################################################
		# d) Designmatrizen fuer die Modelle auf Objektebene
		########################################################
		
		# 	# MW: Wird das zum Problem?
		#	# LB: Auskommentiert, wahre Klasse kennen wir ja nicht
		# 	klasse.1 <- data.frame(klasse=klasse.roh1$art)
		# 	rownames(klasse.1) <- klasse.roh1$ID
		
		
		##################################################################
		# Baseline, SdCorr, Motion & SVD separat berechnen
		##################################################################
		
		################
		# i) Baseline
		################
		
		baseline.1 <- obj.base(data.hot = data.hot.1, 
													 long=long,
													 only.watch=only.watch,
													 FUN=FUN,
													 FUN.ow=FUN.ow)
		
		# Option long=TRUE: Auch die Mittelwerte von Laenge, Breite, Flaeche,
		#		Seitenverhaeltnis und Volumen werden berechnet
		# Option only.watch: Wenn Variablen basierend auf watch hands benutzt
		#		werden fuer Laenge, Breite, Flaeche, Seitenverhaeltnis,
		#		dann wird von Laenge und Breite die Funktion FUN.ow statt FUN 
		#		genommen. (Z.B. Max statt mean)
		#		=> Greift nur, wenn long=TRUE
		
		################
		# ii) SdCorr
		################
		
		if(sdcorr){
			sdcorr.1 <- obj.sdcorr(data.hot = data.hot.1, 
														 long=FALSE)
		}
		
		
		################
		# iii) Motion
		################
		
		motion.1 <- motion.roh.select
		
		
		# 	################
		# 	# iv) SVD
		# 	################
		# 	
		# 	svd.1 <- obj.svd(data.hot = data.hot.1, long=FALSE)
		
		
		########################################################
		# e) "Test"-datensatz fuer die Klassifikation
		########################################################
		
		if(!sdcorr){
			# 		data.gesamt <- data.frame(klasse.1, baseline.1, motion.1)
			data.gesamt <- data.frame(baseline.1, motion.1)
		}else{
			# 		data.gesamt <- data.frame(klasse.1, baseline.1, sdcorr.1, motion.1)
			data.gesamt <- data.frame(baseline.1, sdcorr.1, motion.1)
		}
		
# 		data.gesamt <- data.roh
		
		###########################################################################
		# 3.) Zaehlung
		###########################################################################
		
		# Klassifikationsregel laden
		
		load(file=file.regel)
		
		if(klass.regel=="lda"){
			
			#################
			# LDA
			#################
			
			pred.obj <- predict(lda.obj, data.gesamt)	
			pred.class <- pred.obj$class
			post.obj <- pred.obj$post
			
		}else if(klass.regel=="qda"){
			
			#################
			# QDA
			#################
			
			pred.obj <- predict(qda.obj, data.gesamt)	
			pred.class <- pred.obj$class
			post.obj <- pred.obj$post
			
		}else if(klass.regel=="svm"){	
			
			#################
			# SVM
			#################
			
			pred.obj <- predict(svm.obj, data.gesamt, probability=TRUE)	
			pred.class <- predict(svm.obj, data.gesamt)
			post.obj <- attr(pred.obj,"probabilities")
			
		}
		
		# 			result.obj <- cbind(as.numeric(pred.class), data.gesamt[,1])
		# 			error.obj <- c(sum(result.obj[,1]!=result.obj[,2]), 
		# 								nrow(data.gesamt))
		# 			
		# 			colnames(result.obj) <- c("Prognose","Wahr")
		# 			error.obj
		
		#########################
		# Kreuztabellen
		#########################
		
		# MW: Kreuztabellen koennen im Echtzeitsystem auch nicht erstellt werden?
		# LB: Stimmt, habe es geloescht
		
		# Output:
		# MW: Welche Klassifizierungsmethode?
		# LB: Kann jetzt variabel in Parameter_Global uebergeben werden
		
# 		zaehlung <- factor(as.numeric(pred.class), levels=c(1,2,3),
# 											 labels=c("Eel","Trout","Debris"))		
		# 		table(zaehlung)
		
		# MW: Variablen Pfad angeben, erledigt.
		
# 		save(zaehlung, file=paste(pfad.mult ,"zaehlung.Rdata",sep="") )
		
		# Ende Klassifikation
		do.classification <- TRUE
		out.list <- list(do.classification=do.classification,
										 which.select=which.select,
										 pred.class=pred.class,
										 post.obj=post.obj,
										 sum.prepro.roh=sum.prepro.roh,
										 data.gesamt=data.gesamt)
		
	}else{
		
		# Sonst wurde nichts ausgewaehlt, dann soll auch ein leerer Output
		#		erstellt werden
		do.classification <- FALSE
		out.list <- list(do.classification=do.classification,
										 which.select=NULL,
										 pred.class=NULL,
										 post.obj=NULL,
										 sum.prepro.roh=NULL)
		
	}

	return(out.list)
}

# ###########################################################################
# # Aus dem Buch "Multivariate Analysemethoden" von Handl, bzw. von der 
# #		Homepage dazu, leicht abgeaendert
# ###########################################################################
# 
# plotmst<-function(X,e,data)
# {# erstellt minimal spannenden Baum fuer die Datenmatrix X  und zeichnet ihn in 
# 	# das Streudiagramm der Scores bezueglich der ersten beiden Hauptkomponenten
# 	# in e steht das Ergebnis der Hauptkomponentenanalyse
# 	# Dies kann man zum Beispiel durch e<-princomp(X) erhalten
# 	n<-dim(X)[1]
# 	m<-mst(dist(X))
# 	v<-(1:n)[m[1,]==1]
# 	welche<-cbind(rep(1,length(v)),v)
# 	for(i in 2:(n-1)) {v<-(1:n)[m[i,]==1];
# 										 welche<-rbind(welche,cbind(rep(i,length(v)),v))}
# 	x<-e$scores[,1]
# 	y<-e$scores[,2]
# 	plot(x,y,xlab="1. HK",ylab="2. HK",main="", pch=19, col=as.numeric(data$art)+1)
# # 	text(x,y,dimnames(X)[[1]],pos=1)
# 	segments(x[welche[,1]],y[welche[,1]],x[welche[,2]],y[welche[,2]])
# }