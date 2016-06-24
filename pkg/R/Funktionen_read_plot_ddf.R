################################################################################
# Fischprojekt
# Funktionen zum Einlesen und plotten von ddf-Dateien
# Datum der ersten Version: 15.05.12
# Letzte Änderung: 12.03.13
# Autor: Ludwig Bothmann, basierend auf R-Code von Michael Windmann
################################################################################



require(hexView)

##############################################################
# Funktionen um .ddf-Files einzulesen:
# Gilt für Versionen des Formats DDF_04!
# Update: Jetzt auch Format DDF_03 unterstützt (vers=3)
##############################################################



##############################################################
# Die Funktion get.max_num.frame() ermittelt die Nummer des 
#	letzten Frames einer ddf-Datei. 
#
# Warnung: absolutes Max ist noch  zu bestimmen. Bisher auf 
#	256*256 Frames pro File beschränkt
#
# Übergeben werden muss:
#	- file:				Dateiname der ddf-Datei
#
# Ausgegeben wird:
#	- max.frame.num:	Die Nummer des letzten Frames
##############################################################

#' @import hexView
get.max_num.frame <- function(file){

	# Header einlesen
	file.header <- readRaw(file=file, offset=4  ,nbytes=2 ,human="int"    )
		
	# Zahl auslesen
	anzahl <- as.numeric(file.header$fileRaw)
	max.frame.num <- anzahl[1] + anzahl[2]*256 - 1
	
	# Zahl zurückgeben
	return(max.frame.num)
}

##############################################################
# Die Funktion read.ddf.frame() liest einen bestimmten Frame 
#	einer ddf-Datei ein.
#
# Übergeben werden muss:
#	- file: 		Der Dateiname der ddf-Datei
#	- i:			Der Index des einzulesenden Frames (beginnend bei 0!)
#	- num_beams:	Die Anzahl der Beams: Nur 96 oder 48 möglich
#	- vers:			Version des ddf-Files (3 oder 4)
#
# Ausgegeben wird:
#	- frame1:		Eine Matrix mit den Werten des Signals
#
##############################################################

#' @import hexView
read.ddf.frame <- function(file, i=0, num_beams=96, vers){
	num_sample <- 512
	
	if(vers==4){
		master_header_length <- 1024 
		frame_header_length <- 1024 
	}
	
	if(vers==3){
		master_header_length <- 512
		frame_header_length <- 256
	}

	frame_length <- num_sample*num_beams
	
	my.offset <- master_header_length + frame_header_length+(frame_header_length+frame_length)*i
	bla <- readRaw(file=file,
						offset=my.offset,
						nbytes=frame_length,
						human="int" ,
						machine="binary",
						endian="littel"  )
	this.frame <- as.numeric(bla$fileRaw)
	frame1 <- matrix(this.frame,
						  byrow=T,
						  ncol=num_beams,
						  nrow=num_sample)
	
	return(frame1)
}

##############################################################
# Die Funktion read.ddf.allframes() liest eine komplette 
#	ddf-Datei ein.
#
# Übergeben werden muss:
#	- file: 		Der Dateiname der ddf-Datei
#	- beams: 		Anzahl der Beams (96 oder 48)
#	- win.start:	Start-Entfernung des Ultraschallfensters von der Linse
#	- win.length:	Länge des Ultraschallfensters
#	- winkel:		Gitter mit Winkeln der Beams, eine Hälfte
#	- frames:		Indizes der zu lesenden Frames 
#						(optional, wenn nichts übergeben wird, 
#						wird die ganze Datei eingelesen)
#	- vers:			Version des ddf-Files (3 oder 4)
#
# Ausgegeben wird eine Liste mit den Einträgen:
#	- Num_Frames:	Anzahl der eingelesenen Frames
#	- xs:			Kartesische x-Koordinaten der Signalpunkte
#	- ys:			Kartesische y-Koordinaten der Signalpunkte
#	- F0 ... Fmax:	Matrizen mit Intensitäten der Echos an den 
#						Signalpunkten (Signalstärken)
#
##############################################################



read.ddf.allframes <- function(file,beams=96,
										 win.start=1.250,
										 win.length=5.00,
										 winkel=seq(from=0.15,to=14.25,by=0.3), 
										 frames=NULL, 
										 vers=4){
	
	# Nummer des letzten Frames bestimmen
	if(is.null(frames)) {
		#frame.names <- paste("f",0:max.frame,sep="")
		max.frame <- get.max_num.frame(file=file)
	}
	
	# ... bzw. Anzahl der einzulesenden Frames
	if(!is.null(frames)) {
		#frame.names <- paste("f",frames,sep="")
		max.frame <- length(frames)
	}
	
	
	# Bestimmung der Koordinaten der Sample Points
	
	sample.dist <- seq(from=win.start,to=win.start+win.length,length=512)
	
	hilfe <- sin( c(winkel)*pi/180 )
# 	xs <- -outer(sample.dist, c(-rev(hilfe),hilfe) )
	
	# 17.07.13 Änderung nach Kommentar von Carolin Maier: Vorzeichen umdrehen
	#	=> Dann müssten die Bilder so aussehen wir im DIDSON-Viewer
	# 22.11.13: Wenn doch das Minus davor steht, passt der Output zu dem 
	#		im DIDSON-Viewer. Wenn ich das mache, muss ich aber auch die Regel
	#		noch mal neu lernen, sonst passt die Klassifikation nicht mehr!
	xs <- outer(sample.dist, c(-rev(hilfe),hilfe) )
		
	hilfe<- cos(winkel*pi/180)
	ys <- outer(sample.dist,c(rev(hilfe),hilfe) )
	
	
	# Anzahl Frames & kartesische Koordinaten abspeichern
	if(is.null(frames)) {
		all.frames <- list(Num_Frames=max.frame+1,xs=xs,ys=ys)
	}
	
	if(!is.null(frames)) {
		all.frames <- list(Num_Frames=max.frame,xs=xs,ys=ys)
	}
	
	#######################################
	# Einlesen der einzelnen Frames
	#######################################
	
	if(is.null(frames)) {
		for(i in 0:max.frame){
			all.frames[[paste("F",i,sep="")]] <- read.ddf.frame(file=file,
																				 i=i, 
																				 vers=vers)
		}
	}
	
	if(!is.null(frames)) {
		for(i in frames){
			all.frames[[paste("F",i,sep="")]] <- read.ddf.frame(file=file,
																				 i=i, 
																				 vers=vers)
		}
	}	
	
	return(all.frames)
}

################################################################
# Die Funktion plot.ddf_frame() plottet einen Frame 
#	eines ddf-Files.
#
# Übergeben werden muss:
#	- data:			Ergebnisliste aus der Funktion read.ddf.allframes()
#	- frame.num:	Nummer des zu plottenden Frames
#	- bw:			Schwarz/weiss vertauschen: 1=tauschen, also Hintergrund weiß
#
# Ausgegeben wird: 
#	- nur der Plot
#
################################################################

plot.ddf_frame <- function(data,frame.num=0,bw=1,...){

	# Koordinaten der Punkte
	x <- data[["xs"]]
	y <- data[["ys"]]
	
	plot(x,y,type="n")
	frame <- data[[paste("F",frame.num,sep="")]]
	graustufen <- abs(bw - (frame/255))
	points(x,y,col=grey(graustufen),pch=16,cex=0.5)
}

################################################################
# Die Funktion plot.ddf() plottet ebenfalls eine ddf-Datei, 
#	allerdings wir hier nicht das ddf-"Objekt" übergeben, sondern
# 	direkt die Koordinaten der Punkte und die zu plottenden Werte.
#
# Übergeben werden muss:
# 	- data:	 	Matrix mit den zu plottenden Werten
# 	- xs, ys: 	Koordinaten der Punkte
# 	- bw: 		1: Vertauschen von Schwarz/Weiss  
# 	- bin:		TRUE => übergeben wurde 0/1 - Bild
#	- cex:		Größe der Punkte
# 	- t:			Zeitpunkt, wird nur für den Titel benutzt
#	- scale:		Sollen die Graustufen auf (0,255) skaliert werden? Nur nötig,
#						wenn bin == FALSE
#	- col:		Optional Farbe der Punkte, wenn bin == TRUE
#	- main:		Optional main für den Plot
#
# Ausgegeben wird:
#	- nur der Plot
#
################################################################


plot.ddf <- function(data, 
							bw=1,
							xs,
							ys,
							bin=FALSE,
							cex=0.5,
							t=1,
							scale=TRUE, 
							col=NULL,
							main=NULL,
							cex.main=1,
							cex.axis=1,
							...){
	

	if(!is.null(dim(data))){
		
		if(bin==FALSE){
			if(scale==TRUE){
				# Korrektur, falls die Werte nicht in (0,255) liegen
				
				# if(min(data)<0) data[which(data<0)] <- 0
				# if(max(data)>255) data[which(data>255)] <- 255
				
				# Falls das Bild nur 0er enthält ist das Maximum 0 und nach 
				#	dem nächsten Schritt sind nur noch NaN im Bild, deshalb 
				#	nur, wenn nicht nur 0er drin sind
				if(sum(data)!=0){
			
					# Alternativ nicht abschneiden, sondern skalieren:
					data <- data - min(data) 		# Jetzt ist das Minimum bei 0
					data <- data/max(data)*255		# Jetzt ist das Maximum bei 255
				}
			}
		}
		
		# Bestimmung von xlim und ylim für den Plot, damit
		#	die Einheiten auf der x- und y-Achse immer gleich sind.
		xlim <- c(min(xs), max(xs))
		x.l <- diff(xlim)
		ylim <- c(min(ys), max(ys))
		y.l <- diff(ylim)
		
		if(x.l>y.l)  ylim <- mean(ylim) + c(-1,1) * x.l/2
		if(x.l<=y.l) xlim <- mean(xlim) + c(-1,1) * y.l/2
		
		# Main festlegen
		if(is.null(main)){
			main <- paste("t = ",t,sep="")
		}
				
		# Plot
		plot(xs,ys,type="n", xlim=xlim, ylim=ylim, 
			  main=main, cex.main=cex.main, cex.axis=cex.axis,...)
		
		# Falls 0/1-Bild übergeben wird: Nur schwarz und weiß zeichnen
		# Falls nicht: Normal zeichnen, also Graustufen zwischen 0 und 255
		
		# FRAGE: Max auf 255 und Min auf 0 setzen? 
		#	Also die Werte strecken, so dass alle Graustufen benutzt werden?
		#	=> Wird jetzt oben gemacht, wenn scale=TRUE
		
		if(bin==FALSE){
			graustufen <- abs(bw - (data/255))
			points(xs,ys,col=grey(graustufen),pch=16,cex=cex) # pch=45 ?
		}
				
		
		if(bin==TRUE){
			graustufen <- abs(bw - data)
				if(is.null(col)){
				points(xs,ys,col=grey(graustufen),pch=16,cex=cex)
			}else{
				points(xs,ys,col=graustufen*col,pch=16,cex=cex)
				}
			
		}
	}
	
	# Falls nur ein Skalar übergeben wird: Leerer Plot
	if(is.null(dim(data))){
		
		if(is.null(main)){
			plot(0,0,col="white", main=paste("t = ",t,": No Fish",sep=""),
					 axes=FALSE, xlab="", ylab="",cex.main=cex.main,...)
			box()
		}else{
			plot(0,0,col="white", main=main,
					 axes=FALSE, xlab="xs", ylab="ys",cex.main=cex.main,...)
			box()
		}
			 
	}

}

################################################################
################################################################
# Beispielaufrufe:
# pfad.basis <- "/home/bothmannlu/Documents/"
# my.file <- paste(pfad.basis,"Fischprojekt/Videos_Daten/DDFs/2012-04-04_190001_HF_Frame_5881-5995.ddf",sep="")
# ddf.data <- read.ddf.allframes(file=my.file)
# plot.ddf_frame(data=ddf.data,frame.num=54)
# 
# y <- ddf.data$F54
# 
# image(t(y))
# image(y)
# 
# par(mfrow=c(1,2))
# 	plot.ddf_frame(data=ddf.data,frame.num=54)
# 	image(t(y))


##############################################################