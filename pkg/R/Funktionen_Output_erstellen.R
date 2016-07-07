###########################################################################
# Fischprojekt
# Funktionen zur Erstellung des Outputs am Ende der Fischzaehlung
# Datum der ersten Version: 18.07.13
# Letzte aenderung: 22.11.13
# Autoren: Ludwig Bothmann & Michael Windmann
###########################################################################


####################################################################
# Die Funktion generate.output() erstellt den Output am Ende der
#		Fischzaehlung
#
# uebergeben werden muss:
#	- do.classification:	Wurde etwas gezaehlt? Sonst Zeile mit Nullen
# which.select,
# pred.class,
# post.obj,
# sum.prepro.roh,
# frames.pack,
# track.eval.out.list,
# merkmale.out.list,
# max.frame,
# pfad.mult,
# zeitstempel
# 
# Ausgegeben wird:
#	- nichts, es wird direkt eine .csv-Datei abgespeichert
#
####################################################################

#' Generate output of fish classification
#' 
#' This function generates the output of the entire fish classification
#' 
#' @param do.classification \code{TRUE}: Results of classification are available
#' @param which.select Vector of indices of objects
#' @param pred.class Predicted classes
#' @param post.obj Posterior probabilities
#' @param sum.prepro.roh Summary of the preprocessing, output of 
#'  \code{\link{summary_func}}
#' @param frames.pack Number of frames to be analyzed in one package
#' @param track.eval.out.list Tracking result, output of \code{\link{preprocess.data}}
#' @param merkmale.out.list Features, output of \code{\link{preprocess.data}}
#' @param max.frame Total number of frames \code{max.frame} of the given video, 
#'  possibly extracted before via \code{\link{get.version}} 
#' @param pfad.mult Folder for resulting plots and output
#' @param zeitstempel time stamp of .ddf file
#' @param vars.hot.list Variables on hotspot level, output of \code{\link{preprocess.data}}
#' @param is.schwarm.vec Vector regarding shoal of fish
#' @param only.schwarm \code{TRUE}: Only shoal analysis was carried out, 
#'  default is \code{FALSE}
#' @return Output table with class of fish and all features is saved at 
#' \code{pfad.mult}
#' @export
generate.output <- function(do.classification,
														which.select,
														pred.class,
														post.obj,
														sum.prepro.roh,
														frames.pack,
														track.eval.out.list,
														merkmale.out.list,
														max.frame,
														pfad.mult,
														zeitstempel,
# 														schwarm.test=FALSE,
# 														max.pix.vec,
# 														max.hots.vec,
														vars.hot.list,
														is.schwarm.vec,
														only.schwarm=FALSE
														){
	
	# Nur Schwaerme gefunden
	if(only.schwarm){
		
		#Filename des Outputs
		out.filename <- paste(pfad.mult,"output_",zeitstempel,".csv",sep="")
		
		output <- matrix(0,nrow=length(vars.hot.list),ncol=21 )
		
		colnames(output) <- c("Art","Laenge","Breite","Flaeche",
													"Geschwindigkeit","Bewegungsrichtung",
													"Ausrichtung","ersterFrame", 
													"letzterFrame","ersterSPx","ersterSPy",
													"letzterSPx","letzterSPy",
													"PWS.Aal","PWS.Forelle","PWSTreibgut",
													"max.hots","max.lebend","max.Laenge",
													"max.Breite","max.Flaeche")
	
		for(i in 1:length(vars.hot.list)){
		
			if(is.schwarm.vec[i]){
				output[i,"Art"] <- 4
				output[i,"Flaeche"] <- vars.hot.list[[i]]$max.flaeche
				output[i,"max.hots"] <- vars.hot.list[[i]]$max.hots
			}else{
				output[i,"Art"] <- 99
				output[i,-1] <- 0
			}
		}
				
# 		colnames(output) <- c("Max.pix","Max.hots")		
		
		# Output ausgeben:
		write.table(output, file=out.filename, sep=";",dec=",", row.names=F,col.names=T)
		cat(paste("Output-Datei",out.filename,  "geschrieben.", sep=" "),"\n")
		
		}else{
			
		if(do.classification){
		
			output <- matrix(NA, 
											 nrow=length(which.select), 
											 ncol=21)
			
			colnames(output) <- c("Art","Laenge","Breite","Flaeche",
														"Geschwindigkeit","Bewegungsrichtung",
														"Ausrichtung","ersterFrame", 
														"letzterFrame","ersterSPx","ersterSPy",
														"letzterSPx","letzterSPy",
														"PWS.Aal","PWS.Forelle","PWSTreibgut",
														"max.hots","max.lebend","max.Laenge",
														"max.Breite","max.Flaeche")
			
			rownames(output) <- which.select
			
			# Art
			# 1 = Aal
			# 2 = Fisch_unbestimmt
			# 3 = Treibgut
			# 4 = Schwarm (Kann hier aber nicht vorkommen)
			# 99 = Nichts gefunden (Kann hier aber nicht vorkommen)
			output[,"Art"] <- as.numeric(pred.class)
			
			# PosterioriWS
			output[,c("PWS.Aal","PWS.Forelle","PWSTreibgut")] <- post.obj
			
			
			for(MY.ID in which.select){
				
		# 		cat(MY.ID)	
				
				# Objektnummer in String umwandeln, dann kann ueber den Zeilennamen
				#		direkt auf die richtige Zeile in der Outputmatrix 
				#		zugegriffen werden
				char.ID <- as.character(MY.ID)
				
				# Erster und letzter Frame, in dem das Objekt auftaucht
				my.p <- sum.prepro.roh[char.ID,"p"]
				my.objectInP <- sum.prepro.roh[char.ID,"objectInP"]
				
				output[char.ID, c("ersterFrame" )] <- frames.pack*(my.p-1)  + min((track.eval.out.list[[my.p]]$object.ids[[my.objectInP]])[,"t"]	)
				
				output[char.ID, c("letzterFrame" )] <- frames.pack*(my.p-1) +  max((track.eval.out.list[[my.p]]$object.ids[[my.objectInP]])[,"t"]	)
				
				# Merkmale
				merkmale.matrix <- merkmale.out.list[[my.p]]$merkmale.list[my.objectInP]
				
				# Mean von Laenge und Breite
				output[char.ID, c("Laenge","Breite") ] <- apply(merkmale.matrix[[1]][c("Laenge (Dim.y)","Breite (Dim.x)") ,],
																																		MARGIN=1,
																																		FUN=mean,
																																		na.rm=TRUE)
				
				# LB: Mean oder Median oder Max?
				# => 19.11.13 Muenster: Mean und Max
				
				# Max von Laenge und Breite
				output[char.ID, c("max.Laenge","max.Breite") ] <- apply(merkmale.matrix[[1]][c("Laenge (Dim.y)","Breite (Dim.x)") ,],
																																		MARGIN=1,
																																		FUN=max,
																																		na.rm=TRUE)
				
				
				# Flaeche
				output[char.ID, c("Flaeche")] <- mean(merkmale.matrix[[1]][c("Volumen"),],
																													 na.rm=TRUE)
				
				# LB: Lieber das Volumen nehmen? 
				#		22.11.13: Ab jetzt statt Flaeche das Volumen
				# LB: Mean oder Median oder Max?
				#		22.11.13: Auch hier beides
				output[char.ID, c("max.Flaeche")] <- max(merkmale.matrix[[1]][c("Volumen"),],
																							na.rm=TRUE)
				
				# Ausrichtung: Hauptachse in Grad zur Horizontalen
				output[char.ID, c("Ausrichtung")] <- -mean(merkmale.matrix[[1]][c( "Hauptachse in Grad zur Horiz."),])
				# LB: 22.11.13:
				#		Komischerweise brauche ich hier ein Minus, damit es zu dem
				#		passt, was man auf dem Video sieht. Ich erinnere mich noch
				#		das mit den kartesischen Koordinaten etwas nicht passte, weiss
				#		aber nicht mehr genau, ob ich das jetzt richtig gemacht habe.
				#		Fuer die Klassifikation ist es aber egal und jetzt stimmt
				#		dieser Winkel hier auch.
				
				
				# Geschwindigkeit	
# 				x.d <- merkmale.matrix[[1]][c( "Marginaler Schwerpunkt Dim.x kart") , c(2:ncol(merkmale.matrix[[1]]) )  ] - merkmale.matrix[[1]][c( "Marginaler Schwerpunkt Dim.x kart") , c(1:(ncol(merkmale.matrix[[1]])-1)  )  ]
# 				y.d <- merkmale.matrix[[1]][c( "Marginaler Schwerpunkt Dim.y kart") , c(2:ncol(merkmale.matrix[[1]]) )  ] - merkmale.matrix[[1]][c( "Marginaler Schwerpunkt Dim.y kart") , c(1:(ncol(merkmale.matrix[[1]])-1)  )  ]
# 				
# 				weg <- sum(sqrt(x.d^2 + y.d^2) )
# 				zeit <- (output[char.ID,"letzterFrame"] - output[as.character(MY.ID),"ersterFrame"]+1)/10
# 				
# 				output[char.ID,"Geschwindigkeit"] <- weg/zeit
				
				# LB: Lieber mean.dist nehmen?
				#		22.11.13: rad.mean, also durchschnittliche zurueckgelegte 
				#		Strecke zwischen zwei Zeitpunkten MAL 10
				#		=> Einheit: Meter/Sekunde
				output[char.ID,"Geschwindigkeit"] <- merkmale.out.list[[my.p]]$dir.vel.list[my.objectInP][[1]]["rad.mean"]*10
				
				# Koordinaten-Schwerpunkte
				output[char.ID, c( "ersterSPx" , "letzterSPx") ] <- -merkmale.matrix[[1]][c( "Marginaler Schwerpunkt Dim.x kart") , c(1,ncol(merkmale.matrix[[1]]) )  ]
				output[char.ID, c( "ersterSPy" , "letzterSPy") ] <- merkmale.matrix[[1]][c( "Marginaler Schwerpunkt Dim.y kart") , c(1,ncol(merkmale.matrix[[1]]) )  ]
				
				# LB, 22.11.13: Auch hier brauche ich ein Minus, damit es zum 
				#		Video passt
				#	=> So ist das also alles ok. Etwas sauberer waere es vielleicht,
				#		die kartesischen x-Koordinaten in Funktionen_read_plot_ddf.R
				#		richtig (also mit Minus) zu berechnen, dann muss ich aber auch
				#		die Klassifikationsregel noch mal neu lernen.
				
# 				# Bewegungsrichtung in Grad zur Horizontalen		
				# 29.11.13: Wieder atan2, siehe weiter unten
# 				alpha <- merkmale.out.list[[my.p]]$dir.vel.list[my.objectInP][[1]]["alpha.mean"]
# 				
# 				# Vom Bogenmass in Gradmass umrechnen und so drehen, dass es der
# 				#		Winkel zur y-Achse ist und in (-180,180] ist
# 				alpha.grad <- alpha*180/pi + 90
# 				if(alpha.grad>180){
# 					alpha.grad <- alpha.grad-360
# 				}else if(alpha.grad<=-180){
# 					alpha.grad <- alpha.grad+360
# 				}
# 				
# 				# LB: 22.11.13: Im Moment muss ich den Winkel noch an der y-Achse
# 				#		spiegeln damit er zu dem auf dem Video passt
# 				alpha.grad <- sign(alpha.grad) * 180 - alpha.grad
# 				
# 				output[char.ID,"Bewegungsrichtung"] <- alpha.grad
			}
			

			# Bewegungsrichtung 
			output[,"Bewegungsrichtung"] <-
				atan2((-output[,"ersterSPy"]+output[,"letzterSPy"]) , (-output[,"ersterSPx"]+output[,"letzterSPx"])) * (180/pi)
			
			# LB: ueberpruefen, ob das stimmt
			#		=> Vorzeichen vertauscht, jetzt ist es glaube ich richtig
			#		=> Oder besser alpha.mean? 
			#				(aus merkmale.out.list[[my.p]]$dir.vel.list[my.objectInP])
			#		22.11.13: alpha.mean oben reinprogrammiert
			#		29.11.13: Wieder atan2() benutzen, bei dieser Mitteleung koennen
			#			komische Sachen passieren, es ist schoener, wenn die Bewegungs-
			#			richtung zu den angegebenen Punkten passt

				
			# Maximale Anzahl an Hotspots pro Frame
			n.hots <- rep(0,max.frame)
			
			for(n in 1:nrow(output)){
				
				von <- output[n,8]
				bis <- output[n,9]
				n.hots[c(von:bis)] <- n.hots[c(von:bis)] + 1
			}
			
			max.hots <- max(n.hots) 
			output[,"max.hots"] <- max.hots
			
			# => Hier ist jetzt in jeder Zeile das gleiche, wenn mehrere Objekte
			#		erkannt wurden
			
			# Maximale Anzahl an lebenden Objekten pro Frame 
			#	(Das gleiche wie max.hots, aber nur Aal und Fisch)
			n.lebend <- rep(0,max.frame)
			
			# Welche Zeilen gehoeren zu lebenden Objekten?
			n.grid <- which(is.element(output[,"Art"],c(1,2)))
			
			for(n in n.grid){
				
				von <- output[n,8]
				bis <- output[n,9]
				n.lebend[c(von:bis)] <- n.lebend[c(von:bis)] + 1
			}
			
			max.lebend <- max(n.lebend) 
			output[,"max.lebend"] <- max.lebend
			
			# => Hier ist jetzt in jeder Zeile das gleiche, wenn mehrere Objekte
			#		erkannt wurden

			
			###########################
			output <- round(output,4)
			
			
			## output muss noch weitere Leerspalten bekommen fuer weitere auszugebenden Variablen
			# Filename des Outputs
			out.filename <- paste(pfad.mult,"output_",zeitstempel,".csv",sep="")
			
			# Output ausgeben:
			write.table(output, file=out.filename, sep=";",dec=",", row.names=F,col.names=T,na=" ")
			cat(paste("Output-Datei",out.filename,  "geschrieben.", sep=" "),"\n")
	
		}else{
			
			#Filename des Outputs
			out.filename <- paste(pfad.mult,"output_",zeitstempel,".csv",sep="")
			
			output <- matrix(0,nrow=1,ncol=21 )
			
			# Wenn nichts gefunden wurde: Art = -99, da 0 schon Treibgut ist
			output[1] <- 99
			
			colnames(output) <- c("Art","Laenge","Breite","Flaeche",
														"Geschwindigkeit","Bewegungsrichtung",
														"Ausrichtung","ersterFrame", 
														"letzterFrame","ersterSPx","ersterSPy",
														"letzterSPx","letzterSPy",
														"PWS.Aal","PWS.Forelle","PWSTreibgut",
														"max.hots","max.lebend","max.Laenge",
														"max.Breite","max.Flaeche")
			
			
			# Output ausgeben:
			write.table(output, file=out.filename, sep=";",dec=",", row.names=F,col.names=T)
			cat(paste("Output-Datei",out.filename,  "geschrieben.", sep=" "),"\n")
			
		}
	}

}