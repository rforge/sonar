#====================================================================
# Alternative Funktionen zum Einlesen der ddf-Dateien.
#		Hier wird auch Versionsnummer, win.length & win.start
#		eingelesen
#
# Datum: ?
# Autoren: Carolin Maier, Michael Windmann, Ludwig Bothmann
#====================================================================

require(hexView)


# num_beams<-as.numeric(readRaw(file=my.file, offset=0,nbytes=master_header_length,human="int")$fileNum[5])
# vers<-as.numeric(readRaw(file=my.file, offset=0,nbytes=master_header_length,human="int")$fileNum[4])
# num_sample<-as.numeric(readRaw(file=my.file, offset=0,nbytes=master_header_length,human="int")$fileNum[7])
# a<-readRaw(file=my.file, offset=master_header_length,nbytes=frame_header_length,human="int")$fileNum
# a<-readRaw(file=my.file, offset=master_header_length,nbytes=frame_header_length,human="int")$fileNum
# if(num_beams == 96){ #higher frequency
#   win.start<-a[14]*0.42
#   win.length<-c(1.25,2.5,5,10)[a[15]+1]
# }else{ #lower frequency
#   win.start<-a[14]*0.84
#   win.length<-c(5,10,20,40)[a[15]+1]
# }


# setwd(dir="../../Documents/Fischprojekt/R/")
# source("Funktionen_read_plot_ddf.R")

# 
# read.ddf.frame <- function(file, i=0){
#   # Im Prinzip ist es egal ob hier 512 oder 1024 fuer master_header_length steht, nur die ersten Eintraege werden ausgelesen
#   # waehle min(512,1024)
#   b<-as.numeric(readRaw(file=my.file, offset=0,nbytes=512,human="int")$fileNum) 
#   vers<-as.numeric(readRaw(file=my.file, offset=0,nbytes=512,human="int")$fileRaw)[4]
#   num_beams<-b[5]
#   num_sample<-b[7]
#   
#   
#   if(vers==4){
#     master_header_length <- 1024 
#     frame_header_length <- 1024 
#   }
#   
#   if(vers==3){
#     master_header_length <- 512
#     frame_header_length <- 256
#   }
#   
#   
#   
#   frame_length <- num_sample*num_beams
#   
#   my.offset <- master_header_length + frame_header_length+(frame_header_length+frame_length)*i
#   bla <- readRaw(file=file,offset=my.offset,nbytes=frame_length,human="int" ,machine="binary",endian="littel"  )
#   this.frame <- as.numeric(bla$fileRaw)
#   frame1 <- matrix(this.frame,byrow=T,ncol=num_beams,nrow=num_sample)
#   
#   return(frame1)
# }



#ohne num_beams und version
#read.ddf.frame(file=my.file,i=0)

#' @import hexView
read.ddf.frame.new <- function(file, i=0, num_beams=96, vers, num_sample){

  
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
  bla <- readRaw(file=file,offset=my.offset,nbytes=frame_length,human="int" ,machine="binary",endian="littel"  )
  this.frame <- as.numeric(bla$fileRaw)
  frame1 <- matrix(this.frame,byrow=T,ncol=num_beams,nrow=num_sample)
  
  return(frame1)
}


#' @import hexView
read.ddf.allframes.new <- function(file,
                               winkel=seq(from=0.15,to=14.25,by=0.3), 
                               frames=NULL){
  
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
  
  ###########################################################################
  #neuer Teil
  #vers, beams, win.start,win.length automatisch auslesen aus .ddf Datei
  #globale Parameter aus dem master_header_File
  vers<-as.numeric(readRaw(file=file, offset=0,nbytes=512,human="int")$fileRaw)[4] #Versionsnummer
  b<-as.numeric(readRaw(file=file, offset=0,nbytes=512,human="int")$fileNum) 
  beams<-b[5] #Anzahl an Beams
  num_sample<-b[7] #Anzahl an Samples je Beam
  
  #Schauen ob es funktioniert:
  #print(paste("Version: ",vers))
  #print(paste("Beams: ",beams))
  #print(paste("Samples: ",num_sample))
  #############################################################################
  
  if(vers==4){
    master_header_length <- 1024 
    frame_header_length <- 1024 
  }
  
  if(vers==3){
    master_header_length <- 512
    frame_header_length <- 256
  }
  
  ##############################################################################################
  #Bestimmen von win.start und win.length abhaengig von der Frequenz (hoch oder niedrig)
  #Nicht aus master header sondern aus frame header
  a<-readRaw(file=file, offset=master_header_length,nbytes=frame_header_length,human="int")$fileNum
  if(beams == 96){ #higher frequency
    win.start<-a[14]*0.42
    win.length<-c(1.25,2.5,5,10)[a[15]+1]
  }else{ #lower frequency
    win.start<-a[14]*0.84
    win.length<-c(5,10,20,40)[a[15]+1]
  }
  
  #Schauen ob es funktioniert:
  #print(paste("win.start ",win.start))
  #print(paste("win.length ",win.length))
  # Bestimmung der Koordinaten der Sample Points
  ###############################################################################################
  sample.dist <- seq(from=win.start,to=win.start+win.length,length=512)
  
  hilfe <- sin( c(winkel)*pi/180 )
  xs <- -outer(sample.dist, c(-rev(hilfe),hilfe) )
  
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
      all.frames[[paste("F",i,sep="")]] <- read.ddf.frame.new(file=file,
                                                          i=i, 
                                                          vers=vers,
                                                          num_sample=num_sample)
    }
  }
  
  if(!is.null(frames)) {
    for(i in frames){
      all.frames[[paste("F",i,sep="")]] <- read.ddf.frame.new(file=file,
                                                          i=i, 
                                                          vers=vers,
                                                          num_sample=num_sample)
    }
  }	
  
  return(list(all.frames=all.frames, win.start=win.start,win.length=win.length,vers=vers))
}

#' Extract meta data from .ddf file
#' 
#' This function extracts meta information from a given sonar video .ddf file
#' 
#' @param file File name of sonar video (.ddf)
#' @return List with \code{vers}, the version of the .ddf file format,
#'  \code{win.start}, the start point of the sonar window in meters from the lense,
#'  \code{win.length}, the length of the sonar window in meters, \code{beams},
#'  the number of beams, \code{num_sample}, the number of pixels per beam and
#'  \code{n.frames}, the number of frames/images.
#' @export
#' @import hexView
get.version <- function(file){
	
	vers<-as.numeric(readRaw(file=file, offset=0,nbytes=512,human="int")$fileRaw)[4] #Versionsnummer
	b<-as.numeric(readRaw(file=file, offset=0,nbytes=512,human="int")$fileNum) 
	beams<-b[5] #Anzahl an Beams
	num_sample<-b[7] #Anzahl an Samples je Beam
	
	if(vers==4){
		master_header_length <- 1024 
		frame_header_length <- 1024 
	}
	
	if(vers==3){
		master_header_length <- 512
		frame_header_length <- 256
	}
	
	##############################################################################################
	#Bestimmen von win.start und win.length abhaengig von der Frequenz (hoch oder niedrig)
	#Nicht aus master header sondern aus frame header
	a<-readRaw(file=file, offset=master_header_length,nbytes=frame_header_length,human="int")$fileNum
	if(beams == 96){ #higher frequency
		win.start<-a[14]*0.42
		win.length<-c(1.25,2.5,5,10)[a[15]+1]
	}else{ #lower frequency
		win.start<-a[14]*0.84
		win.length<-c(5,10,20,40)[a[15]+1]
	}
	
	# Anzahl Frames (-1)
	n.frames <- get.max_num.frame(file=file)
	
	out <- data.frame(vers=vers,
										win.start=win.start,
										win.length=win.length,
										beams=beams,
										num_sample=num_sample,
										n.frames=n.frames)
	
	return(out)
}
