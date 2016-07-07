###########################################################################
# Fischprojekt
# Funktion für das Preprocessing der Filme
# Datum der ersten Version: 18.10.12
# Letzte Änderung: 12.03.13
# Autor: Ludwig Bothmann
###########################################################################

require(parallel)

####################################################################
# Die Funktion preprocess.data() soll die Filme für die Klassifikation
#	vorbereiten.
#
# Übergeben werden muss:
#	- my.file:			Dateiname des .ddf-Films
#	- win.start:		Beginn des Sonarfensters
#	- win.length:		Länge des Sonarfensters
#	- vers:				Version des .ddf-Films	
#	- y.lim:				Soll etwas abgeschnitten werden? Dan < 512 setzen
#	- a.1:				Grenzwert a.1
#	- a.2:				Grenzwert a.2
#	- cut:				Minimale Größe einer Punktewolke
#	- c:					Anzahl (weiße) Pixel, die am Rand eines jeden Hotspots
#									mitabgespeichert werden (default = 0)
#	- m.d.cs:			Maximaler Abstand für das Tracking
#	- pfad.mult:		Dateipfad für die Plots
#	- n.angle:			Anzahl Vektoren pro Viertel
#	- signal.neg:		Ist das Signal der Objekte negativ?
#	- do.plot:			Hotspots plotten?
#	- do.plot.vec:		Vektoren plotten?
#	- n.cores:			Anzahl Rechenkerne
#	- which.packs:		NULL => Ganzer Film, sonst die Nummer der Pakete angeben,
#									welche analysiert werden sollen, ein Paket sind per
#									default 100 Frames, also 10 Sekunden
#	- df.t:				Anzahl Freiheitsgrade der B-Spline-Basis in Zeitrichtung,
#									NULL => Gleich Anzahl Zeipunkte
#	- frames.pack:		Wieviele Frames sollen zusammen analysiert werden?
#									def=100
#	- maxdist:			Maximaler Abstand zweier Pixel um zum selben Cluster
#									zu gehören
#	- fastclust:		Hotspots finden mit fast.clust()?
#	- speeditup:		Parameter für fast.clust()
#
# Ausgegeben wird:
#	Nichts, alles wird direkt hart abgespeichert
#
####################################################################

 
#' Preprocess sonar video
#' 
#' This function carries out the entire preprocessing of a sonar video, i.e.,
#' localizes hotspots and tracks objects over time
#' 
#' @param my.file Filename of .ddf file to be analyzed
#' @param win.start Start of the sonar window, possibly extracted before via 
#'  \code{\link{get.version}}
#' @param win.length Length of the sonar window, possibly extracted before via 
#'  \code{\link{get.version}}
#' @param vers Version of the .ddf file, possibly extracted before via 
#'  \code{\link{get.version}}
#' @param y.lim Number of pixels to analyze in each beam, counted from the 
#'  camera
#' @param a.1 Threshold for centered data
#' @param a.2 Threshold for uncentered data, default \code{0} means that no thresholding is done
#' @param cut Minimal size of cluster in number of pixels
#' @param c Number of white pixels as frames for the hotspot. Default \code{0} is a 
#'  good choice
#' @param m.d.cs Maximal distance for which two hotspots can be assigned the 
#'  same tracking number
#' @param pfad.mult Folder for resulting plots and output
#' @param n.angle Number of watch hands per quarter, i.e., total number of watch hands is \code{4 x n.angle}
#' @param signal.neg Is the signal negativ? Default is \code{FALSE}
#' @param do.plot \code{TRUE}: Plots are saved, default is \code{FALSE}
#' @param do.plot.vec \code{TRUE}: Plots are saved, default is \code{FALSE}
#' @param n.cores Number of cores if parallelization is needed, default is \code{1}
#' @param which.packs	Default \code{NULL} results in analysis of the entire 
#'  video. Alternatively, the indices of packages (per default 100 frames) to
#'  be analyzed
#' @param frames.pack Number of frames to be analyzed in one package
#' @param df.t Number of degrees of freedom for the splines. Default \code{NULL} 
#'  results in \code{c(100, 25, round(n.frames/2))} where \code{n.frames} is the 
#'  total number of frames of the analyzed video
#' @param maxdist Maximal distance of two pixels to be assigned to the same cluster, 
#'  not relevant for cluster methods floodclust and floodCpp
#' @param fastclust \code{TRUE}: fastclust method is used for clustering of 
#'  pixels and definition of hotspots, default is \code{FALSE}
#' @param speeditup \code{TRUE}: fast version of fastclust method is used for clustering of 
#'  pixels and definition of hotspots, default is \code{FALSE}
#' @param floodclust \code{TRUE}: floodclust method is used for clustering of 
#'  pixels and definition of hotspots, default is \code{FALSE}
#' @param floodCpp \code{TRUE}: floodCpp method is used for clustering of 
#'  pixels and definition of hotspots, default is \code{TRUE} - recommended method
#' @param plot.hots.mult \code{TRUE}: Plots are saved, default is \code{FALSE}
#' @param only.watch \code{TRUE}: Use watch hands instead of trapezoid variables - 
#'  strongly recommended, default is \code{TRUE}
#' @param schwarm.find \code{TRUE}: A detection of shoal of fish is carried out, 
#'  default is \code{FALSE}
#' @param regel.schwarm File name of classification rule for shoal of fish
#' @param which.regel.schwarm Type of classification rule for shoal of fish, i.e., lda, qda...
#' @param save.preprocess \code{TRUE}: Results of preprocessing are saved, default is \code{FALSE}
#' @param save.jpg \code{TRUE}: jpgs for the creating of the video are saved, 
#'  default is \code{TRUE}
#' @param zeitstempel time stamp of .ddf file
#' @param save.movie \code{TRUE}: Video is created, 
#'  default is \code{TRUE}
#' @param delete.jpg \code{TRUE}: jpgs for the creating of the video are deleted 
#'  after creating the video, default is \code{TRUE}
#' @param each Interval of images for the video, for example default \code{5} 
#'  means, that each 5th image is saved in the video
#' @param wait.jpg \code{TRUE}: Code pauses until creation of video is finished, 
#'  default is \code{FALSE}
#' @param win Specify \code{TRUE} if you are running under windows, then, 
#'  parallelization is not possible
#' @return List with various elements necessary for further analysis steps
#' @export
#' @import parallel
#' 
preprocess.data <- function(my.file,
                            win.start,
                            win.length,
                            vers,
                            y.lim,
                            a.1, 
                            a.2, 
                            cut,
                            c=0,
                            m.d.cs,
                            pfad.mult,
                            n.angle,
                            signal.neg,
                            do.plot,
                            do.plot.vec,
                            n.cores=1,
                            which.packs=NULL,
                            frames.pack=100,
                            df.t=NULL,
                            maxdist=1,
                            fastclust=FALSE,
                            speeditup=FALSE,
                            floodclust=FALSE,
                            floodCpp=TRUE,
                            plot.hots.mult=FALSE,
                            only.watch=FALSE,
                            schwarm.find=TRUE,
                            regel.schwarm,
                            which.regel.schwarm="qda",
                            save.preprocess=TRUE,
                            save.jpg=FALSE,
                            zeitstempel="zeitstempel",
                            save.movie=FALSE, 
                            delete.jpg=FALSE,
                            each=1,
                            wait.jpg=FALSE,
                            win
){
  
  if(floodclust & (n.cores>1)){
    
    # 		stop("Achtung: Verwendung von floodclust bei mehreren Kernen \n aktuell noch nicht möglich \n")
    cat("Achtung: Verwendung von floodclust bei mehreren Kernen \n aktuell mit Vorsicht zu benutzen \n")
  }
  
  # Anzahl Frames bestimmen um Anzahl der Teile zu bestimmen, in der 
  #	Film geteilt wird
  
  # Anzahl Frames - 1
  max.frame <- get.max_num.frame(file=my.file)
  
  # Wieviele Frames sind zu löschen, damit ich es in Pakete zu frames.pack 
  #	Frames aufteilen kann? Achtung, eigentlich beginnt die Zählung bei
  #	0, den Frame lasse ich auf jeden Fall weg
  n.delete <- max.frame%%frames.pack
  
  # Wie viele Pakete bleiben dann übrig?
  n.packs <- (max.frame-n.delete) / frames.pack
  
  if(is.null(which.packs)){
    which.packs <- 1:n.packs
  }
  
  # Anzahl der zu analysierenden Pakete
  n.p.total <- length(which.packs)
  
  # Liste mit allen track.evals
  track.evals <- vector("list",n.p.total)
  
  # Liste mit allen track.eval.out
  track.eval.out.list <- vector("list", n.p.total)
  
  # Liste mit allen track.ddfs
  track.ddf.list <- vector("list", n.p.total)
  
  # Liste mit allen merkmalen
  merkmale.out.list <- vector("list", n.p.total)
  
  # Liste mit allen t.grid
  t.grid.list <- vector("list", n.p.total)
  
  # Liste mit allen Pfaden (wobei die immer gleich sind, aber egal)
  pfad.mult.list <- vector("list", n.p.total)
  
  # Liste mit allen 0/1-Hotspots für alle Pakete und dort für jeden
  #	Zeitpunkt und in jedem Zeitpunkt für jeden Hotspot
  hots.y01.list <- vector("list", n.p.total)
  
  # Liste mit allen Dimensionen aller 0/1-Hotspots für alle Pakete 
  #	und dort für jeden
  #	Zeitpunkt und in jedem Zeitpunkt für jeden Hotspot
  dims.list <- vector("list", n.p.total)
  
  # 	# Vector mit maximaler Anzahl Pixeln
  # 	max.pix.vec <- rep(0,n.p.total)
  # 			
  # 	# Vector mit maximaler Anzahl Hotspots
  # 	max.hots.vec <- rep(0,n.p.total)
  
  # Liste, in die alle Variablen geschrieben werden, die vor dem Tracking
  #		berechnet werden
  vars.hot.list <- vector("list", n.p.total)
  
  # Vector mit Info, ob Schwarm oder nicht
  is.schwarm.vec <- rep(FALSE,n.p.total)
  
  if(n.cores>1){
    
    # 		cat("Multicore aktuell nicht möglich \n")
    
    # Ehemals for-Schleife
    help.res <- mclapply(which.packs, 
                         help.prepro,
                         track.evals = track.evals,
                         my.file = my.file,
                         win.start = win.start,
                         win.length = win.length,
                         vers = vers,
                         a.1 = a.1,
                         a.2 = a.2,
                         cut = cut,
                         maxdist = maxdist,
                         c = c,
                         signal.neg = signal.neg,
                         m.d.cs = m.d.cs,
                         n.angle = n.angle,
                         mc.cores = n.cores,
                         do.plot = do.plot,
                         do.plot.vec = do.plot.vec,
                         y.lim=y.lim,
                         df.t=df.t,
                         pfad.mult=pfad.mult,
                         frames.pack = frames.pack,
                         fastclust=fastclust,
                         speeditup=speeditup,
                         floodclust=floodclust,
                         floodCpp=floodCpp,
                         plot.hots.mult=plot.hots.mult,
                         only.watch=only.watch,
                         schwarm.find=schwarm.find,
                         regel.schwarm=regel.schwarm,
                         which.regel.schwarm=which.regel.schwarm,
                         save.jpg=save.jpg,
                         zeitstempel=zeitstempel,
                         save.movie=save.movie, 
                         delete.jpg=delete.jpg,
                         each=each,
                         wait.jpg=wait.jpg,
                         win=win
    )
  }else{
    
    # Ehemals for-Schleife
    # Variante, wenn n.cores=1. 
    #	Grund: Paket multicore nur für mac und linux
    help.res <- lapply(which.packs, 
                       help.prepro,
                       track.evals = track.evals,
                       my.file = my.file,
                       win.start = win.start,
                       win.length = win.length,
                       vers = vers,
                       a.1 = a.1,
                       a.2 = a.2,
                       cut = cut,
                       maxdist = maxdist,
                       c = c,
                       signal.neg = signal.neg,
                       m.d.cs = m.d.cs,
                       n.angle = n.angle,
                       # 								 mc.cores = n.cores,
                       do.plot = do.plot,
                       do.plot.vec = do.plot.vec,
                       y.lim=y.lim,
                       df.t=df.t,
                       pfad.mult=pfad.mult,
                       frames.pack = frames.pack,
                       fastclust=fastclust,
                       speeditup=speeditup,
                       floodclust=floodclust,
                       floodCpp=floodCpp,
                       plot.hots.mult=plot.hots.mult,
                       only.watch=only.watch,
                       schwarm.find=schwarm.find,
                       regel.schwarm=regel.schwarm,
                       which.regel.schwarm=which.regel.schwarm,
                       save.jpg=save.jpg,
                       zeitstempel=zeitstempel,
                       save.movie=save.movie, 
                       delete.jpg=delete.jpg,
                       each=each,
                       wait.jpg=wait.jpg,
                       win=win
    )
  }
  
  # Alles so abspeichern, wie es gehört
  for(p in 1:n.p.total){
    
    track.eval.out.list[[p]] <- help.res[[p]]$track.eval.out
    track.ddf.list[[p]] <- help.res[[p]]$track.ddf3
    merkmale.out.list[[p]] <- help.res[[p]]$merkmale.out
    t.grid.list[[p]] <- help.res[[p]]$t.grid
    pfad.mult.list[[p]] <- help.res[[p]]$pfad.mult
    hots.y01.list[[p]] <- help.res[[p]]$hots.y01
    dims.list[[p]] <- help.res[[p]]$dims.l
    track.evals[[p]] <- help.res[[p]]$track.ev
    # 		max.pix.vec[p] <- help.res[[p]]$max.pix
    # 		max.hots.vec[p] <- help.res[[p]]$max.hots
    vars.hot.list[[p]] <- help.res[[p]]$vars.hot
    is.schwarm.vec[p] <- help.res[[p]]$is.schwarm
  }
  
  # Kartesische Koordinaten sind immer die gleichen, deshalb nehme ich die
  #		des ersten Frames
  cart.coord <- help.res[[1]]$cart.coord
  
  # 	# Wenn im ersten Paket kein Hotspot gefunden wird ist help.res=NULL
  # 	#		und somit auch cart.coord leer, dann muss man in den nächsten suchen
  # 	ind <- 1
  # 	
  # 	while(is.null(cart.coord) & ind<=n.p.total){
  # 			ind <- ind + 1
  # 			cart.coord <- help.res[[ind]]$cart.coord
  # 		}
  # UPDATE 27.11.12: Brauche ich nicht mehr, da ich in der Funktion jetzt
  #		immer cart.coord ausgeben lasse
  
  if(save.preprocess){
    # Hart abspeichern und bei Klassifikation darauf zurückgreifen
    save(track.eval.out.list,
         track.ddf.list,
         track.evals,
         merkmale.out.list,
         t.grid.list,
         pfad.mult.list,
         hots.y01.list,
         dims.list,
         n.angle,
         cart.coord,
         # 			 max.pix.vec,
         # 			 max.hots.vec,
         vars.hot.list,
         is.schwarm.vec,
         file=paste(pfad.mult, zeitstempel,"Results_n_angle_",n.angle,"_a2_",a.2,".RData",sep="")
    )
  }
  
  out <- list(track.eval.out.list=track.eval.out.list,
              track.ddf.list=track.ddf.list,
              track.evals=track.evals,
              merkmale.out.list=merkmale.out.list,
              t.grid.list=t.grid.list,
              pfad.mult.list=pfad.mult.list,
              hots.y01.list=hots.y01.list,
              dims.list=dims.list,
              n.angle=n.angle,
              cart.coord=cart.coord,
              vars.hot.list=vars.hot.list,
              is.schwarm.vec=is.schwarm.vec)
  return(out)
}

####################################################################
# Die Funktion help.prepro() ist eine Hilfsfunktion für preprocess.data(),
#	dieser Teil wird ausgelagert, um eine Parallelisierung mittels
#	mclapply() möglich zu machen
#
# Übergeben werden muss:
#	- p:		Index der "for"-Schleife
#	- track.evals,
# 	- my.file,
# 	- win.start,
# 	- win.length,
# 	- vers,
# 	- a.1,
# 	- a.2,
# 	- cut,
# 	- maxdist,
# 	- c,
# 	- signal.neg,
# 	- m.d.cs,
# 	- n.angle
#	- df.t
#	- pfad.mult
#	- fastclust
#	- speeditup
#
# Ausgegeben wird:
#	- Liste mit den Ergebnissen 
#
####################################################################


help.prepro <- function(p,
                        track.evals,
                        my.file,
                        win.start,
                        win.length,
                        vers,
                        a.1,
                        a.2,
                        cut,
                        maxdist,
                        c,
                        signal.neg,
                        m.d.cs,
                        n.angle,
                        do.plot,
                        do.plot.vec,
                        y.lim,
                        df.t,
                        pfad.mult,
                        frames.pack,
                        fastclust,
                        speeditup,
                        floodclust,
                        floodCpp,
                        plot.hots.mult,
                        only.watch,
                        schwarm.find,
                        regel.schwarm,
                        which.regel.schwarm,
                        save.jpg,
                        zeitstempel,
                        save.movie, 
                        delete.jpg,
                        each,
                        wait.jpg,
                        win
){
  
  cat(as.character(Sys.time())," help.prepro mit Paket: \n")
  cat(paste("p =",p,"von",length(track.evals)),"\n")
  
  # Interessante Frames
  rep.time <- 1
  from <- (p-1)*frames.pack*rep.time+rep.time
  to <- p*frames.pack*rep.time
  t.grid <- seq(from,to,by=rep.time)	
  
  # ddf einlesen
  
  cat(as.character(Sys.time())," Beginn Daten einlesen \n")
  
  ddf.data <- read.ddf.allframes(file=my.file, 
                                 frames=t.grid, 
                                 win.start=win.start,
                                 win.length=win.length,
                                 vers=vers)
  
  cat(as.character(Sys.time())," Ende Daten einlesen \n")
  
  # 	# Version und win.start und win.length wird ausgelesen
  # 	data.new <- read.ddf.allframes.new(file=my.file, 
  # 																		 frames=t.grid)
  # 	
  # 	get.version(my.file)
  
  # Kartesische Koordinaten der Pixel
  xs <- ddf.data$xs[1:y.lim,]
  ys <- ddf.data$ys[1:y.lim,]
  
  # Matrix mit kartesischen Koordinaten jedes Punktes
  cart.coord <- array(dim=c(y.lim, 96, 2))
  cart.coord[,,1] <- xs
  cart.coord[,,2] <- ys
  
  # Anzahl Frames
  n.frames <- length(t.grid)
  
  # Interessante Zeitpunkte
  t.plot <- 1:n.frames
  
  # Datenarray (also ohne xs und ys und num_frames)	
  Y <- array(dim=c(y.lim, 96, n.frames))
  
  for(i in 1:n.frames){
    Y[,,i] <- ddf.data[[i+3]][1:y.lim,]
  }
  
  # Anzahl Freiheitsgrade der B-Spline-Basis für die 3 Dimensionen
  if(is.null(df.t)){
    df <- c(100, 25, round(n.frames/2))
  }else{
    if(length(df.t) == 1){
      df <- c(100, 25, df.t)	
    }else{
      df <- df.t
    }
  }
  
  #####################################################################
  # 1.) Hotspots finden mit der Funktion find.mult.hotspots()
  #####################################################################
  
  cat(as.character(Sys.time())," Beginn Hotspots finden \n")
  
  hots.mult <- find.mult.hotspots(Y=Y, 
                                  df=df, 
                                  a.1=a.1, 
                                  a.2=a.2, 
                                  cut=cut, 
                                  c=c,
                                  maxdist=maxdist,
                                  signal.neg=signal.neg,
                                  fastclust=fastclust,
                                  speeditup=speeditup,
                                  floodclust=floodclust,
                                  floodCpp=floodCpp,
                                  pack=p,
                                  plot.hots.mult=plot.hots.mult,
                                  save.jpg=save.jpg
  )
  
  cat(as.character(Sys.time())," Ende Hotspots finden \n")
  
  # Ordner erstellen, in den die jpgs kommen
  folder <- paste(pfad.mult, zeitstempel, sep="")
  
  # Falls das Video mehrere Teile hat, 
  #	brauche ich auch mehrere Unterordner
  folderp <- paste(folder,"_Teil_",p,sep="")
  
  # 	system(paste("mkdir",folder, sep=" "))#, wait=FALSE)
  # Geht offensichtlich unter windows nicht, deshalb: 
  
  if(!file.exists(sub("/","",pfad.mult))){
    
    dpm <- dir.create(pfad.mult)
    
    # Hat Erstellen des Ordners geklappt?
    if(dpm){
      cat(as.character(Sys.time())," Ordner", pfad.mult, "erstellt \n")
    }else{
      cat(as.character(Sys.time())," Ordner", pfad.mult, "nicht erstellt \n")
    }
  }else{
    cat(as.character(Sys.time())," Ordner", pfad.mult, "existiert bereits \n")
  }
  
  if(!file.exists(folder)){
    
    folder.true <- dir.create(folder)
    
    # Hat Erstellen des Ordners geklappt?
    if(folder.true){
      cat(as.character(Sys.time())," Ordner", folder, "erstellt \n")
    }else{
      cat(as.character(Sys.time())," Ordner", folder, "nicht erstellt \n")
    }
  }else{
    cat(as.character(Sys.time())," Ordner", folder, "existiert bereits \n")
  }
  
  if(p!=1){
    
    if(!file.exists(folderp)){
      
      folder.true <- dir.create(folderp)
      
      # Hat Erstellen des Ordners geklappt?
      if(folder.true){
        cat(as.character(Sys.time())," Ordner", folderp, "erstellt \n")
      }else{
        cat(as.character(Sys.time())," Ordner", folderp, "nicht erstellt \n")
      }
    }else{
      cat(as.character(Sys.time())," Ordner", folderp, "existiert bereits \n")
    }
    
    folder <- folderp
  }
  
  # Sollen die jpgs gespeichert werden, mit denen dann der Film erstellt
  #		werden kann, wo links Y und rechts Y.01 ist?
  # Mit save.movie=TRUE wird der Film dann auch gleich erstellt
  if(save.jpg){
    
    # Wo soll das.Rout file gespeichert werden?
    rout <- paste(folder,"/make_movie.Rout",sep="")
    
    Y.01 <- hots.mult$Y.01
    
    cat(as.character(Sys.time())," Beginn make_movie.RData speichern  \n")
    
    # Y, Y.01, xs und ys speichern, damit der andere Prozess es wieder 
    #		findet
    save(Y, Y.01, xs, ys, 
         file=paste(folder,"/make_movie.RData",sep=""))
    
    cat(as.character(Sys.time())," Ende make_movie.RData speichern  \n")
    
    # Aufruf, der den Prozess startet, die jpgs zu schreiben. 
    # Bei windows anders als bei Linux und Mac
    if(win){
      jpgaufruf <- paste('Rcmd BATCH --no-save --no-restore "--args ',
                         folder, 
                         save.movie, 
                         delete.jpg, 
                         each,
                         ' " make_movie.R', rout) 
      # unter windows \ statt / in der Konsole? " statt '
    }else{
      jpgaufruf <- paste("R CMD BATCH --no-save --no-restore '--args ",
                         folder, 
                         save.movie, 
                         delete.jpg, 
                         each,
                         " ' make_movie.R", rout)
    }
    
    # 		cat(as.character(Sys.time()), jpgaufruf, "  \n")
    
    system(jpgaufruf, wait=wait.jpg)
  }
  
  cat(as.character(Sys.time())," Beginn hots.merkmale \n")
  
  # Merkmale auf Hotspotebene bestimmen
  vars.hot <- hots.merkmale(hots.mult=hots.mult,
                            cart.coord=cart.coord)
  
  cat(as.character(Sys.time())," Ende hots.merkmale \n")
  
  # Auf FALSE setzen, da es immer mit ausgegeben wird
  is.schwarm <- FALSE
  
  # Falls auf Schwarm getestet werden soll
  if(schwarm.find){
    # 		out <- list(max.pix=max.pix,
    # 								max.hots=max.hots,
    # 								cart.coord=cart.coord)
    
    # Ist der Film ein Schwarmfilm?
    is.schwarm <- is.schwarm.func(vars.hot=vars.hot,
                                  regel.schwarm=regel.schwarm,
                                  which.regel.schwarm=which.regel.schwarm)
    
    # Falls es ein Schwarmfilm ist: abbrechen, sonst weiter
    if(is.schwarm){
      
      out <- list(vars.hot=vars.hot,
                  cart.coord=cart.coord,
                  is.schwarm=is.schwarm)
      
      return(out)
    }
    
  }
  
  if(plot.hots.mult){
    t.plot.hotspots.mult <- round(seq(1,n.frames, length=16))
    plot.hotspots.mult(hotspots=hots.mult,
                       pfad=paste(pfad.mult, zeitstempel,"/", sep=""),
                       t.plot=t.plot.hotspots.mult,
                       mfrow.tplot=n2mfrow(length(t.plot.hotspots.mult)),
                       mfrow.all=n2mfrow(length(t.plot.hotspots.mult)),
                       xs=xs,
                       ys=ys,
                       t.grid=c(1:n.frames),
                       plot.hots=TRUE
    )
  }
  
  
  ###########################################################################
  #	2.) Tracking
  ###########################################################################
  
  cat(as.character(Sys.time())," Beginn Tracking \n")
  
  # nur weitermachen, wenn mindestens ein Hotspot gefunden wurde
  if(sum(hots.mult$n.groups)>0){
    
    track.ddf <- track.func(hots.mult=hots.mult,
                            t.plot=t.plot,
                            # max.dist.cs=0.05, # Standardwert
                            max.dist.cs=m.d.cs,
                            rep.time=rep.time,
                            cart.coord=cart.coord)
    
    
    # 2. Teil Tracking
    track.ddf2 <- track.split.over(track.ddf=track.ddf)
    
    # 3. Teil Tracking
    track.ddf3 <- track.split(track.ddf=track.ddf2)
    
    ##################################################################
    # 2.a) Ergebnis evaluieren und Evaluation abspeichern
    ##################################################################
    
    track.eval.out <- track.eval(track.ddf3)
    
    track.ev <- track.eval.out$track.ev
    
    # 		track.evals[[p]] <- track.ev
    
    ##########################################################################
    # 2.b) Ergebnis plotten
    ##########################################################################	
    
    if(do.plot){
      filename <- paste("Tracking_",from,"_",to,"_by",rep.time,".jpg",sep="")
      
      plot.track(track.ddf=track.ddf3,
                 hots.mult=hots.mult,
                 pfad=paste(pfad.mult, zeitstempel,"/",sep=""),
                 filename=filename,
                 t.grid=t.grid,
                 cart.coord=cart.coord)
    }
    
    
    cat(as.character(Sys.time())," Ende Tracking \n")
    
    ##########################################################################		
    # 3.) Merkmale bestimmen und Vektoren plotten
    ##########################################################################
    
    cat(as.character(Sys.time())," Beginn Merkmale \n")
    
    # Merkmale bestimmen	
    merkmale.out <- merkmale.func(hotspots=hots.mult,
                                  track.ddf=track.ddf3,
                                  track.eval.out=track.eval.out,
                                  n.angle=n.angle,
                                  cart.coord=cart.coord,
                                  t.plot=t.plot,
                                  only.watch=only.watch)
    
    cat(as.character(Sys.time())," Ende Merkmale \n")
    
    # Plots
    
    if(do.plot.vec){
      
      filename <- paste("Tracking_",from,"_",to,"_by",rep.time,"_Vec_Norm.jpg",sep="")
      
      plot.vecs.func(merkmale.out=merkmale.out,
                     track.ddf=track.ddf3,
                     track.eval.out=track.eval.out,
                     n.angle=n.angle, 
                     objects=NULL,
                     t.grid=t.grid, 
                     norm=TRUE,
                     do.save=TRUE,
                     pfad=paste(pfad.mult, zeitstempel,"/", sep=""),
                     filename=filename)
      
      filename <- paste("Tracking_",from,"_",to,"_by",rep.time,"_Vec_True.jpg",sep="")
      
      plot.vecs.func(merkmale.out=merkmale.out,
                     track.ddf=track.ddf3,
                     track.eval.out=track.eval.out,
                     n.angle=n.angle, 
                     objects=NULL,
                     t.grid=t.grid, 
                     norm=FALSE,
                     do.save=TRUE,
                     pfad=paste(pfad.mult, zeitstempel,"/", sep=""),
                     filename=filename)
    }
    
    # Alles abspeichern
    out <- list(track.eval.out = track.eval.out,
                track.ddf3 = track.ddf3,
                merkmale.out = merkmale.out,
                t.grid = t.grid,
                pfad.mult = paste(pfad.mult, zeitstempel,"/", sep=""),
                hots.y01 = hots.mult$hots.y01,
                dims.l = hots.mult$dims.l,
                track.ev = track.ev,
                cart.coord = cart.coord,
                # 								max.pix=max.pix,
                # 								max.hots=max.hots,
                vars.hot=vars.hot,
                is.schwarm=is.schwarm)
    
    return(out)
    
  }else{
    
    # Sonst nur cart.coord ausgeben
    # 		out <- list(max.pix=max.pix,
    # 								max.hots=max.hots,
    # 								cart.coord = cart.coord)
    out <- list(vars.hot=vars.hot,
                cart.coord = cart.coord,
                is.schwarm=is.schwarm)
    
    return(out)
    
  }
}

####################################################################
# Die Funktion plot.prepro.hot() soll die Ergebnisse des Preprocessings
#	plotten, also die Tracks:
#		Großes Schachbrett mit #Zeilen = #Objekte, #Spalten = max.length
#
# Übergeben werden muss:
#	- hot:					TRUE => Hotspots, FALSE => Vektoren
#	- fisch:					Name der Fischart
#	- pfad:					Dateipfad wohin geplottet werden soll
#	- sum.prepro:			Ergebnis aus summary.func
#	- min.length:			Minimale Länge der zu plottenden Tracks
#	- cart.coord:			Kartesische Koordinaten der Pixel
#	- hots.y01.list:		hots.y01.list
#	- dims.list:			dims.list
#	- track.ddf.list: 	track.ddf.list
#	- merkmale.out.list:	merkmale.out.list
#	- res:					Auflösung, default=100
#	- cex:					Grafikparameter, default=1
#	- filename:				Name der Datei, optional
#	- max.row:				Maximale Anzahl Zeilen pro Plot, default=50
#	- n.cores:				Anzahl Rechenkerne
#	- a.2:					Für Dateinamen
#
# Ausgegeben wird:
#	- nichts, nur plot
#
####################################################################

plot.prepro <- function(hot,
                        fisch,
                        pfad,
                        sum.prepro,
                        min.length,
                        cart.coord,
                        hots.y01.list,
                        dims.list,
                        track.ddf.list,
                        merkmale.out.list,
                        res=200,
                        cex=1,
                        filename=NULL,
                        max.row=50,
                        n.cores=1,
                        a.2){
  
  # Gekürztes sum.prepro, mindestens so lang wie min.length und unproble-
  #	matisch.
  sum.precut1 <- sum.prepro[which(sum.prepro$n.hots>=min.length),]
  sum.precut <- sum.precut1[which(sum.precut1[,4]==TRUE),]
  
  # Anzahl der Objekte, die mindestens 5 Hotspots haben
  anz.obj <- nrow(sum.precut)
  max.length <- max(sum.prepro$n.hots) #sum.precut?
  
  # 	# Anzahl Zeilen und Spalten bestimmen
  # 	nRow <- min(100,anz.obj)
  # 	nCol <- max.length
  # 	
  # 	# Grafikparameter
  # 	height <- nRow * 1
  # 	width <-  nCol * 1
  
  if(hot){
    # Kartesische Koordinaten
    xs <- cart.coord[,,1]
    ys <- cart.coord[,,2]
  }
  
  # Dateiname
  if(is.null(filename)){
    filename.basis  <- paste(pfad,"Tracks_min_",min.length, sep="")
  }else{
    filename.basis  <- paste(pfad,filename, sep="")	
  }
  
  # Anzahl Zeilen in letztem Plot
  rest <- anz.obj%%max.row # 
  
  # Anzahl Plots 
  if(rest==0){
    anz.plots <- (anz.obj - rest) /max.row
  }else{
    anz.plots <- (anz.obj - rest) /max.row + 1
  }
  
  
  
  if(hot){
    
    lapply(1:anz.plots,
           plot.prepro.hot.help,
           anz.plots=anz.plots,
           max.row=max.row,
           filename.basis=filename.basis,
           anz.obj=anz.obj,
           fisch=fisch,
           sum.precut=sum.precut,
           hots.y01.list=hots.y01.list,
           dims.list=dims.list,
           track.ddf.list=track.ddf.list,
           res=res,
           xs=xs,
           ys=ys,
           #					mc.cores=n.cores,
           a.2=a.2)
    
  }else{
    
    lapply(1:anz.plots,
           plot.prepro.vec.help,
           anz.plots=anz.plots,
           max.row=max.row,
           filename.basis=filename.basis,
           anz.obj=anz.obj,
           fisch=fisch,
           sum.precut=sum.precut,
           # 					hots.y01.list=hots.y01.list,
           # 					dims.list=dims.list,
           # 					track.ddf.list=track.ddf.list,
           merkmale.out.list=merkmale.out.list,
           res=res,
           #					mc.cores=n.cores,
           a.2=a.2)
  }
  
  # 	for(pl.i in 1:anz.plots){
  # 		
  # 		i.range <- (1 + 100*(pl.i-1)):min(pl.i*100,anz.obj)
  # 		
  # 		if(hot){
  # 			filename_pl <- paste(filename.basis,"_hots_",pl.i,".jpg",sep="")
  # 			
  # 			plot.prepro.hot.help(filename_pl=filename_pl,
  # 										i.range=i.range,
  # 										fisch=fisch,
  # 										sum.precut=sum.precut,
  # 										hots.y01.list=hots.y01.list,
  # 										dims.list=dims.list,
  # 										track.ddf.list=track.ddf.list,
  # 										res=res,
  # 										xs=xs,
  # 										ys=ys)
  # 
  # 		}else{
  # 			filename_pl <- paste(filename.basis,"_vecs_",pl.i,".jpg",sep="")
  # 			
  # 			plot.prepro.vec.help(filename_pl=filename_pl,
  # 										i.range=i.range,
  # 										fisch=fisch,
  # 										sum.precut=sum.precut,
  # # 										hots.y01.list=hots.y01.list,
  # # 										dims.list=dims.list,
  # # 										track.ddf.list=track.ddf.list,
  # 										merkmale.out.list=merkmale.out.list,
  # 										res=res)
  # 		}
  # 	}
}

####################################################################
# Die Funktion plot.prepro.hot.help() ist eine Hilfsfunktion für 
#	plot.prepro(), damit die Plots maximal max.row Zeilen haben. 
#
# Übergeben werden muss:
#	- filename.basis:	Dateiname, Basis
#	- max.row:			Maximale Anzahl Zeilen
#	- fisch:				Name der Fischart
#	- sum.precut:		Ergebnis aus summary.func
#	- hots.y01.list:	hots.y01.list
#	- dims.list:		dims.list
#	- track.ddf.list: track.ddf.list
#	- res:				Auflösung des jpegs
#	- cex:				Grafikparameter, default=1
#
#	ACHTUNG, HIER FEHLT NOCH WAS
#
# Ausgegeben wird:
#	- nichts, nur plot
#
####################################################################

plot.prepro.hot.help <- function(pl.i,
                                 anz.plots,
                                 max.row,
                                 filename.basis,
                                 anz.obj,
                                 fisch,
                                 sum.precut,
                                 hots.y01.list,
                                 dims.list,
                                 track.ddf.list,
                                 xs,
                                 ys,
                                 res,
                                 cex=1,
                                 a.2){
  
  cat(paste("Plot Hotspots", pl.i,"von",anz.plots), "\n")
  
  i.range <- (1 + max.row*(pl.i-1)):min(pl.i*max.row,anz.obj)
  
  filename_pl <- paste(filename.basis,"_hots_a2_",a.2,"_",pl.i,".jpg",sep="")
  
  height <- nRow <- length(i.range)
  
  width <- nCol <- min(100,max(sum.precut$n.hots[i.range]))
  
  jpeg(filename_pl,
       width=width,
       height=height, 
       units="in", 
       res=res, 
       quality=100)
  
  # Plot initialisieren
  par(mfrow=c(nRow, nCol), mar=c(1.7,1.8,1.7,1.2)+0.2)
  
  for(i in i.range){	
    
    # 		print(i)
    
    # Objektname
    id <- sum.precut[i,1]
    
    # Paket des Films, in dem das Objekt auftauchte
    p <- sum.precut[i,2]
    
    # Objektname innerhalb des Pakets
    obj.in.p <- sum.precut[i,3]
    
    # Hotspots in diesem Paket
    hots.y01 <- hots.y01.list[[p]]
    
    # Dimensionen der Hotspots in diesem Paket
    dims.l <- dims.list[[p]]
    
    # Track.ddf zu diesem Paket
    track.ddf <- track.ddf.list[[p]]
    
    # Tracker und Matcher zu diesem Paket
    tracker <- track.ddf$tracker
    matcher <- track.ddf$matcher		
    
    # IDs innerhalb des Pakets p, die zu diesem Objekt gehören
    ids <- which(is.element(tracker,obj.in.p))
    
    # Da ich nur die unproblematischen betrachte, kann ich mir das 
    #	komplizierte hier wohl sparen
{
      # 		ids <- NULL
      # 		
      # 		# n.objects enthält die Anzahl Objekte zu diesem Hotspot
      # 		# Achtung, ist was anderes als Anzahl Hotpots zu einem Objekt
      # 		
      # 		n.objects <- NULL
      # 		
      # 		for(k in 1:length(tracker)){
      # 			if(length(intersect(tracker[[k]],obj.in.p))>0){
      # 				ids <- union(ids, k)	
      # 				n.objects[k] <- length(tracker[[k]])
      # 			}
      # 		}
    }

# Zeitpunkte der Hotspots
ts <- matcher[ids,2]
# IDs des Hotspots innerhalb des Zeitpunkts t
i.ts <- matcher[ids,3]

for(j in 1:nCol){
  
  # Wenn Objekt id j oder mehr Hotspots hat, kann noch ein Hotspot
  #	geplottet werden, sonst leerer Plot
  if(sum.precut[i,5]>=j){
    
    # Zeitpunkt des Hotspots in Paket p
    t <- ts[j]
    
    # Nummer des Hotspots in Paket p zum Zeitpunkt t
    i.t <- i.ts[j]
    
    # Dimensionen dieses Hotspots
    dims <- dims.l[[t]][[i.t]]
    
    # Plot nur, wenn es was zu plotten gibt
    if(!is.null(dims)){
      
      # x und y Koordinaten der zu plottenden Punkte
      xs.plot <- -xs[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
      ys.plot <- ys[dims[1,1]:dims[2,1],dims[1,2]:dims[2,2]]
      
      # Titel für den Plot
      # 					main <- paste(fisch,id, sep=",")
      main <- paste("ID",id,",t",t, sep="")
      
      # Plot
      plot.ddf(hots.y01[[t]][[i.t]], bw=0, xs=xs.plot, ys=ys.plot, 
               bin=TRUE, cex=cex, main=main, col=1)
    }else{
      plot(0,0,col="white", 
           main="",
           axes=FALSE, xlab="", ylab="")
      box()
    }
  }else{
    plot(0,0,col="white", 
         main="",
         axes=FALSE, xlab="", ylab="")
    box()
  }
}
  }
dev.off()
}

####################################################################
# Die Funktion plot.prepro.vec.help() ist eine Hilfsfunktion für 
#	plot.prepro(), damit die Plots maximal 100 Zeilen haben. 
#
# Übergeben werden muss:
#	- filename_pl:			Dateiname
#	- i.range:				Range der IDs
#	- fisch:					Name der Fischart
#	- sum.precut:			Ergebnis aus summary.func
#	- hots.y01.list:		hots.y01.list
#	- dims.list:			dims.list
#	- track.ddf.list: 	track.ddf.list
#	- merkmale.out.list:	merkmale.out.list
#	- res:					Auflösung des jpegs
#	- cex:					Grafikparameter, default=1
#
#
#	ACHTUNG, HIER FEHLT NOCH WAS
#
#
# Ausgegeben wird:
#	- nichts, nur plot
#
####################################################################

plot.prepro.vec.help <- function(pl.i,
                                 anz.plots,
                                 max.row,
                                 filename.basis,
                                 anz.obj,
                                 fisch,
                                 sum.precut,
                                 # 									 		min.length,
                                 # 									 		cart.coord,
                                 # 									 		hots.y01.list,
                                 # 									 		dims.list,
                                 # 									 		track.ddf.list,
                                 merkmale.out.list,
                                 res,
                                 cex=1,
                                 filename=NULL,
                                 a.2,
                                 n.angle=NULL){
  
  cat(paste("Plot Vektorensterne", pl.i,"von",anz.plots), "\n")
  
  i.range <- (1 + max.row*(pl.i-1)):min(pl.i*max.row,anz.obj)
  
  filename_pl <- paste(filename.basis,"_vecs_a2_",a.2,"_",pl.i,".jpg",sep="")
  
  height <- nRow <- length(i.range)
  
  width <- nCol <- min(100,max(sum.precut$n.hots[i.range]))
  
  jpeg(filename_pl,
       width=width,
       height=height, 
       units="in", 
       res=res, 
       quality=100)
  
  # Plot initialisieren
  par(mfrow=c(nRow, nCol), mar=c(1.7,1.8,1.7,1.2)+0.2)
  
  for(i in i.range){	
    
    # 		print(i)
    
    # Objektname
    id <- sum.precut[i,1]
    
    # Paket des Films, in dem das Objekt auftauchte
    p <- sum.precut[i,2]
    
    # Objektname innerhalb des Pakets
    obj.in.p <- sum.precut[i,3]
    
    # Merkmale
    merkmale <- merkmale.out.list[[p]]$merkmale.list[[obj.in.p]]
    
    x.s <- vector(length=4*n.angle)
    y.s <- vector(length=4*n.angle)		
    
    for(j in 1:nCol){
      
      # Wenn Objekt id j oder mehr Hotspots hat, kann noch ein Stern
      #	geplottet werden, sonst leerer Plot
      if(sum.precut[i,5]>=j){
        
        if(!is.null(merkmale[12,j])){
          
          # Normieren, also "Schnauze rechts" und deren Länge = 1				
          # X0.Grad auf 1 setzen
          merkmale[12,j] <- 1
          
          # Winkel (nicht die wahren, sondern die, für die sie sich 
          #	ausgeben) im Gradmaß
          alphas <- c(90/n.angle * 0:((4*n.angle)-1))
          
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
          
          # Titel für den Plot
          main <- paste(fisch,id, sep=",")
          
          # Stern plotten					
          plot(0,0, xlim=xlim, ylim=ylim, xlab="xs", ylab="ys",
               main=main)
          points(x.s, y.s, pch=16)
          for(p in 0:((4*n.angle)-1)){
            lines(c(x.s[p+1],0), c(y.s[p+1],0))
          }
        }else{
          plot(0,0,col="white", 
               main="",
               axes=FALSE, xlab="", ylab="")
          box()
        }
      }
      else{
        plot(0,0,col="white", 
             main="",
             axes=FALSE, xlab="", ylab="")
        box()
      }
    }
  }
  dev.off()
}

###########################################################################
# Die Funktion summary.func() soll die Merkmale aller Filmschnipsel
#	zusammen fassen
#
# Übergeben werden muss:
#
#	- track.ddf.list
#	- track.eval.out.list
#	- merkmale.out.list
#	- n.angle
#
# Ausgegeben wird:
#
#	- klass.out: Matrix mit:
#		- Eindeutige ID für die Objekte
#		- Zeitraum des Objekts
#		- Nummer innerhalb des Zeitraums
#		- Anzahl Tracks
#		- Unproblematisch? T/F
#		- Merkmale: 
#			- Mean der Uhrzeiger
#			- Varianz der Uhrzeiger
#
###########################################################################

#' Summary of preprocessing
#' 
#' @param track.ddf.list Output of \code{\link{preprocess.data}}
#' @param track.eval.out.list Output of \code{\link{preprocess.data}}
#' @param merkmale.out.list Output of \code{\link{preprocess.data}}
#' @param t.grid.list Output of \code{\link{preprocess.data}}
#' @param n.angle Number of watch hands per quarter, i.e., total number of 
#'  watch hands is \code{4 x n.angle}
#' @return Matrix summarizing results of preprocessing for each object
#' @export
summary.func <- function(track.ddf.list,
                         track.eval.out.list,
                         merkmale.out.list,
                         t.grid.list,
                         n.angle){
  
  # N = # Uhrzeiger
  N <- n.angle*4
  
  # Anzahl Objekte insgesamt:
  anz.obj.total <- 0
  
  for(p in 1:length(t.grid.list)){
    
    if(!is.null(track.ddf.list[[p]]$anz.obj)){
      anz.obj.total <- anz.obj.total + track.ddf.list[[p]]$anz.obj
    }
  }
  
  # Track.ev.long sind alle track.evs untereinander
  track.ev.long <- NULL
  for(p in 1:length(t.grid.list)){
    
    if(!is.null(track.ddf.list[[p]]$anz.obj)){
      track.ev.long <- rbind(track.ev.long,
                             cbind( rep(p,track.ddf.list[[p]]$anz.obj), 
                                    track.eval.out.list[[p]]$track.ev))
      
    }
  }
  
  # klass.out soll die Endmatrix sein
  klass.out <- cbind(1:anz.obj.total, track.ev.long)
  
  colnames(klass.out)[1:6] <- c("ObjectID", "p", "objectInP", "noprob", 
                                "n.hots", "problem")
  
  # Merkmale dazu
  # 	klass.out[7:(N*2+6] <- 0
  anz.col <- (N*2+6+(N*(N-1)/2) + 4)
  klass.out[7:anz.col] <- 0
  colnames(klass.out)[7:(N+6)] <- paste("Mean_",90/n.angle * 0:((4*n.angle)-1), "_Grad", sep="")
  colnames(klass.out)[(N+7):(N*2+6)] <- paste("Stdev_",90/n.angle * 0:((4*n.angle)-1), "_Grad", sep="")
  namen <- paste(rep(1:N,N),"_",rep(1:N,each=N),sep="")
  S <- matrix(namen, ncol=N)
  w <- which(upper.tri(S)==TRUE, arr.ind=TRUE)
  colnames(klass.out)[(N*2+7):(anz.col-4)] <- S[w]
  
  colnames(klass.out)[(anz.col-3):anz.col] <- c("alpha.mean", "alpha.sd", 
                                                "rad.mean", "rad.sd")
  
  vel <- vector(length=anz.obj.total)
  
  # Alle Objekte durchgehen
  for(i in 1:anz.obj.total){
    
    # 		print(i)
    
    # Nur rechnen, wenn kein Problem und mehr als ein Hotspots
    if(klass.out[i,4]==TRUE & klass.out[i,5]>1){
      
      # Zeitraum
      p <- klass.out[i,2]
      
      # Objekt innerhalb Zeitraum
      j <- klass.out[i,3]
      
      merks <- matrix(merkmale.out.list[[p]][[1]][[j]][12:(11+4*n.angle),], 
                      nrow=N)
      
      # Für alle Hotspots dieses Objekts...
      for(k in 1:ncol(merks)){
        # ... Uhrzeiger zurückrechnen auf echte Werte
        if(merks[1,k]>0){
          merks[-1,k] <- merks[-1,k]*merks[1,k]
        }
      }
      
      # Mean für jeden Uhrzeiger
      means <- apply(X=merks, MARGIN=1, mean)
      klass.out[i,7:(N+6)] <- means
      
      if(ncol(merks)>1){
        # Standardabweichungen für jeden Uhrzeiger
        stdevs <- apply(X=merks, MARGIN=1, sd)
        klass.out[i,(N+7):(N*2+6)] <- stdevs
        
        # Korrelationen der Uhrzeiger
        X <- cor(t(merks))				
        klass.out[i,(N*2+7):(anz.col-4)] <- X[w]
      }
      
      # Winkel und Radius, mean und sd
      klass.out[i,(anz.col-3):anz.col] <- merkmale.out.list[[p]]$dir.vel.list[[j]]
      
      
      # Geschwindigkeit
      vel[i] <- merkmale.out.list[[p]]$v.vec[j]
      
    }
  }
  
  klass.out$vel <- vel
  
  return(klass.out)
  
}

###########################################################################
# Die Funktion hots.merkmale() bestimmt die Merkmale auf Hotspotsebene,
#		bevor das Tracking durchgeführt wird. Die Hoffnung ist, dass man 
#		mit diesen Merkmalen Schwärme erkennen kann
#
# Übergeben werden muss:
#
#	- hots.mult
#	- cart.coord
#
# Ausgegeben wird:
#
#	- vars.hot: Data.frame mit den Merkmalen (nur eine Zeile)
#
###########################################################################

hots.merkmale <- function(hots.mult,
                          cart.coord){
  
  # Anzahl Hotspots pro Zeitpunkt
  # 		sort(hots.mult$n.groups)
  max.hots <- max(hots.mult$n.groups, na.rm=TRUE)
  
  # 95 % - Quantil
  q95.hots <- quantile(hots.mult$n.groups, 0.95, na.rm=TRUE)
  
  # 	# Mean der 20 % größten (10 %)
  # 	# sort(hots.mult$n.groups)[(round(length(hots.mult$n.groups)*0.85):length(hots.mult$n.groups))]
  mean.q90.hots <- mean(sort(hots.mult$n.groups)[(round(length(hots.mult$n.groups)*0.9):length(hots.mult$n.groups))])
  
  
  # Anzahl Pixel pro Zeitpunkt
  # sort(apply(X=hots.mult$Y.01, MARGIN=3, sum, na.rm=TRUE))
  # 	max.pix <- max(apply(X=hots.mult$Y.01, MARGIN=3, sum, na.rm=TRUE))
  
  # Statt Anzahl Pixel besser die Fläche:
  
  ##################################
  # Fläche eines Pixels bestimmen	
  ##################################
  
  Y.01 <- hots.mult$Y.01
  
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
  
  ################################# 
  # Fläche fertig
  #################################
  
  # Gesamtfläche zu jedem einzelnen Zeitpunkt
  hots.flaeche <- Y.01* fl.pix
  
  # Maximale Fläche pro Film
  max.flaeche <- max(apply(X=hots.flaeche, MARGIN=3, sum, na.rm=TRUE))
  # 	sort(apply(X=hots.flaeche, MARGIN=3, sum, na.rm=TRUE))
  
  ##################################
  # Fläche jedes Hotspots bestimmen, dann Mean, Median, Varianz der
  #		Fläche pro Zeitpunkt
  ##################################
  
  n.frames <- length(hots.mult$n.groups)
  
  mean.h <- rep(0, n.frames)
  med.h <- rep(0, n.frames)
  sd.h <- rep(0, n.frames)
  
  for(i in 1:n.frames){
    
    if(hots.mult$n.groups[i]>0){
      
      # Vektor, in die Fläche jedes Hotspots zu diesem Zeitpunkt geschrieben
      #		wird
      fl.hot <- hots.mult$n.groups[i]
      
      for(j in 1:hots.mult$n.groups[i]){
        
        fl.hot[j] <- sum(hots.mult$hots[[i]][[j]]$y01.all*fl.pix)
        
      }
      
      mean.h[i] <- mean(fl.hot, na.rm=TRUE)
      med.h[i] <- median(fl.hot, na.rm=TRUE)
      sd.h[i] <- sd(fl.hot, na.rm=TRUE)
    }
  }
  
  # 	sort(mean.h)
  # 	sort(med.h)
  # 	sort(sd.h)
  
  # Daraus weitere Merkmale bestimmen (Anzahl erstmal egal, ich schaue
  #		nachher bei der Klassifikation, was mir was bringt und was nicht)
  max.mean.h <- max(mean.h, na.rm=TRUE)
  max.med.h <- max(med.h, na.rm=TRUE)
  max.sd.h <- max(sd.h, na.rm=TRUE)
  # 	med.mean.h <- median(mean.h, na.rm=TRUE)
  # 	med.med.h <- median(med.h, na.rm=TRUE)
  # 	med.sd.h <- median(sd.h, na.rm=TRUE)
  # 	q75.mean.h <- quantile(mean.h, 0.75, na.rm=TRUE)
  # 	q75.med.h <- quantile(med.h, 0.75, na.rm=TRUE)
  # 	q75.sd.h <- quantile(sd.h, 0.75, na.rm=TRUE)
  # 	if(length(which(sd.h>0)) == 0){
  # 		min.sd.h <- 0
  # 	}else{
  # 		min.sd.h <- min(sd.h[which(sd.h>0)], na.rm=TRUE)
  # 	}
  
  vars.hot <- round(data.frame(max.hots=max.hots,
                               q95.hots=q95.hots,
                               mean.q90.hots=mean.q90.hots,
                               max.flaeche=max.flaeche,
                               max.mean.h=max.mean.h,
                               max.med.h=max.med.h,
                               max.sd.h=max.sd.h
                               # 															 ,med.mean.h=med.mean.h,
                               # 															 med.med.h=med.med.h,
                               # 															 med.sd.h=med.sd.h,
                               # 															 q75.mean.h=q75.mean.h,
                               # 															 q75.med.h=q75.med.h,
                               # 															 q75.sd.h=q75.sd.h,
                               # 															 min.sd.h=min.sd.h
  ),4)
  
  return(vars.hot)
}


###########################################################################
# Die Funktion is.schwarm.func() entscheidet, ob der analysierte Film
#		ein Schwarmfilm ist, oder nicht
#
# Übergeben werden muss:
#
#	- vars.hot:				output von hots.merkmale
#	- regel.schwarm:	Dateipfad der Regel zur Schwarmbestimmung
#	- which.regel.schwarm:	Welche Regel? (lda oder qda)
#
# Ausgegeben wird:
#
#	- is.schwarm:			TRUE/FALSE ob Schwarm oder nicht
#
###########################################################################

is.schwarm.func <- function(vars.hot,	
                            regel.schwarm,
                            which.regel.schwarm,
                            schwarm.qda=NULL,
                            schwarm.lda=NULL){
  
  # Regel laden
  load(regel.schwarm)
  
  # Welche Variablen wurden benutzt, um die Regel zu bestimmen?
  col.names <- attr(schwarm.qda$terms,"term.labels")
  
  # Nur die Variablen auswählen, auf denen die Regel bestimmt wurde
  vars.hot.class <- vars.hot[,c(col.names)]
  
  if(which.regel.schwarm=="qda"){
    is.schwarm <- predict(schwarm.qda, newdata=vars.hot.class)$class==1
  }else{
    is.schwarm <- predict(schwarm.lda, newdata=vars.hot.class)$class==1
  }
  
  return(is.schwarm)
}

###########################################################################
# Die Funktion movie.function() schreibt die jpgs, macht den Film und 
#	löscht die jpgs wieder
#
# Übergeben werden muss:
#
#	- 
#
# Ausgegeben wird:
#
#	- nichts
#
###########################################################################


movie.function <- function(Y,
                           Y.01,
                           xs,
                           ys,
                           save.movie=FALSE,
                           delete.jpg=FALSE,
                           folder,
                           t.grid=NULL,
                           each=1,
                           cex.main,
                           cex.axis,
                           y.lim=512,
                           main_Y="Raw signal",
                           main_Y01="Hotspots",
                           movie_name="movie",
                           movie_type=".mp4",
                           mar=c(1.7,1.8,1.7,0.2)+0.2,
                           wait=FALSE){
  
  if(is.null(t.grid)){
    
    # Welche Bilder?
    t.grid <- seq(1, dim(Y)[3], by=each)
    
  }
  
  # Wie viele Bilder werden es?
  #   => Leere Datei mit Namen endxxx.txt speichern, wobei xxx Nummer des
  #       letzten Bildes ist
  
  maxnum <- length(t.grid)
  save(maxnum, file = paste(folder,"/end",maxnum,".bin", sep=""))
  
  t <- 0
  
  for(i in t.grid){
    
    t <- t+1
    
    cat("i=",i,", t=",t,"\n")
    
    filename <- paste(folder,"/",t,".jpg", sep="") 	
    
    jpeg(filename, width=960, height=480, quality=100) 
    
    par(mfrow=c(1,2), mar=mar)
    plot.ddf(Y[,,i], bw=1, xs=-xs, ys=ys, 
             main=paste(main_Y,", t=",i,sep=""), 
             cex.main=cex.main,
             cex.axis=cex.axis,
             xlab="",ylab="")
    lines(as.vector(xs[1,]),as.vector(ys[1,]),lwd=3)
    lines(as.vector(xs[y.lim,]),as.vector(ys[y.lim,]),lwd=3)
    lines(as.vector(xs[c(1,y.lim),1]),as.vector(ys[c(1,y.lim),1]),lwd=3)
    lines(as.vector(xs[c(1,y.lim),96]),as.vector(ys[c(1,y.lim),96]),lwd=3)
    
    # 	plot.ddf(Y.hat.zentr[,,i], bw=1, xs=-xs, ys=ys, 
    # 					 main=paste("Zentriertes Signal, t=",i,sep=""), 
    # 					 cex.main=3)
    # 	lines(as.vector(xs[1,]),as.vector(ys[1,]),lwd=3)
    # 	lines(as.vector(xs[512,]),as.vector(ys[512,]),lwd=3)
    # 	lines(as.vector(xs[c(1,512),1]),as.vector(ys[c(1,512),1]),lwd=3)
    # 	lines(as.vector(xs[c(1,512),96]),as.vector(ys[c(1,512),96]),lwd=3)
    # 	
    # 	plot.ddf(Y.hat.cut[,,i], bw=1, xs=-xs, ys=ys, 
    # 					 main=paste("Abgeschnittenes Signal, t=",i,sep=""), 
    # 					 cex.main=3)
    # 	lines(as.vector(xs[1,]),as.vector(ys[1,]),lwd=3)
    # 	lines(as.vector(xs[512,]),as.vector(ys[512,]),lwd=3)
    # 	lines(as.vector(xs[c(1,512),1]),as.vector(ys[c(1,512),1]),lwd=3)
    # 	lines(as.vector(xs[c(1,512),96]),as.vector(ys[c(1,512),96]),lwd=3)
    
    if(all(is.na(Y.01[,,i]))){
      plot.ddf(1,
               main=paste(main_Y01,", t=",i,sep=""),
               cex.main=cex.main,
               cex.axis=cex.axis,xlab="",ylab="")
      # 		lines(as.vector(xs[1,]),as.vector(ys[1,]),lwd=3)
      # 		lines(as.vector(xs[512,]),as.vector(ys[512,]),lwd=3)
      # 		lines(as.vector(xs[c(1,512),1]),as.vector(ys[c(1,512),1]),lwd=3)
      # 		lines(as.vector(xs[c(1,512),96]),as.vector(ys[c(1,512),96]),lwd=3)
    }else{	
      plot.ddf(Y.01[,,i], bw=1, xs=-xs, ys=ys, 
               main=paste(main_Y01,", t=",i,sep=""), 
               cex.main=cex.main,
               cex.axis=cex.axis,xlab="",ylab="")
      lines(as.vector(xs[1,]),as.vector(ys[1,]),lwd=3)
      lines(as.vector(xs[y.lim,]),as.vector(ys[y.lim,]),lwd=3)
      lines(as.vector(xs[c(1,y.lim),1]),as.vector(ys[c(1,y.lim),1]),lwd=3)
      lines(as.vector(xs[c(1,y.lim),96]),as.vector(ys[c(1,y.lim),96]),lwd=3)
    }
    dev.off()
  }
  
  if(save.movie){
    
    # Alten Film löschen, falls vorhanden
    # 		system(paste("rm ", folder, "/movie.mp4",sep=""))
    
    file.remove(paste0(folder, "/",movie_name,movie_type))
    
    # Neuen Film schreiben
    framerate <- max(1,round(10/each))
    
    system(paste("ffmpeg -r ",framerate," -i ",folder,"/%d.jpg -b 2M ",folder,"/",movie_name,movie_type, sep=""),
           wait=wait)
    # -b gibt die Qualität an (bitrate)
  }
  
  if(delete.jpg){
    for(t in 1:length(t.grid)){
      
      filename <- paste(folder,"/",t,".jpg", sep="") 	
      
      system(paste("rm ",filename, sep=""))
    }
  }
  
}