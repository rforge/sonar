###########################################################################
# Fischprojekt
# Funktion, die die gesamte Analyse durchfuehrt, also Preprocessing mit
#		Hotspots finden und Tracking, Daten fuer Klassifikation vorbereiten, 
#		klassifizieren und Output erstellen.
# Datum der ersten Version: 25.11.13
# Letzte Aenderung: 25.11.13
# Autor: Ludwig Bothmann
###########################################################################

####################################################################
# Die Funktion analyze.ddf() analysiert einen ddf Film komplett, also
# 	Preprocessing mit
#		Hotspots finden und Tracking, Daten fuer Klassifikation vorbereiten, 
#		klassifizieren und Output erstellen.
#
# Uebergeben werden muss:
#	- ...
#
# Ausgegeben wird:
#	Nichts, alles wird direkt hart abgespeichert
#
####################################################################

#' Analyze complete .ddf-file
#'
#' @param my.file Filename of .ddf file to be analyzed
#' @param win.start Start of the sonar window, possibly extracted before via 
#'  \code{\link{get.version}}
#' @param win.length Length of the sonar window, possibly extracted before via 
#'  \code{\link{get.version}}
#' @param vers Version of the .ddf file, possibly extracted before via 
#'  \code{\link{get.version}}
#' @param max.frame Total number of frames \code{max.frame} of the given video, 
#'  possibly extracted before via \code{\link{get.version}} 
#' @param y.lim Number of pixels to analyze in each beam, counted from the 
#'  camera
#' @param a.1 Threshold for centered data
#' @param a.2 Threshold for uncentered data, default \code{0} means that no thresholding is done
#' @param cut Minimal size of cluster in number of pixels
#' @param m.d.cs Maximal distance for which two hotspots can be assigned the 
#'  same tracking number
#' @param pfad.mult Base folder for resulting plots and output
#' @param n.angle Number of watch hands per quarter, i.e., total number of watch hands is \code{4 x n.angle}
#' @param signal.neg Is the signal negativ? Default is \code{FALSE}
#' @param do.plot \code{TRUE}: Plots are saved, default is \code{FALSE}
#' @param do.plot.vec \code{TRUE}: Plots are saved, default is \code{FALSE}
#' @param do.plot.hots \code{TRUE}: Plots are saved, default is \code{FALSE}
#' @param do.plot.hots.vecs \code{TRUE}: Plots are saved, default is \code{FALSE}
#' @param plot.hots.mult \code{TRUE}: Plots are saved, default is \code{FALSE}
#' @param n.cores Number of cores if parallelization is needed, default is \code{1}
#' @param df.t Number of degrees of freedom for the splines. Default \code{NULL} 
#'  results in \code{c(100, 25, round(n.frames/2))} where \code{n.frames} is the 
#'  total number of frames of the analyzed video
#' @param frames.pack Number of frames to be analyzed at once, possibly the 
#'  total number of frames \code{max.frame} of the given video is 
#'  extracted before via \code{\link{get.version}}
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
#' @param only.watch \code{TRUE}: Use watch hands instead of trapezoid variables - 
#'  strongly recommended, default is \code{TRUE}
#' @param schwarm.find \code{TRUE}: A detection of shoal of fish is carried out, 
#'  default is \code{FALSE}
#' @param regel.schwarm File name of classification rule for shoal of fish
#' @param which.regel.schwarm Type of classification rule for shoal of fish, i.e., lda, qda...
#' @param save.preprocess \code{TRUE}: Results of preprocessing are saved, default is \code{FALSE}
#' @param zeitstempel time stamp of .ddf file
#' @param min.length Minimal length of tracks. Tracks with less hotspots are deleted.
#' @param fisch Name of fish for plots
#' @param max.row Maximal number of rows in plot of hotspots
#' @param center \code{TRUE}: Center variables on hotspot level, default is \code{FALSE}
#' @param long \code{TRUE}: Average of all hotspot variables are used, default is \code{TRUE}
#' @param file.regel File name of classification rule
#' @param klass.regel Type of classification rule for fish species, i.e., lda, qda...
#' @param sdcorr \code{TRUE}: Standard deviations and correlations of variables 
#'  on hotspot level are used, default is \code{FALSE}
#' @param FUN Function for gathering hotspot information on object level, default is \code{mean}
#' @param FUN.ow Function for gathering hotspot information on object level 
#'  when not using the trapezoid variables, default is \code{mean}
#' @param save.jpg \code{TRUE}: jpgs for the creating of the video are saved, 
#'  default is \code{TRUE}
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
#' @return As output, a table is generated which lists the found objects with 
#'  predicted class and computed features at \code{pfad.mult}
#' @export
#'
analyze.ddf <- function(my.file, 
                        win.start=0.83, 
                        win.length=5, 
                        vers=3, 
                        max.frame,
                        y.lim=512, 
                        a.1=18, 
                        a.2=0, 
                        cut=50, 
                        m.d.cs=0.2, 
                        pfad.mult, 
                        n.angle=3, 
                        signal.neg=FALSE, 
                        do.plot=FALSE, 
                        do.plot.vec=FALSE, 
                        do.plot.hots=FALSE, 
                        do.plot.hots.vecs=FALSE, 
                        plot.hots.mult=FALSE,
                        n.cores=1, 
                        df.t=NULL, 
                        frames.pack=600, 
                        maxdist=5, 
                        fastclust=FALSE, 
                        speeditup=FALSE, 
                        floodclust=FALSE, 
                        floodCpp=TRUE,
                        only.watch=TRUE, 
                        schwarm.find=FALSE, 
                        regel.schwarm=NULL, 
                        which.regel.schwarm=NULL, 
                        save.preprocess=FALSE, 
                        zeitstempel, 
                        min.length=10, 
                        fisch="Object", 
                        max.row=50, 
                        center=FALSE, 
                        long=TRUE, 
                        file.regel, 
                        klass.regel, 
                        sdcorr=FALSE, 
                        FUN=mean, 
                        FUN.ow=mean, 
                        save.jpg=TRUE,
                        save.movie=TRUE, 
                        delete.jpg=TRUE,
                        each=5,
                        wait.jpg=FALSE,
                        win=FALSE){
  
  cat(as.character(Sys.time())," Beginn Preprocessing \n")
  
  # Preprocessing durchfuehren
  preprocessed <- preprocess.data(my.file=my.file,
                                  win.start=win.start,
                                  win.length=win.length,
                                  vers=vers,
                                  y.lim=y.lim,
                                  a.1=a.1, 
                                  a.2=a.2, 
                                  cut=cut,
                                  m.d.cs=m.d.cs,
                                  pfad.mult=pfad.mult,
                                  n.angle=n.angle,
                                  signal.neg=signal.neg,
                                  do.plot=do.plot,
                                  do.plot.vec=do.plot.vec,
                                  n.cores=n.cores,
                                  df.t=df.t,
                                  frames.pack=frames.pack,
                                  maxdist=maxdist,
                                  fastclust=fastclust,
                                  speeditup=speeditup,
                                  floodclust=floodclust,
                                  floodCpp=floodCpp,
                                  only.watch=only.watch,
                                  schwarm.find=schwarm.find,
                                  regel.schwarm=regel.schwarm,
                                  which.regel.schwarm=which.regel.schwarm,
                                  save.preprocess=save.preprocess,
                                  plot.hots.mult=plot.hots.mult,
                                  save.jpg=save.jpg,
                                  zeitstempel=zeitstempel,
                                  save.movie=save.movie, 
                                  delete.jpg=delete.jpg,
                                  each=each,
                                  wait.jpg=wait.jpg,
                                  win=win
  )
  
  cat(as.character(Sys.time())," Preprocessing fertig \n")
  cat("Preprocessing fertig \n")
  
  ########################################################################
  # ACHTUNG: Das ist noch unschoen, fuer das Echtzeit-Warnsystem brauche ich
  #	hier kein .RData auf die Festplatte legen, das kostet nur Rechenzeit
  # LB 25.11.13: Erledigt
  ########################################################################
  
  # Ergebnisse laden
  # LB 25.11.13: Brauche ich nicht mehr, Ergebnisse des Preprocessing
  #		werden jetzt an das Objekt preprocessed uebergeben
  # load(file=paste(pfad.mult,"Results_n_angle_",n.angle,"_a2_",a.2,".RData",sep=""))
  
  # Falls nur Schwaerme gefunden wurden
  if(any(preprocessed$is.schwarm.vec)){
    
    # Auch Output generieren, falls zwar keine Objekte, aber Schwaerme 
    #		erkannt wurden
    
    cat("Ausschliesslich Schwaerme erkannt \n")
    
    only.schwarm <- TRUE
    
    generate.output(pfad.mult=pfad.mult,
                    zeitstempel=zeitstempel,
                    vars.hot.list=preprocessed$vars.hot.list,
                    is.schwarm.vec=preprocessed$is.schwarm.vec,
                    only.schwarm=only.schwarm)
    
    
  }else{
    
    # Nur weiter mit Klassifikation, wenn mindestens ein Hotspot 
    #		erkannt wurde
    if(!is.null(preprocessed$track.evals[1][[1]])){
      
      # Evaluierung gesamt & Plot eines Histogramms der Pfadlaengen
      tracks.prob <- track.eval.total(track.evals=preprocessed$track.evals, 
                                      min.length=min.length,
                                      save.plot=FALSE,
                                      pfad.mult=paste(pfad.mult,zeitstempel,"/", sep=""),
                                      # 											  fisch=fisch, 
                                      # MW: Wird das benoetigt? Filename?
                                      # LB: Nein, hab es auskommentiert
                                      a.2=a.2)
      
      print(tracks.prob)
      
      # 	if(max(sum.prepro$n.hots)>= min.length){
      if(tracks.prob$tracks > 0 & tracks.prob$unproblematisch > 0){
        
        do.classification <- TRUE
        
        if(do.plot.hots){
          
          # Plot der gefunden Hotspots
          sum.prepro <- summary.func(track.ddf.list=preprocessed$track.ddf.list,
                                     track.eval.out.list=preprocessed$track.eval.out.list,
                                     merkmale.out.list=preprocessed$merkmale.out.list,
                                     t.grid.list=preprocessed$t.grid.list,
                                     n.angle=n.angle)
          # MW: Fehlermeldung "In cor(t(merks)) : Standardabweichung ist Null", Unwichtig?
          # 	Sollte aber abgefangen werden.
          
          cat("Sum.prepro fertig \n")
          
          plot.prepro(hot=TRUE,
                      fisch=fisch,
                      pfad=paste(pfad.mult, zeitstempel,"/", sep=""),
                      sum.prepro=sum.prepro,
                      min.length=min.length,
                      cart.coord=preprocessed$cart.coord,
                      hots.y01.list=preprocessed$hots.y01.list,
                      dims.list=preprocessed$dims.list,
                      track.ddf.list=preprocessed$track.ddf.list,
                      max.row=max.row,
                      n.cores=n.cores,
                      a.2=a.2
          )
          
          cat("Plot Hotspots fertig \n")
          
          if(do.plot.hots.vecs){
            
            plot.prepro(hot=FALSE,
                        fisch=fisch,
                        pfad=paste(pfad.mult, zeitstempel,"/", sep=""),
                        sum.prepro=sum.prepro,
                        min.length=min.length,
                        merkmale.out.list=preprocessed$merkmale.out.list,
                        max.row=max.row,
                        n.cores=n.cores,
                        a.2=a.2
            )
            
            cat("Plot Vektoren fertig")
          }
        }
        
        
      }else{
        cat("Kein Objekt mit mehr als", min.length, "Hotspots erkannt \n")
        
        do.classification <- FALSE
      }
      
      
      
    }else{ 
      cat("Kein Objekt erkannt \n")
      
      do.classification <- FALSE
      
    }
    
    # Klassifikation hier abbrechen, wenn keine Tracks gefunden wurden
    # => Dann nur eine Zeile mit 0 ausgeben
    # => Oder Schwarmoutput
    
    cat(do.classification, "\n")
    
    if(do.classification){
      
      
      ###########################################################################
      # 2.) Daten fuer Klassifikation vorbereiten
      # 	und
      # 3.) Zaehlung
      ###########################################################################
      
      class.out.list <- prepare.and.count(track.ddf.list=preprocessed$track.ddf.list,
                                          track.eval.out.list=preprocessed$track.eval.out.list,
                                          merkmale.out.list=preprocessed$merkmale.out.list,
                                          t.grid.list=preprocessed$t.grid.list,
                                          n.angle=n.angle,
                                          min.length=min.length,
                                          center=center,
                                          long=long,
                                          only.watch=only.watch,
                                          file.regel=file.regel,
                                          klass.regel=klass.regel,
                                          sdcorr=sdcorr,
                                          FUN=FUN,
                                          FUN.ow=FUN.ow)
      
      ###########################################################################
      # 4.) Output-Datei erstellen
      ###########################################################################
      
      generate.output(do.classification=class.out.list$do.classification,
                      which.select=class.out.list$which.select,
                      pred.class=class.out.list$pred.class,
                      post.obj=class.out.list$post.obj,
                      sum.prepro.roh=class.out.list$sum.prepro.roh,
                      frames.pack=frames.pack,
                      track.eval.out.list=preprocessed$track.eval.out.list,
                      merkmale.out.list=preprocessed$merkmale.out.list,
                      max.frame=max.frame,
                      pfad.mult=pfad.mult,
                      zeitstempel=zeitstempel)
      
    }else{
      
      generate.output(do.classification=do.classification,
                      pfad.mult=pfad.mult,
                      zeitstempel=zeitstempel)
      
    }
    
    
  }
  
  cat(as.character(Sys.time())," Fischerkennung abgeschlossen \n")
  cat("Fischerkennung abgeschlossen. \n")
  #End
}


####################################################################
# Die Funktion analyze_ddf() wrapped analyze.ddf damit dort alles englisch ist 
#   und nur die Parameter aus der Diss/dem Paper übergeben werden müssen/können
####################################################################

#' Analyze complete .ddf file
#' 
#' This function analyzes a given .ddf file, i.e., carries out the entire 
#'  preprocessing with identification of hotspots, tracking, extraction of 
#'  features and finally classifies the identified objects.
#'
#' @param ddf.file Filename of .ddf file to be analyzed
#' @param win.start Start of the sonar window, possibly extracted before via 
#'  \code{\link{get.version}}
#' @param win.length Length of the sonar window, possibly extracted before via 
#'  \code{\link{get.version}}
#' @param vers Version of the .ddf file, possibly extracted before via 
#'  \code{\link{get.version}}
#' @param max.frame Total number of frames of the given video, 
#'  possibly extracted before via \code{\link{get.version}} 
#' @param y.lim Number of pixels to analyze in each beam, counted from the 
#'  camera
#' @param a Threshold for centered data
#' @param cut Minimal size of clusters in number of pixels. Smaller clusters are 
#'  deleted before the tracking and hence not considered as hotspots.
#' @param m.d.cs Maximal distance for which two hotspots can be assigned the 
#'  same tracking number
#' @param folder.output Folder for resulting plots and output
#' @param n.angle Number of watch hands per quarter, i.e., total number of 
#'  watch hands is \code{4 x n.angle}
#' @param do.plots \code{TRUE}: Plots of preprocessing are saved, default is \code{FALSE}
#' @param n.cores Number of cores if parallelization is needed, default is \code{1}
#' @param frames.pack Number of frames to be analyzed at once, possibly the 
#'  total number of frames \code{max.frame} of the given video. If \code{frames.pack} is 
#'    less than \code{max.frame}, the video is splitted in parts of 
#'    \code{frames.pack} frames and each part is analyzed separately. 
#' @param df.t Number of degrees of freedom for the splines. Default \code{NULL} 
#'  results in \code{c(100, 25, round(n.frames/2))} where \code{n.frames} is the 
#'  total number of frames of the analyzed video, i.e. \code{length(frames.pack)}.
#' @param save.preprocess \code{TRUE}: Results of preprocessing are saved, 
#'  default is \code{FALSE}
#' @param timestamp Time stamp of .ddf file
#' @param min.length Minimal length of tracks. Tracks with less hotspots are deleted.
#' @param file.classrule File name of classification rule
#' @param type.classrule Type of classification rule, i.e., lda, qda...
#' @param save.movie \code{TRUE}: Video of raw signal and cleaned signal is created, 
#'  default is \code{TRUE}
#' @param each Interval of images for the video, for example default \code{5} 
#'  means, that each \code{5}th image is saved in the video
#' @param win Specify \code{TRUE} if you are running under windows, then 
#'  parallelization is not possible
#' @return As output, a table is generated which lists the found objects with 
#'  predicted class and computed features at \code{folder.output}.
#' @export
#'
analyze_ddf <- function(ddf.file, 
                        win.start=0.83, 
                        win.length=5, 
                        vers=3, 
                        max.frame,
                        y.lim=512, 
                        a=18, 
                        cut=50, 
                        m.d.cs=0.2, 
                        folder.output, 
                        n.angle=3, 
                        do.plots=FALSE, 
                        n.cores=1, 
                        frames.pack=600, 
                        df.t=NULL, 
                        save.preprocess=FALSE, 
                        timestamp, 
                        min.length=10, 
                        file.classrule, 
                        type.classrule, 
                        save.movie=TRUE, 
                        each=5,
                        win=FALSE){
  
  analyze.ddf(my.file=ddf.file, 
              win.start=win.start, 
              win.length=win.length, 
              vers=vers, 
              y.lim=y.lim, 
              max.frame=max.frame,
              a.1=a, 
              a.2=0, 
              cut=cut, 
              m.d.cs=m.d.cs, 
              pfad.mult=folder.output, 
              n.angle=n.angle, 
              signal.neg=FALSE, 
              do.plot=do.plots, 
              do.plot.vec=do.plots, 
              do.plot.hots=do.plots, 
              do.plot.hots.vecs=do.plots, 
              plot.hots.mult=do.plots,
              n.cores=n.cores, 
              df.t=df.t, 
              frames.pack=frames.pack, 
              maxdist=5, 
              fastclust=FALSE, 
              speeditup=FALSE, 
              floodclust=FALSE, 
              floodCpp=TRUE,
              only.watch=TRUE, 
              schwarm.find=FALSE, 
              regel.schwarm=NULL, 
              which.regel.schwarm=NULL, 
              save.preprocess=save.preprocess, 
              zeitstempel=timestamp, 
              min.length=min.length, 
              fisch="Object", 
              max.row=50, 
              center=FALSE, 
              long=TRUE, 
              file.regel=file.classrule, 
              klass.regel=type.classrule, 
              sdcorr=FALSE, 
              FUN=mean, 
              FUN.ow=mean, 
              save.jpg=TRUE,
              save.movie=save.movie, 
              delete.jpg=TRUE,
              each=each,
              wait.jpg=FALSE,
              win=win)
  
}