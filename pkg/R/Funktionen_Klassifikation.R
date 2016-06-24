###########################################################################
# Fischprojekt
# Funktionen zur Klassifikation
# Datum der ersten Version: 10.09.12
# Letzte Aenderung: 19.02.13
# Autor: Ludwig Bothmann
###########################################################################

require(MASS)
# library(pls)
# library(genridge)
# require(descr)
# require(nnet)
# library(VGAM)
require(tree)
require(e1071)

###########################################################################
# Die Funktion cv.func() soll eine K-fache Kreuzvalidierung durchfuehren.
#	Die Funktion kann benutzt werden, wenn direkt klassifiziert werden
#	soll, also fuer Objekt- und Hotspotebene. Wenn kumuliert klassifiziert
#	werden soll, also aus mehreren Hotspots eine Vorhersage fuer
#	ein Objekt erstellt werden soll, muss die Funktion cv.func.cumul()
#	benutzt werden. 
#
# Uebergeben werden muss:
#	- K:			Anzahl Schichten (default = 5)
#	- seed:		Seed
#	- daten:		Der Datensatz fuer die Klassifikation
#	- regel:		Zu benutzende Klassifikationsregel
#	- groups:	Anzahl Gruppen (default = 2)
#
# Ausgegeben wird:
#	- error:		Anzahl falsch klassifizierter Objekte in jeder Schicht
#	- result:	Vorhergesagt und wahre Klasse fuer jedes Objekt
#	- post:		Posteriori-Wahrscheinlichkeiten pro Objekt
#
###########################################################################

cv.func <- function(K=5,
						  seed=as.numeric(Sys.time()),	
						  daten,
						  regel,
						  groups=2
						  ){
	
	set.seed(seed)
	
	# Anzahl der Objekte
	n <- nrow(daten)
	
	# Cutpoints fuer die Schichten des Datensatzes
	cuts <- quantile(1:n,probs=seq(0,1,length=K+1))
	cut.points <- c(0,round(cuts[-1]))
	
	# Zufaellige Ordnung der Objekt-IDs
	id.random <- sample(1:n, n)	
	
	# Matrizen, in die die Anzahl der fehlerhaft klassifizierten Objekte
	#	geschrieben werden.
	error.total <- matrix(0, nrow=K, ncol=2, dimnames=list(1:K,c("Falsch","Gesamt")))
	
	# Matrizen mit prognostizierten und wahren Klassen
	result.total<- NULL
	
	# Matrizen mit Posteriori-Wahrscheinlichkeiten
	post.total<- NULL
	
	for(i in 1:K){
		
		# Datensatz teilen
		id.test <- id.random[(cut.points[i]+1):cut.points[i+1]]
		data.tr <- daten[-id.test,]
		data.te <- daten[id.test,]
		
		# Neu: Klassifikation mit Hilfsfunktion, an die nur Trainings- und
		#		Testdatensatz sowie Anzahl Gruppen und Regel uebergeben wird.
		class.out <- class.func(data.tr = data.tr,
														data.te = data.te,
														regel = regel,
														groups = groups)
		
		# Alt: Per Hand
# 		# Regel bestimmen ...
# 		if(regel == "lda"){
# 			class.obj <- lda(klasse~., data=data.tr)
# 		}
# 		if(regel == "qda"){
# 			class.obj <- qda(klasse~., data=data.tr)
# 		}
# 		if(regel == "logit"){
# 			glm.obj <- glm(klasse~., data=data.tr, family="binomial")
# 		}
# 		if(regel == "mnlogit"){
# 			mnl.obj <- multinom(klasse~., data=data.tr)
# 		}
# 		if(regel == "tree"){
# 			tre.obj <- tree(klasse~., data=data.tr)
# 		}
# 		if(regel == "svm"){
# 			svm.obj <- svm(klasse~., data=data.tr, probability=TRUE)
# 		}
# 		
# 		# ... und auswerten
# 		if(is.element(regel, c("lda", "qda"))){
# 			pred.obj <- predict(class.obj, data.te)	
# 			pred.class <- pred.obj$class
# 			if(groups==2){
# 				pred.post <- pred.obj$post[,1]
# 				post <- matrix(pred.post, ncol=1, 
# 									dimnames=list(names(pred.post),"Post-W' fuer 1"))
# 			}
# 			if(groups==3){
# 				post <- pred.obj$post
# 			}
# 		}
# 		if(regel=="logit"){
# 			pred.obj <- predict(glm.obj, data.te, type="response")	
# 			pred.class <- (pred.obj>=0.5)*1 + 1
# 			post <- matrix(1-pred.obj, ncol=1, 
# 								dimnames=list(names(pred.class),"Post-W' fuer 1"))
# 		}
# 		if(regel=="mnlogit"){
# 			pred.obj <- predict(mnl.obj, data.te, type="probs")	
# 			pred.class <- predict(mnl.obj, data.te, type="class")
# 			post <- pred.obj
# 		}
# 		if(regel=="tree"){
# 			pred.obj <- predict(tre.obj, data.te)	
# 			pred.class <- predict(tre.obj, data.te, type="class")
# 			post <- pred.obj
# 		}
# 		if(regel=="svm"){
# 			pred.obj <- predict(svm.obj, data.te, probability=TRUE)
# 			pred.class <- predict(svm.obj, data.te)
# 			if(groups==2){
# 				pred.post <- attr(pred.obj,"probabilities")[,1]
# 				post <- matrix(pred.post, ncol=1, 
# 									dimnames=list(names(pred.post),"Post-W' fuer 1"))
# 			}
# 			if(groups==3){
# 				post <- attr(pred.obj,"probabilities")
# 			}
# 		}
# 		
# 		
# 		result <- cbind(as.numeric(pred.class), data.te[,1])
# 		error <- c(sum(result[,1]!=result[,2]), nrow(data.te))
# 		error.total[i,] <- error
# 		result.total <- rbind(result.total,result)
# 		post.total <- rbind(post.total, post)
		
		error.total[i,] <- class.out$error[1:2]
		result.total <- rbind(result.total,class.out$result)
		post.total <- rbind(post.total, class.out$post)
		
	}
	
	colnames(result.total) <- c("Prognose","Wahr")
	
	return(list(error = error.total,
					result = result.total,
					post = post.total))
}

###########################################################################
# Die Funktion cv.func.all() fuehrt die K-fache Kreuzvalidierung fuer
#		alle 5 Klassifizierungsmethoden durch und schreibt die Fehlerraten
#		raus
#
# Uebergeben werden muss:
#	- K:			Anzahl Schichten (default = 5)
#	- seed:		Seed
#	- daten:		Der Datensatz fuer die Klassifikation
#	- regel:		Zu benutzende Klassifikationsregeln als Vektor
#	- groups:	Anzahl Gruppen (default = 2)
# - cumul:	TRUE => Hotspotebene, sonst Objektebene. def = FALSE
# - ID.match:			Bei Hotspotebene: Welcher Hotspot gehoert zu welchem Objekt?
# - match.klasse:	Bei Hotspotebene: Wahre Klasse
#
# Ausgegeben wird:
#	- error:		Matrix mit Fehlklassifikationsraten
#
###########################################################################

cv.func.all <- function(K=5,
												seed=as.numeric(Sys.time()),	
												daten,
												regel=c("lda", "qda", "mnlogit", "tree", "svm"),
												groups=2,
												cumul=FALSE,
												ID.match=NULL,
												match.klasse=NULL
){

	zs <- NULL
	
	for(i in 1:length(regel)){
	
		# Kreuzvalidierung mit der i-ten Regel durchfuehren, entweder 
		#		normal (Objeketeben) oder
		#		kumuliert (Hotspotebenw)
		if(cumul==FALSE){	
			cv.i <- cv.func(K=K, 
											seed=seed, 
											daten=daten, 
											regel=regel[i],
											groups=groups)
		}else{	
			cv.i <- cv.func.cumul(K=K, 
														seed=seed, 
														daten=daten, 
														regel=regel[i],
														ID.match=ID.match,
														match.klasse=match.klasse,
														groups=groups)
		}
		
		# Fehler speichern
		zs <- rbind(zs,apply(X=cv.i$error, MARGIN=2, sum))
		
	}
	
	# Zeilennamen: Benutzte Regel
	rownames(zs) <- regel
	
	# Fehlklassifikaitonsrate ausrechnen
	zs <- cbind(zs, round(zs[,1]/zs[,2],3))
	
	# Spaltennamen
	colnames(zs) <- c("Missclassified", "Total","Ratio")
	
	return(zs)
}

###########################################################################
# Die Funktion cv.func.cumul() soll eine K-fache Kreuzvalidierung 
#	durchfuehren.
#	Die Funktion kann benutzt werden, wenn aus mehreren Hotspots eine 
#	Vorhersage fuer ein Objekt erstellt werden soll.
#
# Uebergeben werden muss:
#	- K:					Anzahl Schichten (default = 5)
#	- seed:				Seed
#	- daten:				Der Datensatz fuer die Klassifikation
#	- regel:				Zu benutzende Klassifikationsregel
#	- ID.match:			Matrix in der die ID-Namen zu den Zeilennummern 
#								gematcht werden.
#	- match.klasse:	ID-Namen wird Klasse zugeordnet
#	- groups:			Anzahl der Gruppen (default = 2)
#
# Ausgegeben wird:
#	- error:		Anzahl falsch klassifizierter Objekte in jeder Schicht
#	- result:	Vorhergesagt und wahre Klasse fuer jedes Objekt
#	- post:		Posteriori-Wahrscheinlichkeiten pro Objekt
#
###########################################################################

cv.func.cumul <- function(K=5,
								  seed=as.numeric(Sys.time()),
								  daten,
								  regel,
								  ID.match,
								  match.klasse,
								  groups=2
){
	
	set.seed(seed)
	
	# Anzahl der Objekte
	n <- nrow(match.klasse)
	
	# Cutpoints fuer die Schichten des Datensatzes
	cuts <- quantile(1:n,probs=seq(0,1,length=K+1))
	cut.points <- c(0,round(cuts[-1]))
	
	# Zufaellige Ordnung der Objekt-IDs
	id.random <- sample(1:n, n)	
	
	# Matrizen, in die die Anzahl der fehlerhaft klassifizierten Objekte
	#	geschrieben werden.
	error.total <- matrix(0, nrow=K, ncol=2, dimnames=list(1:K,c("Falsch","Gesamt")))
	
	# Matrizen mit prognostizierten und wahren Klassen
	result.total<- NULL
	
	# Matrizen mit Posteriori-Wahrscheinlichkeiten
	post.total<- NULL
	
	for(i in 1:K){
		
		####################
		# Datensatz teilen
		####################
		
		# IDs der Testobjekte (ungeordnet)
		id.test.u <- match.klasse[id.random[(cut.points[i]+1):cut.points[i+1]],1]
		
		# Zeilennummern der Hotspots zu den Testobjekten
		hots.test <- which(is.element(ID.match[,1],id.test.u))
		
		# Hier werden die IDs zu den Zeilennamen gematcht
		ID.match.test <- ID.match[hots.test,]
		
		# IDs der Testobjekte (geordnet)
		#####################################################################
		# ACHTUNG: Nicht nach Alphabet sondern wie in data.hot, ist das ok?
		# 	=> Ja, das passt, am Ende stimmt alles
		# ;-) Vermutlich, weil ich es gar nicht mehr benutze
		#	=> Kann raus
		#####################################################################
# 		id.test <- unique(ID.match.test[,1])
		
		# Trainings- und Testdatensatz
		data.tr <- daten[-hots.test,]
		data.te <- daten[hots.test,]
		
		# Regel bestimmen ...
		if(regel == "lda"){
 			class.obj <- lda(klasse~., data=data.tr)
		}
		if(regel == "qda"){
			class.obj <- qda(klasse~., data=data.tr)
		}
		if(regel == "logit"){
			glm.obj <- glm(klasse~., data=data.tr, family="binomial")
		}
		if(regel == "mnlogit"){
			mnl.obj <- multinom(klasse~., data=data.tr)
		}
		if(regel == "tree"){
			tre.obj <- tree(klasse~., data=data.tr)
		}
		if(regel == "svm"){
			svm.obj <- svm(klasse~., data=data.tr, probability=TRUE)
		}
		
		#####################
		# ... und auswerten
		#####################
		
		# 3 Gruppen
		if(groups==3){
			
			if(is.element(regel, c("lda", "qda"))){
				pred.post <- predict(class.obj, newdata=data.te)$post
			}
			
			if(regel=="mnlogit"){
				pred.post <- predict(mnl.obj, data.te, type="probs")
			}
			
			if(regel == "tree"){
				pred.post <- predict(tre.obj, data.te)
			}
			if(regel=="svm"){
				pred.obj <- predict(svm.obj, data.te, probability=TRUE)
				pred.post <- attr(pred.obj,"probabilities")
			}
			
			pred.post.1 <- tapply(pred.post[,1], INDEX=ID.match.test[,1], mean)
			pred.post.2 <- tapply(pred.post[,2], INDEX=ID.match.test[,1], mean)
			pred.post.3 <- tapply(pred.post[,3], INDEX=ID.match.test[,1], mean)
			pred.response <- cbind(pred.post.1, pred.post.2, pred.post.3)
			pred.class <- apply(pred.response,1,which.max)
			post <- pred.response
			
			if(regel=="logit") cat("Logit bei 2 Gruppen nicht sinnvoll! \n")
		}
		
		# 2 Gruppen
		if(groups==2){
			
			if(is.element(regel, c("lda", "qda"))){					
				pred.post <- predict(class.obj, data.te)$post[,1]
				pred.response <- tapply(pred.post, INDEX=ID.match.test[,1], mean)
				pred.class <- (pred.response<0.5)*1 + 1
				post <- matrix(pred.response, ncol=1, 
									dimnames=list(names(pred.class),"Post-W' fuer 1"))
			}	
			
			if(regel=="logit"){
				pred.post <- predict(glm.obj, data.te, type="response")	
				pred.response <- tapply(pred.post, INDEX=ID.match.test[,1], mean)
				pred.class <- (pred.response>=0.5)*1 + 1
				# Hier >=, da 1 = nein, bei lda und qda 1 = ja
				# und ich brauche 1 und 2, statt 1 und 0
				post <- matrix(1-pred.response, ncol=1, 
									dimnames=list(names(pred.class),"Post-W' fuer 1"))
			}

			if(regel=="svm"){ 
				pred.obj <- predict(svm.obj, data.te, probability=TRUE)
				pred.post <- attr(pred.obj,"probabilities")[,1]
				pred.response <- tapply(pred.post, INDEX=ID.match.test[,1], mean)
				pred.class <- (pred.response<0.5)*1 + 1
				post <- matrix(pred.response, ncol=1, 
									dimnames=list(names(pred.class),"Post-W' fuer 1"))
			}
			
			if(regel=="tree") cat("Baeume fuer 2 Gruppen nicht implementiert. \n")
			if(regel=="mnlogit") cat("MN-Logit bei 2 Gruppen sinnvoll? \n")
			
		}
			
		# Vergleich von prognostizierter und wahrer Klasse
		result <- cbind(pred.class, 
							 as.numeric(match.klasse[which(is.element(match.klasse[,1],
																	  names(pred.class))),2]))		
		
		# Wie viele falsch klassifiziert?
		error <- c(sum(result[,1]!=result[,2]),length(id.test.u))
		
		# Alles abspeichern
		error.total[i,] <- error
		result.total <- rbind(result.total,result)
		post.total <- rbind(post.total, post)
	}
	
	colnames(result.total) <- c("Prognose","Wahr")
	
	return(list(error = error.total,
					result = result.total,
					post = post.total))
}


###########################################################################
# Die Funktion cv.func.kombi() soll eine K-fache Kreuzvalidierung 
#	durchfuehren. Hier wird die Information von Objekt- und Hotspotebene
#	miteinander verknuepft, was hoffentlich eine bessere Prognose erlaubt.
#
# Uebergeben werden muss:
#	- K:					Anzahl Schichten
#	- seed:				Seed
#	- daten.obj:		Der Datensatz auf Objektebene
#	- daten.hot:		Der Datensatz auf Hotspotebene
#	- regel.obj:		Zu benutzende Klassifikationsregel O-Ebene
#	- regel.hot:		Zu benutzende Klassifikationsregel H-Ebene
#	- ID.match:			Matrix in der die ID-Namen zu den Zeilennummern 
#								gematcht werden.
#	- match.klasse:	ID-Namen wird Klasse zugeordnet
#	- groups:			Anzahl der Gruppen (default = 2)
#	- weights:			Gewichte der Objektebene
#
# Ausgegeben wird:
#	- error:		Anzahl falsch klassifizierter Objekte in jeder Schicht
#	- result:	Vorhergesagt und wahre Klasse fuer jedes Objekt
#	- post:		Posteriori-Wahrscheinlichkeiten pro Objekt
#
###########################################################################
 
cv.func.kombi <- function(K,
								  seed=as.numeric(Sys.time()),
								  daten.obj,
								  daten.hot,
								  regel.obj,
								  regel.hot,
								  ID.match,
								  match.klasse,
								  groups=2,
								  weights=0.5
){

	# CV fuer Objektebene duchfuehren
	cv.obj <- cv.func(K=K, 
							seed=seed, 
							daten=daten.obj, 
							regel=regel.obj,
							groups=groups)
	
	# CV fuer Hotspotebene durchfuehren
	cv.hot <- cv.func.cumul(K=K, 
									seed=seed, 
									daten=daten.hot, 
									regel=regel.hot, 
									ID.match=ID.match, 
									match.klasse=match.klasse, 
									groups=groups)
	
	if(groups==2){
		
		# Posteriori-W' gleich ordnen
		post.obj <- cv.obj$post[order(row.names(cv.obj$post))]
		post.hot <- cv.hot$post[order(row.names(cv.hot$post))]
		
		names(post.obj) <- row.names(cv.obj$post)[order(row.names(cv.obj$post))]
		names(post.hot) <- row.names(cv.hot$post)[order(row.names(cv.hot$post))]
		
		# Posteriori-W' miteinander verrechnen
		komb <- data.frame(post.obj, post.hot)
		komb$post <- post.obj * weights + post.hot * (1-weights)
		komb$prognose <- (komb$post<0.5)*1 + 1
		komb$wahr <- as.numeric(match.klasse[order(match.klasse[,1]),2])
	}
	
	if(groups==3){
		
		# Posteriori-W' gleich ordnen
		post.obj <- cv.obj$post[order(row.names(cv.obj$post)),]
		post.hot <- cv.hot$post[order(row.names(cv.hot$post)),]
		
		# Posteriori-W' miteinander verrechnen
		post.komb <- post.obj * weights + post.hot * (1-weights)
		
		# Prognose
		komb.class <- apply(post.komb,1,which.max)
		
		# In einem data.frame
		komb <- data.frame(post.komb, prognose=komb.class)
		komb$wahr <- as.numeric(match.klasse[order(match.klasse[,1]),2])				
	}
	
	# Gesamtfehler
	error <- c(sum(komb$prognose!=komb$wahr),nrow(daten.obj))
	names(error) <- c("Falsch","Gesamt")
	
	return(list(error=error,
					komb=komb))
}


###########################################################################
# Die Funktion class.func() ist eine Hilfsfunktion fuer die Klassifizierung
#
# Uebergeben werden muss:
#	- K:			Anzahl Schichten (default = 5)
#	- seed:		Seed
#	- daten:		Der Datensatz fuer die Klassifikation
#	- regel:		Zu benutzende Klassifikationsregel
#	- groups:	Anzahl Gruppen (default = 2)
#
# Ausgegeben wird:
#	- error:		Anzahl falsch klassifizierter Objekte in jeder Schicht
#	- result:	Vorhergesagt und wahre Klasse fuer jedes Objekt
#	- post:		Posteriori-Wahrscheinlichkeiten pro Objekt
#
###########################################################################


#' @import MASS
#' @import e1071
#' @import tree
#' @import nnet
class.func <- function(data.tr,
											 data.te,
											 regel,
											 groups=2
){
	
# 	prior <- rep(1, groups)/groups
	
	# Regel bestimmen ...
	if(regel == "lda"){
		class.obj <- lda(klasse~., data=data.tr
# 										 ,prior=prior
										 )
	}
	if(regel == "qda"){
		class.obj <- qda(klasse~., data=data.tr
# 										 ,prior=prior
		)
	}
	if(regel == "logit"){
		glm.obj <- glm(klasse~., data=data.tr, family="binomial")
	}
	if(regel == "mnlogit"){
		mnl.obj <- multinom(klasse~., data=data.tr)
	}
	if(regel == "tree"){
		tre.obj <- tree(klasse~., data=data.tr)
	}
	if(regel == "svm"){
		svm.obj <- svm(klasse~., data=data.tr, probability=TRUE)
	}
	
	# ... und auswerten
	if(is.element(regel, c("lda", "qda"))){
		pred.obj <- predict(class.obj, data.te)	
		pred.class <- pred.obj$class
		if(groups==2){
			pred.post <- pred.obj$post[,1]
			post <- matrix(pred.post, ncol=1, 
										 dimnames=list(names(pred.post),"Post-W' fuer 1"))
		}
		if(groups>=3){
			post <- pred.obj$post
		}
	}
	if(regel=="logit"){
		pred.obj <- predict(glm.obj, data.te, type="response")	
		pred.class <- (pred.obj>=0.5)*1 + 1
		post <- matrix(1-pred.obj, ncol=1, 
									 dimnames=list(names(pred.class),"Post-W' fuer 1"))
	}
	if(regel=="mnlogit"){
		pred.obj <- predict(mnl.obj, data.te, type="probs")	
		pred.class <- predict(mnl.obj, data.te, type="class")
		post <- pred.obj
	}
	if(regel=="tree"){
		pred.obj <- predict(tre.obj, data.te)	
		pred.class <- predict(tre.obj, data.te, type="class")
		post <- pred.obj
	}
	if(regel=="svm"){
		pred.obj <- predict(svm.obj, data.te, probability=TRUE)
		pred.class <- predict(svm.obj, data.te)
		if(groups==2){
			pred.post <- attr(pred.obj,"probabilities")[,1]
			post <- matrix(pred.post, ncol=1, 
										 dimnames=list(names(pred.post),"Post-W' fuer 1"))
		}
		if(groups>=3){
			post <- attr(pred.obj,"probabilities")
		}
	}
	
	
	result <- cbind(as.numeric(pred.class), data.te[,1])
	error.1 <- c(sum(result[,1]!=result[,2]), nrow(data.te))
	error <- c(error.1, error.1[1]/error.1[2])
	
	out <- list(result=result,
							error=error,
							post=post)
	
	return(out)
}

###########################################################################
# Die Funktion subsample.func() soll den Prognosefehler subsampeln
#
# Uebergeben werden muss:
# - B:				Anzahl Subsample-Iterationen
# - sub.test:	Anteil Testdaten an Gesamtdaten
#	- seed:			Seed
#	- daten:		Der Datensatz fuer die Klassifikation
#	- regel:		Zu benutzende Klassifikationsregel
#	- groups:		Anzahl Gruppen (default = 2)
#
# Ausgegeben wird:
#	- error:		Anzahl falsch klassifizierter Objekte in jeder Schicht
#	- result:		Vorhergesagt und wahre Klasse fuer jedes Objekt
#	- post:			Posteriori-Wahrscheinlichkeiten pro Objekt
#
###########################################################################

subsample.func <- function(B=100,
													 sub.test=1/4,
													 seed=as.numeric(Sys.time()),	
													 daten,
													 regel,
													 groups=2
){
	
	
	#Seed setzen
	set.seed(seed)
	
	# Anzahl Objekte
	n <- nrow(daten)

	# Matrizen, in die die Anzahl der fehlerhaft klassifizierten Objekte
	#	geschrieben werden.
	error.total <- matrix(0, 
												nrow=B, 
												ncol=3, 
												dimnames=list(1:B,c("Falsch","Gesamt","Fehler")))
	
	# Matrix, in der steht, welches Objekt in welcher Iteration im 
	#		Test oder Training war
	#	1: Training
	# 0: Test
	set <- matrix(0, ncol=B, nrow=n)
	
	# Array, in das die Kreuztabellen geschrieben werden
	cross <- array(dim=c(3,3,B))
	
	# B mal Trainingsdatensatz und Testdatensatz erstellen und 
	# klassifizieren
	for(b in 1:B){
	
		# Welche Objekte kommen in den Testdatensatz?
		test <- sample(x=1:n, 
									 round(n*sub.test), 
									 replace=FALSE)
		
		# Trainings- und Testdatensatz erstellen
		data.tr <- daten[-test,]
		data.te <- daten[test,]
		
		# Klassifikation durchfuehren
		class.out <- class.func(data.tr = data.tr,
														data.te = data.te,
														groups = groups,
														regel = regel)
		
		# Merken, welches Objekt in welchem Datensatz war
		set[-test,b] <- 1
		
		# Fehlklassifizierung
		error.total[b,] <- class.out$error
		
		# Kreuztabelle wahre / vorhergesagte Klasse
		# => relative Haeufigkeit
		cross[,,b] <- table(class.out$result[,1], 
												class.out$result[,2]) / nrow(class.out$result)
		
	}
	
	# Gemittelte Kreuztabelle
	cross.tab <- round(apply(cross, 1:2, mean),4)
	
	error <- apply(error.total, 2, mean)
	prog.fehler <- error[1]/error[2]
	names(prog.fehler) <- "Prognosefehler"
	
	out <- list(error.total = error.total,
							set = set,
							prog.fehler=prog.fehler,
							cross.tab = cross.tab)
	
	return(out)

}

###########################################################################
# Die Funktion boot.func() soll den Prognosefehler bootstrappen
#
# Uebergeben werden muss:
# - B:				Anzahl Bootstrap-Iterationen
#	- seed:			Seed
#	- daten:		Der Datensatz fuer die Klassifikation
#	- regel:		Zu benutzende Klassifikationsregel
#	- groups:		Anzahl Gruppen (default = 2)
#
# Ausgegeben wird:
#	- error:		Anzahl falsch klassifizierter Objekte in jeder Schicht
#	- result:		Vorhergesagt und wahre Klasse fuer jedes Objekt
#	- post:			Posteriori-Wahrscheinlichkeiten pro Objekt
#
###########################################################################

boot.func <- function(B=100,
											seed=as.numeric(Sys.time()),	
											daten,
											regel,
											groups=2
){
	
	
	#Seed setzen
	set.seed(seed)
	
	# Anzahl Objekte
	n <- nrow(daten)
	
	# Matrizen, in die die Anzahl der fehlerhaft klassifizierten Objekte
	#	geschrieben werden.
	error.total <- matrix(0, 
												nrow=B, 
												ncol=3, 
												dimnames=list(1:B,c("Falsch","Gesamt","Fehler")))
	
	# Matrix, in der steht, welches Objekt in welcher Iteration im 
	#		Test oder Training war
	#	1-x: Anzahl im Trainingsdatensatz
	# 0: Test
	set <- matrix(0, ncol=B, nrow=n)
	
	# Array, in den die Kreuztabellen geschrieben werden
	cross <- array(dim=c(3,3,B))
	
	# Gewicht fuer das gewichtete arithemtische Mittel
	# => Nur bei Bootstrap noetig
	weight <- vector(length=B)
	
	# B mal Trainingsdatensatz und Testdatensatz erstellen und 
	# klassifizieren
	
	# ToDo: boot() benutzen!
	for(b in 1:B){
		
		# Welche Objekte kommen in den Trainingsdatensatz?
		# => Es koennen mehrere Objekte auch oefter vorkommen
		train <- sample(x=1:n, 
										n,
										replace=TRUE)
		
		# Welche Objekte kommen in den Testdatensatz?
		test <- which(!is.element(1:n, train))
		
		# Trainings- und Testdatensatz erstellen
		data.tr <- daten[train,]
		data.te <- daten[test,]
		
		# Klassifikation durchfuehren
		class.out <- class.func(data.tr = data.tr,
														data.te = data.te,
														groups = groups,
														regel = regel)
		
		# Merken, welches Objekt in welchem Datensatz war
		set[sort(unique(train)),b] <- table(train)
		
		# Fehlklassifizierung
		error.total[b,] <- class.out$error
		
		
		# Kreuztabelle wahre / vorhergesagte Klasse
		# => relative Haeufigkeit
		cross[,,b] <- table(class.out$result[,1], 
												class.out$result[,2]) / nrow(class.out$result)
		
		weight[b] <- nrow(class.out$result)
		
	}
	
	# Gemittelte Kreuztabelle
	cross.tab <- round(apply(cross, 1:2, weighted.mean, w=weight),4)
	
	error <- apply(error.total, 2, mean)
	prog.fehler <- error[1]/error[2]
	names(prog.fehler) <- "Prognosefehler"
	
	out <- list(error.total = error.total,
							set = set,
							prog.fehler=prog.fehler
							,cross.tab = cross.tab
							)
	
	return(out)
	
}

###########################################################################
# Die Funktion cv.B.func() soll den Prognosefehler B mal K mal
#		kreuzvalidieren
#
# Uebergeben werden muss:
# - B:				Anzahl Iterationen aussen
# - K: 				Anzahl Iterationen innen
#	- seed:			Seed
#	- daten:		Der Datensatz fuer die Klassifikation
#	- regel:		Zu benutzende Klassifikationsregel
#	- groups:		Anzahl Gruppen (default = 2)
#
# Ausgegeben wird:
#	- error:		Anzahl falsch klassifizierter Objekte in jeder Schicht
#	- result:		Vorhergesagt und wahre Klasse fuer jedes Objekt
#	- post:			Posteriori-Wahrscheinlichkeiten pro Objekt
#
###########################################################################

cv.B.func <- function(B=100,
											K,
											seed=as.numeric(Sys.time()),	
											daten,
											regel,
											groups=2
){
	
	
	#Seed setzen
	set.seed(seed)
	
	# Anzahl Objekte
	n <- nrow(daten)
	
	# Matrizen, in die die Anzahl der fehlerhaft klassifizierten Objekte
	#	geschrieben werden.
	error.total <- matrix(0, 
												nrow=B, 
												ncol=3, 
												dimnames=list(1:B,c("Falsch","Gesamt","Fehler")))
	
	# Matrix, in der steht, welches Objekt in welcher Iteration in 
	#		welchem Datensatz was
	# 1 - K: Schicht, in der das Objekt war
	set <- matrix(0, ncol=B, nrow=n)
	
	# Array, in den die Kreuztabellen geschrieben werden
	cross <- array(dim=c(3,3,B))
	
	# B mal Trainingsdatensatz und Testdatensatz erstellen und 
	# klassifizieren
	for(b in 1:B){
		
		# Objekte in Schichten einteilen
		# ACHTUNG: Passt das so? Groesse der Schichten ist so etwas
		#		variabel
		schichten <- sample(x=1:K, 
												n,
												replace=TRUE)
		
		error.K <- NULL
		cross.k <- matrix(0, nrow=3, ncol=3)
		
		for(k in 1:K){
		
			# Testdaten sind die Objekte der Schicht k
			test.k <- which(schichten==k)
			
			# Trainings- und Testdatensatz erstellen
			data.tr <- daten[-test.k,]
			data.te <- daten[test.k,]
			
			# Klassifikation durchfuehren
			class.out <- class.func(data.tr = data.tr,
															data.te = data.te,
															groups = groups,
															regel = regel)
			
			# Fehler in Iterationen k = 1, ..., K
			error.K <- rbind(error.K, class.out$error)
		
			# Kreuztabelle wahre / vorhergesagte Klasse
			# => relative Haeufigkeiten
			cross.k <- cross.k + table(class.out$result[,1], 
																 class.out$result[,2])
		}
		
	
		# Fehlklassifizierung
		error.total[b,1:2] <- apply(error.K,2,sum)[1:2]
		error.total[b,3] <- error.total[b,1] / error.total[b,2]
		
		# Kreuztablle relative Haeufigkeiten Iteration b
		cross[,,b] <- cross.k / sum(cross.k)
		
	}
	
	# Kreuztabelle relative Haeufigkeiten gesamt
	cross.tab <- round(apply(cross, 1:2, mean),4)
	
	error <- apply(error.total, 2, mean)
	prog.fehler <- error[1]/error[2]
	names(prog.fehler) <- "Prognosefehler"
	
	out <- list(error.total = error.total,
							set = set,
							prog.fehler=prog.fehler,
							cross.tab = cross.tab)
	
	return(out)
	
}


########################################################################
# Hinweis: Die folgenden 3 Funktionen mit "svd" im Namen sind veraltet.
#		Die neuere Version ist, dass ich die Designmatrix fuer svd bestimme
#		und dann die normale Klassifizierungsfunktion benutze.
########################################################################


###########################################################################
# Die Funktion cv.svd() soll eine K-fache Kreuzvalidierung durchfuehren.
#		In dieser Funktion ist die Klassifikation mittels SVD implementiert. 
#
# Uebergeben werden muss:
#	- K:						Anzahl Schichten (default = 5)
#	- seed:					Seed
#	- daten:				Der Datensatz fuer die Klassifikation
#	- regel:				Zu benutzende Klassifikationsregel
#	- ID.match:			Matrix in der die ID-Namen zu den Zeilennummern 
#										gematcht werden.
#	- match.klasse:	ID-Namen wird Klasse zugeordnet
#	- groups:				Anzahl der Gruppen (default = 3)
#
# Ausgegeben wird:
#	- error:				Anzahl falsch klassifizierter Objekte in jeder Schicht
#	- result:				Vorhergesagt und wahre Klasse fuer jedes Objekt
#	- post:					Posteriori-Wahrscheinlichkeiten pro Objekt
#
###########################################################################

cv.svd <- function(K=5,
									 seed=as.numeric(Sys.time()),
									 daten,
									 regel,
									 ID.match,
									 match.klasse,
									 groups=3,
									 n.angle,
									 nv=3,
									 nu=3,
									 delete=NULL
){
	
	set.seed(seed)
	
	# Zeilen loeschen?
	if(!is.null(delete)){
		ID.match <- ID.match[-delete,]
		daten <- daten[-delete,]
	}
	
	# Anzahl der Objekte
	n <- nrow(match.klasse)
	
	# Cutpoints fuer die Schichten des Datensatzes
	cuts <- quantile(1:n,probs=seq(0,1,length=K+1))
	cut.points <- c(0,round(cuts[-1]))
	
	# Zufaellige Ordnung der Objekt-IDs
	id.random <- sample(1:n, n)	
	
	# Matrizen, in die die Anzahl der fehlerhaft klassifizierten Objekte
	#	geschrieben werden.
	error.total <- matrix(0, nrow=K, ncol=2, dimnames=list(1:K,c("Falsch","Gesamt")))
	
	# Matrizen mit prognostizierten und wahren Klassen
	result.total<- NULL
	
	# Matrizen mit Posteriori-Wahrscheinlichkeiten
	post.total<- NULL
	
	# Spalten in daten, die die Uhrzeiger enthalten
	index.var <- c(3:(2+n.angle*4))
	
# 	# Winkel, zu denen die Uhrzeiger gehoeren
# 	winkel <-  seq(0, (360-90/n.angle), length=4*n.angle)
	
	for(k in 1:K){
		
# 		cat("k=",k)
		####################
		# Datensatz teilen
		####################
		
		# IDs der Testobjekte (ungeordnet)
		id.test.u <- match.klasse[id.random[(cut.points[k]+1):cut.points[k+1]],1]
		
		# Zeilennummern der Hotspots zu den Testobjekten
		hots.test <- which(is.element(ID.match[,1],id.test.u))
		
		# Hier werden die IDs zu den Zeilennamen gematcht
		ID.match.test <- ID.match[hots.test,]
		
		# IDs der Testobjekte (geordnet)
		#####################################################################
		# ACHTUNG: Nicht nach Alphabet sondern wie in data.hot, ist das ok?
		# 	=> Ja, das passt, am Ende stimmt alles
		# 		;-) Vermutlich, weil ich es gar nicht mehr benutze
		#		=> Kann raus
		#####################################################################
		# 		id.test <- unique(ID.match.test[,1])
		
		# Trainings- und Testdatensatz
		data.tr <- daten[-hots.test,]
		data.te <- daten[hots.test,]
		
		# Datensaetze fuer die Klassifikation erstellen mit Hilfsfunktion
		#	svd.data.func()
		class.data <- svd.data.func(data.tr, 
																data.te, 
																index.var,
																nv=nv,
																nu=nv)
		
		diff.tr <- class.data$diff.tr
		diff.te <- class.data$diff.te
		
		#####################################################################
		# Ist das so sinnvoll? Jetzt haette ich auf jeden Fall einen 
		#	Trainings- und Testdatensatz fuer die Klassifikation
		#####################################################################
		
		# Regel bestimmen ...
		if(regel == "lda"){
			class.obj <- lda(klasse~., data=diff.tr)
		}
		if(regel == "qda"){
			class.obj <- qda(klasse~., data=diff.tr)
		}
		if(regel == "mnlogit"){
			mnl.obj <- multinom(klasse~., data=diff.tr)
		}
		if(regel == "tree"){
			tre.obj <- tree(klasse~., data=diff.tr)
		}
		if(regel == "svm"){
			svm.obj <- svm(klasse~., data=diff.tr, probability=TRUE)
		}
		
		#####################
		# ... und auswerten
		#####################
		
		# 3 Gruppen
		if(groups==3){
			
			if(is.element(regel, c("lda", "qda"))){
				pred.post <- predict(class.obj, newdata=diff.te)$post
			}
			
			if(regel=="mnlogit"){
				pred.post <- predict(mnl.obj, diff.te, type="probs")
			}
			
			if(regel == "tree"){
				pred.post <- predict(tre.obj, diff.te)
			}
			if(regel=="svm"){
				pred.obj <- predict(svm.obj, diff.te, probability=TRUE)
				pred.post <- attr(pred.obj,"probabilities")
			}
			
# 			pred.post.1 <- tapply(pred.post[,1], INDEX=ID.match.test[,1], mean)
# 			pred.post.2 <- tapply(pred.post[,2], INDEX=ID.match.test[,1], mean)
# 			pred.post.3 <- tapply(pred.post[,3], INDEX=ID.match.test[,1], mean)
# 			pred.response <- cbind(pred.post.1, pred.post.2, pred.post.3)
			pred.class <- apply(pred.post,1,which.max)
			post <- pred.post
			
		}
		
		# Vergleich von prognostizierter und wahrer Klasse
		result <- cbind(pred.class, 
							 as.numeric(match.klasse[which(is.element(match.klasse[,1],
							 													  names(pred.class))),2]))		
		
		# Wie viele falsch klassifiziert?
		error <- c(sum(result[,1]!=result[,2]),length(id.test.u))
		
		# Alles abspeichern
		error.total[k,] <- error
		result.total <- rbind(result.total,result)
		post.total <- rbind(post.total, post)
	}
	
	colnames(result.total) <- c("Prognose","Wahr")
	
	return(list(error = error.total,
					result = result.total,
					post = post.total))
}


###########################################################################
# Die Funktion svd.data.func() soll den Training- und Testdatensatz
#		fuer die Klassifikation mit der SVD-Methode berechnen und laeuft 
#		innerhalb von cv.svd()
#
# Uebergeben werden muss:
#	- data.tr:			Rohdaten Training
#	- data.te:			Rohdaten Test
# - index.var:		Vektor, der angibt in welchen Spalten der Rohdaten die
#										relevanten Kovariablen fuer die Designmatrix stehen
#	- nv:						Anzahl Spalteneffekte
#	- nu:						Anzahl Zeileneffekte (ist glaube ich unwichtig)
#
# Ausgegeben wird:
#	- diff.tr:			Datensatz fuer Klassifikation Training
#	- diff.te:			Datensatz fuer Klassifikation Test
#
###########################################################################

svd.data.func <- function(data.tr, 
													data.te,
													index.var,
													nv=3,
													nu=3){
	
	# Welche Zeilen sind Aale bzw. Forellen?
	index.tr.aal <-  data.tr$art=="Aal"
	index.tr.for <-  data.tr$art=="Forelle"
	
# 	index.te.aal <-  data.te$art=="Aal"
# 	index.te.for <-  data.te$art=="Forelle"
	
	# Singulaerwertzerlegung der Trainings-Designmatrix
	svd.aal <-  svd(data.tr[index.tr.aal,index.var],nu=nu,nv=nv)
	svd.for <- svd(data.tr[index.tr.for,index.var],nu=nu,nv=nv)
	
	# Rechte Singulaervektoren
	v.aal <-   svd.aal$v
	v.for <-  svd.for$v
	
	# 		# Regressionskoeffizienten
	# 		a.aal <-  as.matrix(data.tr[,index.var]) %*% v.aal
	# 		a.for <-  as.matrix(data.tr[,index.var]) %*% v.for
	
	# Approximation der Trainigs-Designmatrix als Aal- und Forellenmodell		
	X.aal <- as.matrix(data.tr[,index.var]) %*% v.aal %*% t(v.aal)
	X.for <- as.matrix(data.tr[,index.var]) %*% v.for %*% t(v.for)
	
	diff.tr <- svd.data.help(X.aal = X.aal,
													 X.for = X.for,
													 daten = data.tr,
													 index.var = index.var)
	
	##########################
	# Fuer Test-Datensatz:
	##########################
	
	# Approximation der Test-Designmatrix als Aal- und Forellenmodell		
	X.aal.te <- as.matrix(data.te[,index.var]) %*% v.aal %*% t(v.aal)
	X.for.te <- as.matrix(data.te[,index.var]) %*% v.for %*% t(v.for)
	
	diff.te <- svd.data.help(X.aal = X.aal.te,
													 X.for = X.for.te,
													 daten = data.te,
													 index.var = index.var)
	
	
	# Zurueckgeben:
	out <- list(diff.tr=diff.tr,
							diff.te=diff.te)
	return(out)
}


###########################################################################
# Die Funktion svd.data.help() laeuft innerhalb von svd.data.func().
#		Ausgelagert, damit ich nicht Code doppelt habe.
#
# Uebergeben werden muss:
#	- X.aal:				Designmatrix als Aal-Modell
#	- X.for:				Designmatrix als Forellen-Modell
#	- daten:				Daten Original
# - index.var:		Vektor, der angibt in welchen Spalten von Daten die
#										relevanten Kovariablen fuer die Designmatrix stehen
#
# Ausgegeben wird:
#	- diff.out:			Datensatz fuer Klassifikation 
#
###########################################################################
svd.data.help <- function(X.aal,
													X.for,
													daten,
													index.var){
	
	diffX <- c()
	diff.einfach <- c()
	steuer <- c()
	
	# Alle Objekte durchgehen
	for(i in 1:length(unique(daten$ID))){
		
		# 			print(i)
		
		# Welche Zeilen gehoeren zu Objekt i?
		index.i <-  daten$ID == unique(daten$ID)[i]
		
		# Nur diese Zeilen der Designmatrix
		data.i <- as.matrix(daten[index.i,index.var])
		
		colnames(data.i) <- NULL
		
		# Diese Zeilen der Prognosematrizen X.aal und X.for
		X.aal.i <- X.aal[index.i,]
		X.for.i <- X.for[index.i,]
		
		# Berechne fuer jeden Uhrzeiger die mittleren quadrierten Differenzen
		#		der Werte fuer Objekt i in der Designmatrix und in den
		#		Prognosematrizen.
		# => Vektor mit 2*p Eintraegen
		
		diffX.i <- c(apply(FUN=mean, MARGIN=2, ((data.i - X.aal.i)**2)),
								 apply(FUN=mean, MARGIN=2, ((data.i - X.for.i)**2)))
		
		# Alles in eine Matrix
		diffX <- rbind(diffX, diffX.i)
		
		# ID und Art des Objekts
		steuer <-  rbind(steuer, c(unique(daten$ID)[i], daten$art[index.i][1]))
		
		# Arithmetisches Mittel der p Uhrzeiger in der Designmatrix
		diff.einfach <- rbind(diff.einfach,
														 apply(FUN=mean, MARGIN = 2, data.i))
	}
	
	rownames(diffX) <- unique(daten$ID)
	diffX <- as.data.frame(diffX)
	
	
	# Datensatz fuer Klassifikation
	diff.class <- diffX[,1:length(index.var)]-diffX[,length(index.var) + 1:length(index.var)]
	diff.out <- data.frame(diff.class, diff.einfach[,-1], klasse=steuer[,2])
	
	return(diff.out)
}