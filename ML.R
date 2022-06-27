  #
  ## MACHINE LEARNING AVEC R
  #
  #install.packages("caret")
  library("lattice")
  library(caret)
  #install.packages("glmnet")
  library("Matrix")
  library(glmnet)
  library(stringr)
  #install.packages("randomForest", repos="http://R-Forge.R-project.org")
  library(randomForest)
  library(tidyr)
  #install.packages("kernlab")
  library(kernlab)
  library(tidyr)
  library(e1071)
  #if (!requireNamespace("remotes", quietly = TRUE)) {
  # install.packages("remotes")}
  #remotes::install_github("robingenuer/VSURF")
  library(VSURF)
  library(dplyr)
  library(ggplot2)
  #install.packages("cowplot")
  library(cowplot)
  
  #fixation de l'aléatoire 
  set.seed(154)
  
  # récupération du noms des fichiers échantillons-condition pour les 5 couples de sous-matrices 
  # de test et d'entrainement
  #files_train_condition <- list.files(path = "~/Documents/k-mers_level/samples-conditions", pattern = "sampleshuf.train", all.files = FALSE,
  #                                    full.names = FALSE, recursive = FALSE,
  #                                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  
  #files_test_condition <- list.files(path = "~/Documents/k-mers_level/samples-conditions", pattern = "sampleshuf.test", all.files = FALSE,
  #                                   full.names = FALSE, recursive = FALSE,
  #                                   ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  # indique le chemin vers le répertoire contenant les sous-matrices
  pathtrain = "~/Urine/train/mat/6mer"
  pathtest = "~/Urine/train/mat/6mer"
  
  patterntrain = "train"
  patterntest = "test"
  
  # récupération de la liste des fichiers présents dans le répertoire correspondant au type choisi
  # pour les sous-matrices d'entrainement
  files_train <- list.files(path = pathtrain, pattern = patterntrain, all.files = FALSE,
                            full.names = FALSE, recursive = FALSE,
                            ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  # pour les sous-matrices de test
  files_test <- list.files(path = pathtest, pattern = patterntest, all.files = FALSE,
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  # Fonction qui permet d'attribuer 5 noms de variables aux 5 fichiers qui souhaitent être lu sous R
  svg_var_sub <- function(files_names){
    nb_files <- 5
    data_names <- vector("list",length=nb_files)
    
    for (i in 1 : nb_files) {
      data_names[i] <- strsplit(files_names[i], split=".tsv")
    }
    
    return(list(data_names = data_names, files_names = files_names))
  }
  
  # appel de la fonction pour les 5 sous-matrices d'entrainement
  setwd(pathtrain)
  res <- svg_var_sub(files_train)
  data_names <- res$data_names
  files_names <- res$files_names
  # enregistrement du contenu des fichiers des les variables R
  for (i in 1:5) {
    print(files_names[i])
    assign(data_names[[i]], 
           read.table(files_names[i], header = TRUE, row.names = 1))
  }
  
  # appel de la fonction pour les 5 sous-matrices de test
  setwd(pathtest)
  res <- svg_var_sub(files_test)
  data_names <- res$data_names
  files_names <- res$files_names
  # enregistrement du contenu des fichiers dans les variables R
  for (i in 1:5) {
    print(files_names[i])
    assign(data_names[[i]], 
           read.table(files_names[i], header = TRUE, row.names = 1))
  }
  
  # idem pour les fichiers échantillons-condition pour les 5 sous-matrices de test et d'entrainement
  #setwd("~/Documents/k-mers_level/samples-conditions")
  #test
  res <- svg_var_sub(files_test)#_condition)
  data_names <- res$data_names
  files_names <- res$files_names
  for (i in 1:5) {
    print(files_names[i])
    assign(data_names[[i]], 
           read.table(files_names[i], header = TRUE))
  }
  #train
  res <- svg_var_sub(files_train)#_condition)
  data_names <- res$data_names
  files_names <- res$files_names
  for (i in 1:5) {
    print(files_names[i])
    assign(data_names[[i]], 
           read.table(files_names[i], header = TRUE))
  }
  
  # fonction pour nommer la ou les colonnes
  rename_colnames <- function(df){
    colnames(df) <- c("id", "condition")
    return(df)
  }
  
  #li_samples_conditions_files <- c("sampleshuf.test0", "sampleshuf.test1", "sampleshuf.test2", "sampleshuf.test3", "sampleshuf.test4", "sampleshuf.train0", "sampleshuf.train1", "sampleshuf.train2", "sampleshuf.train3", "sampleshuf.train4")
  #li_contenu_samples_conditions <- list(sampleshuf.test0, sampleshuf.test1, sampleshuf.test2, sampleshuf.test3, sampleshuf.test4, sampleshuf.train0, sampleshuf.train1, sampleshuf.train2, sampleshuf.train3, sampleshuf.train4)
  
  # appel de la fonction pour renommer sur les fichiers éhantillons-conditions
  #for (i in 1:10){
  #  assign(li_samples_conditions_files[i], 
  #         rename_colnames(li_contenu_samples_conditions[[i]]))
  #}
  
  
  #li_train = list(`matrix_merge_counts.train0`,`matrix_merge_counts.train1`,`matrix_merge_counts.train2`,`matrix_merge_counts.train3`,`matrix_merge_counts.train4`)
  #li_trainl = c('matrix_merge_counts.train0','matrix_merge_counts.train1','matrix_merge_counts.train2','matrix_merge_counts.train3','matrix_merge_counts.train4')
  #li_test = list(`contigEvaluate.test0`,`contigEvaluate.test1`,`contigEvaluate.test2`,`contigEvaluate.test3`,`contigEvaluate.test4`)
  
  
  #for (i in 1:5){
   # assign(li_trainl[i], suppr_col(li_train[[i]]))
  #}
  
  # maj du contenu de la liste
  #li_train = list(`matrix_merge_counts.train0`,`matrix_merge_counts.train1`,`matrix_merge_counts.train2`,`matrix_merge_counts.train3`,`matrix_merge_counts.train4`)
  
  
  # Etape de transposition des matrices qui sont sous la forme : k-mers en ligne et echantillons en colonne. 
  # Cette etape est nécessaire puisque l'apprentissage automatique fonctionne avec les variables en colonne donc les k-mers en colonne         
  #transpose_matrixResSig_matrixML <- function(datas){
   # data <- data.frame(t(datas))
    #coldata <- data.frame(features = rownames(data), condition = as.factor(str_split(rownames(data), "_", simplify = TRUE)[,1]), row.names = 1 )
    #return(cbind(data, coldata))
  #}
  
  
  #
  ## Renomage mais pas necessaire pour les train et test
  #
  #names=c()
  #for (i in 1:length(train0)){
   # names=c(names,colnames(train0[1]))
    #colnames(train0[i])=colnames(train0[1])
  #}
  #colnames(train0)=names
  #rownames(test3_mer)=test3_mer[,1]
  #test3_mer=test3_mer[,-1]
  
  
  #rownames(train3_mer)=train3_mer[,1]
  #train3_mer=train3_mer[,-1]
  # Etape de filtration sur le nombre de counts minimum (nbminc) 
  # et le nombre de patients minimum chez lequels ont retrouves les k-mers (nbminp)
  # ces 2 variables peuvent être changées
  
  #nbminc = 10
  #nbminp = 5
  
  #ligne <- which(rowSums(`matrix_merge_counts.train0` >= nbminc) >= nbminp)
  #`matrix_merge_counts.train0` = `matrix_merge_counts.train0`[ligne,]
  #contigEvaluate.test0 = contigEvaluate.test0[ligne,]
  
  #ligne <- which(rowSums(`matrix_merge_counts.train1` >= nbminc) >= nbminp)
  #`matrix_merge_counts.train1` = `matrix_merge_counts.train1`[ligne,]
  #contigEvaluate.test1 = contigEvaluate.test1[ligne,]
  
  #ligne <- which(rowSums(`matrix_merge_counts.train2` >= nbminc) >= nbminp)
  #`matrix_merge_counts.train2` = `matrix_merge_counts.train2`[ligne,]
  #contigEvaluate.test2 = contigEvaluate.test2[ligne,]
  
  #ligne <- which(rowSums(`matrix_merge_counts.train3` >= nbminc) >= nbminp)
  #`matrix_merge_counts.train3` = `matrix_merge_counts.train3`[ligne,]
  #contigEvaluate.test3 = contigEvaluate.test3[ligne,]
  
  #ligne <- which(rowSums(`matrix_merge_counts.train4` >= nbminc) >= nbminp)
  #`matrix_merge_counts.train4` = `matrix_merge_counts.train4`[ligne,]
  #contigEvaluate.test4 = contigEvaluate.test4[ligne,]
  
  # maj du contenu de la liste
  #li_train = list(`matrix_merge_counts.train0`,`matrix_merge_counts.train1`,`matrix_merge_counts.train2`,`matrix_merge_counts.train3`,`matrix_merge_counts.train4`)
  #li_test = list(`contigEvaluate.test0`,`contigEvaluate.test1`,`contigEvaluate.test2`,`contigEvaluate.test3`,`contigEvaluate.test4`)
  
  
  # Appel de la fonction de transposition de la matrice
  # FOLD
  #0
  test0 <- t(test0)#transpose_matrixResSig_matrixML(li_test[[1]])
  train0 <- t(train0)#transpose_matrixResSig_matrixML(li_train[[1]])
  
  #1
  test1=t(test1)#test1 <- transpose_matrixResSig_matrixML(li_test[[2]])
  train1=t(train1)#train1 <- transpose_matrixResSig_matrixML(li_train[[2]])
  
  #2
  test2=t(test2)#test2 <- transpose_matrixResSig_matrixML(li_test[[3]])
  train2=t(train2)#train2 <- transpose_matrixResSig_matrixML(li_train[[3]])
  
  #3
  test3=t(test3)#test3 <- transpose_matrixResSig_matrixML(li_test[[4]])
  train3=t(train3)#train3 <- transpose_matrixResSig_matrixML(li_train[[4]])
  
  #4
  test4=t(test4)#test4 <- transpose_matrixResSig_matrixML(li_test[[5]])
  train4=t(train4)#train4 <- transpose_matrixResSig_matrixML(li_train[[5]])
  
  
  #
  ## Aprentissage avec Ridge Regression, Random Forest et Lasso Regression
  #
  
  fit_pred_confusionMatrix_for_1_kfold_glmnet <- function(seed, train, test, method){
    set.seed(seed)
    train.x <- as.matrix(train[, -ncol(train)])
    train.y <- sub("[[:digit:]]*$", "", rownames(train))
    train.y=as.factor(train.y)
    test.x <- as.matrix(test[, -ncol(test)])
    test.y <- sub("[[:digit:]]*$", "", rownames(test)) #suppression des chiffres en fin de caractère/définis conditions
    test.y=as.factor(test.y)
    
    if (method=="ridge"){
      set.seed(seed)
      mdl <- cv.glmnet(x = train.x, y = train.y, family = "binomial", alpha = 0, type.measure = "deviance")
      pred <- predict(mdl, test.x ,type="class",s=c(0))
      pred <- factor(c(pred), levels = c("normal_", "tumeur_"))
      mat <- confusionMatrix(data=pred, reference=test.y, positive = 'tumeur_')
    } else if (method=="rf"){
      set.seed(seed)
      mdl <- randomForest(x = train.x, y = train.y)
      pred <- predict(mdl, test.x)
      mat <- confusionMatrix(data=as.factor(pred), reference=test.y, positive = 'tumeur_')
    }else if (method=="lasso"){
      set.seed(seed)
      mdl <- cv.glmnet(x = train.x, y = train.y, family = "binomial", alpha = 1, type.measure = "deviance")
      pred <- predict(mdl, test.x ,type="class",s=c(0))
      pred <- factor(c(pred), levels = c("normal_", "tumeur_"))
      mat <- confusionMatrix(data=pred, reference=test.y, positive = 'tumeur_')
    }else if (method=="vsurf"){
      set.seed(seed)
      thers = VSURF_thres(train.x, train.y, ntree = 2000,mtry = max(floor(ncol(train.x)/3), 1), nfor.thres = 50, nmin = 1,
                          RFimplem = "randomForest", parallel = TRUE, clusterType = "PSOCK",
                          ncores = parallel::detectCores() - 1, verbose = TRUE)
      
      interp = VSURF_interp(train.x, train.y, ntree = 2000, thers$varselect.thres,
                            nfor.interp = 25, nsd = 1, RFimplem = "randomForest",
                            parallel = TRUE, ncores = parallel::detectCores() - 1,
                            clusterType = "PSOCK", verbose = TRUE)
      
      pred = VSURF_pred(train.x, train.y, ntree = 2000, interp$err.interp,
                        interp$varselect.interp, nfor.pred = 25, nmj = 1,
                        RFimplem = "randomForest", parallel = TRUE, ncores = parallell::detectCores()
                        - 1, verbose = TRUE)
      
      train.xrf = train.x[,pred$varselect.pred]
      test.xrf = test.x[,pred$varselect.pred]
      mdl <- randomForest(x = train.xrf, y = train.y)
      pred <- predict(mdl, test.xrf)
      mat <- confusionMatrix(data=as.factor(pred), reference=test.y, positive = 'tumeur_')
    }
    
    return(mat)
  }         
  
  
  # fold0
  mat_0 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train0 , test0, "ridge")
  
  graphe <- data.frame(balancedaccuracy = c(mat_0$byClass["Balanced Accuracy"]), fold = "fold1" , method="Ridge")
  graphe_Sensitivity <- data.frame(Sensitivity  = c(mat_0$byClass["Sensitivity"]), fold = "fold1", method="Ridge")
  graphe_Specificity <- data.frame(Specificity  = c(mat_0$byClass["Specificity"]), fold = "fold1", method="Ridge")
  
  # fold1
  mat_1 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train1 , test1, "ridge")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_1$byClass["Balanced Accuracy"]), fold = "fold2" , method="Ridge"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_1$byClass["Sensitivity"]), fold = "fold2", method="Ridge"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_1$byClass["Specificity"]), fold = "fold2", method="Ridge"))
  
  # fold2
  mat_2 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train2, test2, "ridge")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_2$byClass["Balanced Accuracy"]), fold = "fold3", method="Ridge"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_2$byClass["Sensitivity"]), fold = "fold3", method="Ridge"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_2$byClass["Specificity"]), fold = "fold3", method="Ridge"))
  
  # fold3
  mat_3 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train3, test3,  "ridge")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_3$byClass["Balanced Accuracy"]), fold = "fold4", method="Ridge"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_3$byClass["Sensitivity"]), fold = "fold4", method="Ridge"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_3$byClass["Specificity"]), fold = "fold4", method="Ridge"))
  
  # fold4 
  mat_4 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train4, test4, "ridge") 
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_4$byClass["Balanced Accuracy"]), fold = "fold5", method="Ridge"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_4$byClass["Sensitivity"]), fold = "fold5", method="Ridge"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_4$byClass["Specificity"]), fold = "fold5", method="Ridge"))
  
  #############################
  
  # fold0
  mat_0 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train0 , test0, "rf")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_0$byClass["Balanced Accuracy"]), fold = "fold1" , method="Random Forest"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_0$byClass["Sensitivity"]), fold = "fold1", method="Random Forest"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_0$byClass["Specificity"]), fold = "fold1", method="Random Forest"))
  
  # fold1
  mat_1 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train1 , test1, "rf")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_1$byClass["Balanced Accuracy"]), fold = "fold2" , method="Random Forest"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_1$byClass["Sensitivity"]), fold = "fold2", method="Random Forest"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_1$byClass["Specificity"]), fold = "fold2", method="Random Forest"))
  
  # fold2
  mat_2 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train2, test2, "rf")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_2$byClass["Balanced Accuracy"]), fold = "fold3", method="Random Forest"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_2$byClass["Sensitivity"]), fold = "fold3", method="Random Forest"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_2$byClass["Specificity"]), fold = "fold3", method="Random Forest"))
  
  # fold3
  mat_3 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train3, test3,  "rf")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_3$byClass["Balanced Accuracy"]), fold = "fold4", method="Random Forest"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_3$byClass["Sensitivity"]), fold = "fold4", method="Random Forest"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_3$byClass["Specificity"]), fold = "fold4", method="Random Forest"))
  
  # fold4
  mat_4 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train4, test4, "rf") 
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_4$byClass["Balanced Accuracy"]), fold = "fold5", method="Random Forest"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_4$byClass["Sensitivity"]), fold = "fold5", method="Random Forest"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_4$byClass["Specificity"]), fold = "fold5", method="Random Forest"))
  
  ###############################
  
  # fold0
  mat_0 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train0 , test0, "lasso")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_0$byClass["Balanced Accuracy"]), fold = "fold1" , method="lasso"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_0$byClass["Sensitivity"]), fold = "fold1", method="lasso"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_0$byClass["Specificity"]), fold = "fold1", method="lasso"))
  
  # fold1
  mat_1 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train1 , test1, "lasso")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_1$byClass["Balanced Accuracy"]), fold = "fold2" , method="lasso"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_1$byClass["Sensitivity"]), fold = "fold2", method="lasso"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_1$byClass["Specificity"]), fold = "fold2", method="lasso"))
  
  # fold2
  mat_2 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train2, test2, "lasso")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_2$byClass["Balanced Accuracy"]), fold = "fold3", method="lasso"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_2$byClass["Sensitivity"]), fold = "fold3", method="lasso"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_2$byClass["Specificity"]), fold = "fold3", method="lasso"))
  
  # fold3
  mat_3 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train3, test3,  "lasso")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_3$byClass["Balanced Accuracy"]), fold = "fold4", method="lasso"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_3$byClass["Sensitivity"]), fold = "fold4", method="lasso"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_3$byClass["Specificity"]), fold = "fold4", method="lasso"))
  
  # fold4 
  mat_4 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train4, test4, "lasso") 
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_4$byClass["Balanced Accuracy"]), fold = "fold5", method="lasso"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_4$byClass["Sensitivity"]), fold = "fold5", method="lasso"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_4$byClass["Specificity"]), fold = "fold5", method="lasso"))
  
  # fold0
  mat_0 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train0 , test0, "vsurf")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_0$byClass["Balanced Accuracy"]), fold = "fold1" , method="vsurf"))
  graphe_Sensitivity  <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_0$byClass["Sensitivity"]), fold = "fold1", method="vsurf"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_0$byClass["Specificity"]), fold = "fold1", method="vsurf"))
  
  # fold1
  mat_1 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train1 , test1, "vsurf")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_1$byClass["Balanced Accuracy"]), fold = "fold2" , method="vsurf"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_1$byClass["Sensitivity"]), fold = "fold2", method="vsurf"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_1$byClass["Specificity"]), fold = "fold2", method="vsurf"))
  
  # fold2
  mat_2 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train2, test2, "vsurf")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_2$byClass["Balanced Accuracy"]), fold = "fold3", method="vsurf"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_2$byClass["Sensitivity"]), fold = "fold3", method="vsurf"))
  graphe_Specificity  <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_2$byClass["Specificity"]), fold = "fold3", method="vsurf"))
  
  # fold3
  mat_3 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train3, test3,  "vsurf")
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_3$byClass["Balanced Accuracy"]), fold = "fold4", method="vsurf"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_3$byClass["Sensitivity"]), fold = "fold4", method="vsurf"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_3$byClass["Specificity"]), fold = "fold4", method="vsurf"))
  
  # fold4
  mat_4 <- fit_pred_confusionMatrix_for_1_kfold_glmnet(154, train4, test4, "vsurf") 
  
  graphe <- rbind(graphe, data.frame(balancedaccuracy = c(mat_4$byClass["Balanced Accuracy"]), fold = "fold5", method="vsurf"))
  graphe_Sensitivity <- rbind(graphe_Sensitivity, data.frame(Sensitivity  = c(mat_4$byClass["Sensitivity"]), fold = "fold5", method="vsurf"))
  graphe_Specificity <- rbind(graphe_Specificity, data.frame(Specificity  = c(mat_4$byClass["Specificity"]), fold = "fold5", method="vsurf"))
  
  
  
  
  #
  ## GRAPHE
  #
  
  ggplot(data=graphe, aes(x=fold, y=balancedaccuracy, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=round(balancedaccuracy,2)), vjust=-0.25, color="black", position =  position_dodge(0.9) ,size=3.5)+
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() + ggtitle("K-fold = 5, Balanced accuracy") +  ylim(0,1)
  
  
  
  ggplot(data=graphe_Sensitivity, aes(x=fold, y=Sensitivity, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=round(Sensitivity,2)), vjust=-0.25, color="black", position =  position_dodge(0.9) ,size=3.5)+
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() + ggtitle("K-fold = 5, Sensitivity")+ ylim(0,1)
  
  
  
  ggplot(data=graphe_Specificity, aes(x=fold, y=Specificity, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=round(Specificity,2)), vjust=-0.25, color="black", position =  position_dodge(0.9) ,size=3.5)+
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() + ggtitle("K-fold = 5, Specificity") + ylim(0,1)
  
  
  
  df.summary_accuracy <- graphe %>%
    group_by(method) %>%
    summarise(
      sd = sd(balancedaccuracy, na.rm = TRUE),
      balancedaccuracy = mean(balancedaccuracy)
    )
  df.summary_sensitivity <- graphe_Sensitivity %>%
    group_by(method) %>%
    summarise(
      sd = sd(Sensitivity, na.rm = TRUE),
      Sensitivity = mean(Sensitivity)
    )
  df.summary_Specificity <- graphe_Specificity %>%
    group_by(method) %>%
    summarise(
      sd = sd(Specificity, na.rm = TRUE),
      Specificity = mean(Specificity)
    )
  
  
  
  accu <- ggplot(data=df.summary_accuracy, aes(x=method, y=balancedaccuracy, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=round(balancedaccuracy,2)), vjust=-0.25, color="black", position =  position_dodge(0.9) ,size=3.5)+
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() + ggtitle("Balanced accuracy") + 
    geom_errorbar(aes(ymin=balancedaccuracy-sd, ymax=balancedaccuracy+sd), width=.2, position=position_dodge(0.9)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + ylim(0, 1)
  
  sensi <- ggplot(data=df.summary_sensitivity, aes(x=method, y=Sensitivity, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=round(Sensitivity,2)), vjust=-0.25, color="black", position =  position_dodge(0.9) ,size=3.5)+
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() + ggtitle("Sensitivity") + 
    geom_errorbar(aes(ymin=Sensitivity-sd, ymax=Sensitivity+sd), width=.2, position=position_dodge(0.9)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + ylim(0,1)
  
  specif <- ggplot(data=df.summary_Specificity, aes(x=method, y=Specificity, fill=method)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=round(Specificity,2)), vjust=-0.25, color="black", position =  position_dodge(0.9) ,size=3.5)+
    scale_fill_brewer(palette = "Paired") +
    theme_minimal() + ggtitle("Specificity") + 
    geom_errorbar(aes(ymin=Specificity-sd, ymax=Specificity+sd), width=.2, position=position_dodge(0.9)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + ylim(0, 1)
  
  library(cowplot)
  p <- plot_grid(accu, sensi, specif, labels=c("A", "B", "C"), ncol = 3, nrow = 1 ) 
  
  title <- ggdraw() + 
    draw_label(
      "K-fold = 5, index-filter-merge, contigs exclusifs tumor, filter : nb patients >= 5, counts >= 10",
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) 
  
  plot_grid(
    title, p,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  
  print("balanced accuracy")
  print(df.summary_accuracy)
  print("sensitivity")
  print(df.summary_sensitivity)
  print("specificty")
  print(df.summary_Specificity)
  
  
