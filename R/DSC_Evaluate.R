#' DSC_Evaluate
#' @export
## SFMiC - Evaluation
DSC_Evaluate <- function(){

  if(!is.null(description)) desc <- description
  else desc <-"Multiclass Fuzzy Nolvelty Detector- Results"
  DSC_Evaluate <- fminaseval_refClass$new()
  structure(list(description = desc, RObj = DSC_Evaluate), class = c("DSC_R","DSC"))
}

fminaseval_refClass <-
  setRefClass("fuzzyMinasEval",
              fields = list(
                ## args
                accuracy = "numeric",
                macro_fscore = "numeric",
                cer = "numeric",
                unknown_rate = "numeric",
                pn = "numeric",
                allconfusion_matrix = "list",
                allconfusion_matrix_exp = "list"
              ),

              methods = list(
                initialize = function() {
                  accuracy <<- numeric()
                  macro_fscore <<- numeric()
                  cer <<- numeric()
                  pn <<- numeric()
                  unknown_rate <<- numeric()
                  allconfusion_matrix <<- list()
                  allconfusion_matrix_exp <<- list()
                .self
                }
              )
  )

fminaseval_refClass$methods(
  eval = function(model){
   for(moment in 1:length(model$RObj$evaluation_hist)){
      #================Calculating confusin matrix===============
      y_explained <- model$RObj$evaluation_hist[[moment]]$predictions
      init <- initial_npoints+1
      end <- initial_npoints+length(y_explained)
      ground_truth <- dataset[init:end,ncol(dataset)]
      ground_truth[is.na(ground_truth)] <- "unknown"
      confusion_matrix <- table(ground_truth,y_explained)
      cm_colnames <- colnames(confusion_matrix)
      cm_rowlnames <- rownames(confusion_matrix)
      pns_index <- which(str_detect(colnames(confusion_matrix),"PN")==TRUE)
      index_class <- sapply(pns_index, function(index){which.max(confusion_matrix[,index])})

      foreach(pn_index = pns_index, cl = index_class) %do%{
        y_explained[which(y_explained == cm_colnames[pn_index])] <- cm_rowlnames[cl]
      }

      levels <- sort(c(unique(ground_truth),"unknown"))

      confusion_matrix <- table("Actual" = factor(ground_truth,levels = levels), "Predicted" = factor(y_explained, levels = levels))
      confusion_matrix_exp <- confusion_matrix[1:(ncol(confusion_matrix)-1),-ncol(confusion_matrix)]

      allconfusion_matrix <<- c(allconfusion_matrix, list(confusion_matrix))
      allconfusion_matrix_exp <<- c(allconfusion_matrix_exp, list(confusion_matrix_exp))
      #===========================================================

      #=====================Calculating Metrics===================
      n = sum(confusion_matrix_exp) # number of instances
      nc = nrow(confusion_matrix_exp) # number of classes
      diag = diag(confusion_matrix_exp) # number of correctly classified instances per class
      rowsums = apply(confusion_matrix_exp, 1, sum) # number of instances per class
      colsums = apply(confusion_matrix_exp, 2, sum) # number of predictions per class
      p = rowsums / n # distribution of instances over the actual classes
      q = colsums / n # distribution of instances over the predicted classes

      #======================Overal Accuracy==============================
      acc = sum(diag) / n
      accuracy <<- c(accuracy, acc)
      #================Precision, Recall, F1 by Class=====================

      precision = diag / colsums
      precision[is.nan(precision)] <- 1
      recall = diag / rowsums
      recall[is.nan(recall)] <- 1
      f1 = 2 * precision * recall / (precision + recall)

      #===============Macro Precision, Recal, F1 ========================
      macroPrecision = mean(precision)
      macroRecall = mean(recall)
      macroF1 = mean(f1)
      macro_fscore <<- c(macro_fscore,macroF1)
      #===================CER============================================
      dumb_cm <- confusion_matrix_exp
      diag(dumb_cm) <- 0
      FN <- apply(dumb_cm, 1, sum) # number of instances per class
      FP <- apply(dumb_cm, 2, sum) # number of predictions per class

      FNR <- FN/rowsums
      FNR[is.nan(FNR)] <- 0
      FPR <- FP/colsums
      FPR[is.nan(FPR)] <- 0
      rate <- rowsums/n

      eval_cer <- sum(rate*FNR + rate*FPR)/2
      cer <<- c(cer,eval_cer)

      #=================UnkR==============================================
      unk <- confusion_matrix[-(nc+1),(nc+1)]
      exc <- apply(confusion_matrix, 1, sum)
      eval_unkr <- sum(unk/exc[-(nc+1)])/nc
      unknown_rate <<- c(unknown_rate, eval_unkr)

      #=====================PN===========================================
      pn <<- c(pn,model$RObj$evaluation_hist[[moment]]$PN)

    }
  },

  get_accuracy = function(...) {accuracy}
)
