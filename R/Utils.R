
#=======Create experiment folder======================
#' create dataset folder
#' @export
create_dsfolder = function(dataset_name){
  dsname <- gsub(".csv","",dataset_name)
  dir.create("Results", showWarnings = FALSE)
  results_path = paste0("Results/",dsname)
  dir.create(results_path, showWarnings = FALSE)
  return(results_path)
}

#==============Check Experiment Number================
#' create experiment folder
#' @export
create_expfolder = function(ds_path){
  experiments <- list.files(ds_path)
  ds_path <- paste0(ds_path,"/")
  exp_range <- sapply(experiments, function(name){gsub("Exp","",name)})
  new_exp <- length(exp_range) + 1
  exp_folder <- paste0("EXP",new_exp)
  exp_folder <- paste0(ds_path,exp_folder)
  dir.create(exp_folder, showWarnings = FALSE)
  return(exp_folder)
}

#=================Printing Parameters=================
#' Save model settings
#' @import gridExtra
#' @export
save_model = function(model, path){
  param_table <- data.frame(initial_npoints,
                            initial_k,
                            model$RObj$maxMiC,
                            model$RObj$minNumSFMiC,
                            model$RObj$m,
                            model$RObj$theta,
                            model$RObj$thetaClass,
                            model$RObj$timeWindow,
                            model$RObj$P,
                            model$RObj$ts,
                            model$RObj$thresholdTemp,
                            model$RObj$k_temp,
                            model$RObj$silhouetteThreshold,
                            model$RObj$fr_threshold,
                            model$RObj$minWeight)
  colnames(param_table) <- c("#Initial points",
                             "#Initial k",
                             "#Max SFMiCS by Class",
                             "#Min SFMiCs total",
                             "Fuzzification M",
                             "Theta",
                             "Theta_class",
                             "Time window",
                             "Sleep threshold P",
                             "Timestamp threshold ts",
                             "#Temporary threshold",
                             "#Temporary k",
                             "silhouette threshold",
                             "FR_Threshold",
                             "#Min Weight NP"
  )
  rownames(param_table) <- "Value"
  param_table <- t(param_table)

  path <- paste0(path,"/Parameters.pdf")
  pdf(path)
  grid.table(param_table)
  dev.off()
}

calculate_overall_cer = function(results_by_execution){
  return(rowMeans(sapply(results_by_execution,function(exec){exec$RObj$cer})))
}


calculate_overall_accuracy = function(results_by_execution){
  return(rowMeans(sapply(results_by_execution,function(exec){exec$RObj$accuracy})))
}

calculate_overall_f1 = function(results_by_execution){
  return(rowMeans(sapply(results_by_execution,function(exec){exec$RObj$macro_fscore})))
}

calculate_overall_unkr = function(results_by_execution){
  return(rowMeans(sapply(results_by_execution,function(exec){exec$RObj$unknown_rate})))
}

calculate_overall_pn = function(results_by_execution){
  return(rowMeans(sapply(results_by_execution,function(exec){exec$RObj$pn})))
}

calculate_sd_cer = function(results_by_execution){
  return(rowSds(sapply(results_by_execution,function(exec){exec$RObj$cer})))
}

calculate_sd_accuracy = function(results_by_execution){
  return(rowSds(sapply(results_by_execution,function(exec){exec$RObj$accuracy})))
}

calculate_sd_f1 = function(results_by_execution){
  return(rowSds(sapply(results_by_execution,function(exec){exec$RObj$macro_fscore})))
}

calculate_sd_unkr = function(results_by_execution){
  return(rowSds(sapply(results_by_execution,function(exec){exec$RObj$unknown_rate})))
}

calculate_sd_pn = function(results_by_execution){
  return(rowSds(sapply(results_by_execution,function(exec){exec$RObj$pn})))
}


plot_metric = function(metric_values,sd,title_name,color,type){
  metric_values <- data.frame(values = metric_values,moments = c(1:length(metric_values)), sd = sd)
  pd <- position_dodge(0.1)
  title <- ggtitle(title_name)
  axis_name <- labs(x="Evaluation Moment", y="Values")
  theme <- theme_classic()+theme(plot.title = element_text(hjust = 0.5))
  line <- geom_line(linetype = type, colour = color)
  point <- geom_point(shape = 2, colour = color)
  sd <- geom_ribbon(aes(ymin=values-sd, ymax=values+sd), fill = "grey70")
  ylim <- ylim(0,1)
  print(ggplot(data = metric_values, aes(y = values, x = moments))+line+theme+title+axis_name+ylim)
}

plot_2metrics = function(metric1,sd1,metric2,sd2,title_name1,title_name2){
  title_name <- paste0(title_name1," x ",title_name2)
  title <- ggtitle(title_name)
  theme <- theme_classic()+theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(),legend.position='bottom')
  axis_name <- labs(x="Evaluation Moment", y="values")
  metric1 <- data.frame(values = metric1,moments = c(1:length(metric1)), Metrics = rep(title_name1,length(metric1)), sd = sd1)
  metric2 <- data.frame(values = metric2,moments = c(1:length(metric2)), Metrics = rep(title_name2,length(metric2)), sd = sd2)
  all_metrics <- rbind(metric1,metric2)

  all_metrics$Metrics <- factor(all_metrics$Metrics)
  line <- geom_line()
  point <- geom_point()
  pd <- position_dodge(0.1)
  sd <- geom_ribbon(aes(ymin=values-sd, ymax=values+sd), fill = "grey70")
  ylim <- ylim(0,1)
  print(ggplot(data = all_metrics, aes(y = values, x = moments, color = Metrics,linetype = Metrics, shape = Metrics))+line+theme+title+axis_name+ylim)
}


#' Plot all metrics
#' @import ggplot2 matrixStats
#' @export
plot_metrics = function(results, path){

  cer <- calculate_overall_cer(results)
  f1 <- calculate_overall_f1(results)
  unkr <- calculate_overall_unkr(results)
  pn <- calculate_overall_pn(results)
  acc <- calculate_overall_accuracy(results)

  sd_cer <- calculate_sd_cer(results)
  sd_f1 <- calculate_sd_f1(results)
  sd_unkr <- calculate_sd_unkr(results)
  sd_pn <- calculate_sd_pn(results)
  sd_acc <- calculate_sd_accuracy(results)

  overall_metrics <- data.frame(Metric = c(mean(cer),mean(f1), mean(unkr), pn[length(pn)], mean(acc)), SD = c(mean(sd_cer), mean(sd_f1), mean(sd_unkr), sd_pn[length(sd_pn)], mean(sd_acc)))
  rownames(overall_metrics) <- c("CER","F1","Unkr","#PN","Accuracy")

  path_overall <- paste0(path,"/Overall_Metrics+SD.pdf")
  pdf(path_overall)
  grid.table(overall_metrics)
  dev.off()

  path_plots <- paste0(path,"/metrics_plots.pdf")
  pdf(path_plots)

  plot_metric(metric_values = cer, sd = sd_cer, color = "darkblue", type = 1, title_name = "CER")
  plot_metric(metric_values = f1, sd = sd_f1, color = "darkblue", type = 1, title_name = "Macro F-Score")
  plot_metric(metric_values = unkr, sd = sd_unkr, color = "darkblue", type = 1, title_name = "UnkR")
  plot_metric(metric_values = acc, sd = sd_acc, color = "darkblue", type = 1, title_name = "Accuracy")
  plot_metric(metric_values = pn, sd = sd_pn, color = "darkblue", type = 1, title_name = "#NP")

  plot_2metrics(metric1 = cer, sd1 = sd_cer, metric2 = unkr, sd2 = sd_unkr, title_name1 = "CER", title_name2 = "UnkR")
  plot_2metrics(metric1 = f1, sd1 = sd_f1, metric2 = unkr, sd2 = sd_unkr, title_name1 = "Macro F-Score", title_name2 = "UnkR")
  plot_2metrics(metric1 = acc, sd1 = sd_acc, metric2 = unkr, sd2 = sd_unkr, title_name1 = "Accuracy", title_name2 = "UnkR")


  dev.off()


}

