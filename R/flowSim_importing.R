#' Import files to analyze
#' 
#' function to import all datasets in the input directory
#' @param n_samples Number of files to analyze. If NULL, all files are analyzed Default to NULL.
#' @param paths_dir directory path of the input files to analyze
#' @param paths_plots Mandatory if paths_dir==NULL, it is possible to feed a list of the paths of the input files to analyze
#' @param n_cores The number of cores to use,default to 1
#' @return list of all input dfs
#' @export
#' @examples 
#' import_all_dfs(paths_dir="path_to_input_directory",paths_plots=NULL,n_cores=1)

import_all_dfs<-function(n_samples=NULL,paths_dir=NULL,paths_plots=NULL,n_cores=1){
  start<-Sys.time()
  if(is.null(paths_dir)==T && is.null(paths_plots)==T){
    stop("Path error: specify a folder or a list of paths")
  }else if(is.null(paths_dir)==F){
    paths_all_plots<-list.files(paths_dir,full.names = T,pattern = "*.csv",recursive = T)
  }else if(is.null(paths_dir)==T){
    paths_all_plots<-paths_plots
  }
  if(is.null(n_samples)==F){
    paths_all_plots<-paths_all_plots[n_samples]
  }
  print("--- all plot names extraction ----")
  list_all_input_plot_names<-mclapply(1:length(paths_all_plots),function(i){
    stringsplitted<-strsplit(paths_all_plots[i],"/")[[1]]
    plot_name<-tail(stringsplitted,1)
    return(c(plot_name,paths_all_plots[i]))
  },mc.cores = n_cores)
  df_names_vs_paths<-as.data.frame(do.call(rbind,list_all_input_plot_names),stringsAsFactors=F)
  colnames(df_names_vs_paths)<-c("name","path")
  print("-----Pre-allocation of df in list-----")
  vec_all_input_plot_names<-df_names_vs_paths$name
  list_all_input_dfs<-mclapply(setNames(vec_all_input_plot_names,vec_all_input_plot_names),function(name){
    ind_name<-which((vec_all_input_plot_names %in% name)==T)
    path_current_plot<-df_names_vs_paths$path[ind_name]
    if(length(path_current_plot)!=1){
      stop("path_current_plot has lenght != 1")
    }
    print(path_current_plot)
    df<-read.csv(path_current_plot)
    return(df)
  },mc.cores = n_cores)
  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
  return(list_all_input_dfs)
}


#' Import features set of all files
#' 
#' function to import the features set of all the files under analysis
#' @param n_samples number of files to analyze
#' @param paths_dir directory path of the input files to analyze
#' @param n_cores The number of cores to use, default to 1
#' @param progress_bar If False,progress bar is disabled. Default to True
#' @return a dataframe of the features set (columns) for each file (row)
#' @export
#' @examples 
#' import_df_features(n_samples=1:5,paths_dir="path_to_input_directory",n_cores=1)
import_df_features<-function(n_samples=NULL,paths_dir,n_cores=1,progress_bar=T){
  start<-Sys.time()
  print("--- get paths all plots -----")
  paths_all_plots<-list.files(paths_dir,full.names = T,pattern = "*.csv",recursive = T)
  if(is.null(n_samples)==F){
    paths_all_plots<-paths_all_plots[n_samples]
  }
  print("---------- features extraction from each plot ----------")
  if(progress_bar==T){
    cl <- makeSOCKcluster(n_cores)
    registerDoSNOW(cl)
    pb<-txtProgressBar(min=0,max=length(paths_all_plots),initial = 0,style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    list_all_plots_info<-foreach(i=1:length(paths_all_plots),.options.snow=opts) %dopar% {
      library(pracma)
      
      stringsplitted<-strsplit(paths_all_plots[i],"/")[[1]]
      plot_name<-tail(stringsplitted,1)
      df<-read.csv(paths_all_plots[i])
      if(is.numeric(df[,1])==F || is.numeric(df[,2])==F){
        return(NULL)
      }
      if(any(is.na(df[,1]))==T || any(is.na(df[,2]))==T){
        return(NULL)
      }
      if(length(df[,1])<6 || length(df[,2])<6){
        return(NULL)
        
      }
      
      density_m1<-(density(df[,1]))$y
      density_m2<-(density(df[,2]))$y
      density_m1<-density_m1/max(density_m1)
      density_m2<-density_m2/max(density_m2)
      # features extraction step
      mean_m1<-round(mean(df[,1]),2)
      mean_m2<-round(mean(df[,2]),2)
      matrix_peaks_m1<-findpeaks(density_m1,minpeakheight = 0.1,minpeakdistance = 50)
      n_peaks_m1<-nrow(matrix_peaks_m1)
      matrix_peaks_m2<-findpeaks(density_m2,minpeakheight = 0.1,minpeakdistance = 50)
      n_peaks_m2<-nrow(matrix_peaks_m2)
      min_peak_m1<-min(matrix_peaks_m1[,2])
      min_peak_m2<-min(matrix_peaks_m2[,2])
      max_peak_m1<-max(matrix_peaks_m1[,2])
      max_peak_m2<-max(matrix_peaks_m2[,2])
      h_min_peak_m1<-min(matrix_peaks_m1[,1])
      h_min_peak_m2<-min(matrix_peaks_m2[,1])
      h_max_peak_m1<-max(matrix_peaks_m1[,1])
      h_max_peak_m2<-max(matrix_peaks_m2[,1])
      final_vec<-c(plot_name,mean_m1,mean_m2,n_peaks_m1,n_peaks_m2,
                   min_peak_m1,min_peak_m2,max_peak_m1,max_peak_m2,
                   h_min_peak_m1,h_min_peak_m2,h_max_peak_m1,h_max_peak_m2)
      return(final_vec)
      
    }
    close(pb)
    stopCluster(cl)
  }else{
    list_all_plots_info<-mclapply(1:length(paths_all_plots),function(i){
      #print(i)
      stringsplitted<-strsplit(paths_all_plots[i],"/")[[1]]
      plot_name<-tail(stringsplitted,1)
      df<-read.csv(paths_all_plots[i])
      if(is.numeric(df[,1])==F || is.numeric(df[,2])==F){
        return(NULL)
      }
      if(any(is.na(df[,1]))==T || any(is.na(df[,2]))==T){
        return(NULL)
      }
      if(length(df[,1])<6 || length(df[,2])<6){
        return(NULL)
        
      }
      
      density_m1<-(density(df[,1]))$y
      density_m2<-(density(df[,2]))$y
      density_m1<-density_m1/max(density_m1)
      density_m2<-density_m2/max(density_m2)
      # features extraction step
      mean_m1<-round(mean(df[,1]),2)
      mean_m2<-round(mean(df[,2]),2)
      matrix_peaks_m1<-findpeaks(density_m1,minpeakheight = 0.1,minpeakdistance = 50)
      n_peaks_m1<-nrow(matrix_peaks_m1)
      matrix_peaks_m2<-findpeaks(density_m2,minpeakheight = 0.1,minpeakdistance = 50)
      n_peaks_m2<-nrow(matrix_peaks_m2)
      min_peak_m1<-min(matrix_peaks_m1[,2])
      min_peak_m2<-min(matrix_peaks_m2[,2])
      max_peak_m1<-max(matrix_peaks_m1[,2])
      max_peak_m2<-max(matrix_peaks_m2[,2])
      h_min_peak_m1<-min(matrix_peaks_m1[,1])
      h_min_peak_m2<-min(matrix_peaks_m2[,1])
      h_max_peak_m1<-max(matrix_peaks_m1[,1])
      h_max_peak_m2<-max(matrix_peaks_m2[,1])
      final_vec<-c(plot_name,mean_m1,mean_m2,n_peaks_m1,n_peaks_m2,
                   min_peak_m1,min_peak_m2,max_peak_m1,max_peak_m2,
                   h_min_peak_m1,h_min_peak_m2,h_max_peak_m1,h_max_peak_m2)
      return(final_vec)
    },mc.cores = n_cores)
  }
  print("---- combine all plots features in one df-------")
  df_features<-as.data.frame(do.call(rbind,list_all_plots_info),stringsAsFactors=F)
  colnames(df_features)<-c("plot_name","mean_m1","mean_m2",
                           "n_peaks_m1","n_peaks_m2","min_peak_m1","min_peak_m2","max_peak_m1",
                           "max_peak_m2","h_min_peak_m1","h_min_peak_m2","h_max_peak_m1","h_max_peak_m2")
  df_features$mean_m1<-as.numeric(df_features$mean_m1)
  df_features$mean_m2<-as.numeric(df_features$mean_m2)
  df_features$n_peaks_m1<-as.numeric(df_features$n_peaks_m1)
  df_features$n_peaks_m2<-as.numeric(df_features$n_peaks_m2)
  df_features$min_peak_m1<-as.numeric(df_features$min_peak_m1)
  df_features$min_peak_m2<-as.numeric(df_features$min_peak_m2)
  df_features$max_peak_m1<-as.numeric(df_features$max_peak_m1)
  df_features$max_peak_m2<-as.numeric(df_features$max_peak_m2)
  df_features$h_min_peak_m1<-as.numeric(df_features$h_min_peak_m1)
  df_features$h_min_peak_m2<-as.numeric(df_features$h_min_peak_m2)
  df_features$h_max_peak_m1<-as.numeric(df_features$h_max_peak_m1)
  df_features$h_max_peak_m2<-as.numeric(df_features$h_max_peak_m2)
  row.names(df_features)<-df_features$plot_name
  df_features<-df_features[,2:ncol(df_features)]
  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
  return(df_features)
}
