

get_similarity_selected_plot<-function(list_all_input_dfs,name_ref_plot,n_cores,
                                       names_analyzed_vec,n_boots,test="sm"){
  #---------- select expression reference plot name --------
  #print(sprintf("----- import expression reference plot selected: %s --------",name_ref_plot))
  ind<-which((names(list_all_input_dfs) %in% name_ref_plot)==T)
  if(length(ind)==0){
    stop("ref plot not found")
  }else if(length(ind)>1){
    print("ref plot:")
    print(name_ref_plot)
    print("duplication:")
    print(names(list_all_input_dfs)[ind])
    stop("ref plot multiple matches")
  }
  df_ref<-list_all_input_dfs[[ind]]
  #---------- get similarity scores -------------------
  #print("--------- get similarity scores ---------")
  list_scores<-mclapply(1:length(list_all_input_dfs),function(i){
    df<-list_all_input_dfs[[i]]
    plot_name_to<-names(list_all_input_dfs)[i]
    check_repetition<-plot_name_to %in% names_analyzed_vec
    if(check_repetition==F){
      # calculate distance between selected ref data and test 
      sink("/dev/null")
      on.exit(sink())
      scores_markers<-get_distance_loc_vs_test(loc_df = df_ref,test_df = df,nboot = n_boots)
      score_marker_1<-scores_markers$pvalue_mark_1
      score_marker_2<-scores_markers$pvalue_mark_2
      #------ report scores
      final_score<-(score_marker_1+score_marker_2)/2
      final_score<-round(final_score,2)
      final_vec<-c(name_ref_plot,plot_name_to,final_score)
    }else{
      final_vec<-NULL
    }
    return(final_vec)
  },mc.cores = n_cores)
  df_scores<-as.data.frame(do.call(rbind,list_scores))
  colnames(df_scores)<-c("from","to","weight")
  if(test=="js"){
    df_scores$weight<-as.numeric(as.character(df_scores$weight))
    df_scores$weight<- 1-(df_scores$weight/max(df_scores$weight))
    df_scores$weight<- round(df_scores$weight,2)
    
  }
  return(df_scores)
}


get_similarity_all_plots<-function(n_samples=NULL,path_dir=NULL,thr_score=0.6,
                                   n_cores=1,nboots=NULL,paths_plots=NULL){
  start<-Sys.time()
  #---------- Pre-allocation of input dfs-------------------
  print("------- Pre-allocation of data --------")
  list_all_input_dfs<-import_all_dfs(n_samples=n_samples,paths_dir = path_dir,n_cores = n_cores,paths_plots = paths_plots)
  #------------ Get all plot names ---------------
  # plot names will be  the nodes of our network
  print(sprintf("Total number of plots: %s",length(list_all_input_dfs)))
  plot_names<-names(list_all_input_dfs)
  #-------- set number of simulations
  if(length(list_all_input_dfs)<30){
    n_boots<-50
  }else if(length(list_all_input_dfs)<150){
    n_boots<-40
  }else if(length(list_all_input_dfs)<500){
    n_boots<-20
  }else{
    warning("Memory warning: Too many files, use flowSim batches approach: get_similarity_all_plots_v2",immediate. = T)
    n_boots<-20
  }
  if(is.null(nboots)==F){
    n_boots<-nboots
  }
  #---------- get similarity scores -------------------
  print("--------- get similarity scores all plots ---------")

  names_already_analyzed<-rep("a",length(plot_names))
  list_df_scores<-list()
  pb<- txtProgressBar(min=0,max=length(plot_names),initial = 0,style=3)
  for(i in 1:length(plot_names)){
    current_plot_name<-plot_names[i]
    # print("-------------reference name:")
    # print(current_plot_name)
    df_scores<-get_similarity_selected_plot(list_all_input_dfs = list_all_input_dfs,name_ref_plot = current_plot_name,
                                            n_cores = n_cores, names_analyzed_vec=names_already_analyzed,n_boots = n_boots)
    names_already_analyzed[i]<-current_plot_name
    list_df_scores[[i]]<-df_scores
    setTxtProgressBar(pb,i)
  }
  close(pb)
  gc()
  df_scores_all_plots<-do.call(rbind,list_df_scores)
  df_scores_all_plots$weight<-as.numeric(as.character(df_scores_all_plots$weight))
  df_scores_all_plots$from<-as.character(df_scores_all_plots$from)
  df_scores_all_plots$to<-as.character(df_scores_all_plots$to)
  # remove loops (connections to same node)
  print("----- remove loops")
  check_loop<-df_scores_all_plots$from == df_scores_all_plots$to
  inds<-which(check_loop==T)
  if(length(inds)!=0){
    df_scores_all_plots<-df_scores_all_plots[-inds,]
  }
  # remove edges with weight < thr_score
  if(thr_score!="all"){
    print("----- remove edges with weight <= threshold")
    inds<-which(df_scores_all_plots$weight<=thr_score)
    if(length(inds)!=0){
      df_scores_all_plots<-df_scores_all_plots[-inds,]
    }
  }
  # make nodes df
  print("----- make nodes df")
  df_nodes<-as.data.frame(plot_names)
  colnames(df_nodes)<-"id"
  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
  return(list(edges=df_scores_all_plots,nodes=df_nodes))
}
# kmeans function
myKmeans <- function(dist, k){
  return(kmeans(dist, k, iter.max = 50, nstart = 5)$cluster)
}
# function to get an approximate version of the similarity (for very big datasets)
get_similarity_all_plots_v2<-function(n_samples=NULL,path_dir,thr_score=0.9,
                                      n_cores=1,n_batches=NULL,show_tsne_plot=T,
                                      progress_bar=T,size_points_tsne=5,nboots=6){
  set.seed(40)
  start<-Sys.time()
  #---------- Features df -------------------
  print("------- Features df generation --------")
  df_features<-import_df_features(n_samples=n_samples,paths_dir =  path_dir,n_cores = n_cores,
                                  progress_bar = progress_bar)
  #------------ Get all plot names ---------------
  # plot names will be the nodes of our network
  print(sprintf("Total number of plots: %s",nrow(df_features)))
  #----- clustering the multidimensional feature set
  print("----- clustering on the multidimensional features set -----")
  if(nrow(df_features)<=3){
    stop("V2 mode can be executed only with >3 files")
  }else if(nrow(df_features)<=10){
    df_reduced<-as.data.frame((Rtsne(df_features,perplexity = 1,check_duplicates = F))$Y)
    print("perplexity: 1")
    max_clust_value<-nrow(df_features)-1
  }else if(nrow(df_features)<=30){
    warning("low number of files, use V1 mode")
    df_reduced<-as.data.frame((Rtsne(df_features,perplexity = 1,check_duplicates = F))$Y)
    print("perplexity: 1")
    max_clust_value<-10
  }else if(nrow(df_features)<=150){
    #warning("low number of files, use V1 mode")
    df_reduced<-as.data.frame((Rtsne(df_features,perplexity = 10,check_duplicates = F))$Y)
    print("perplexity: 10")
    max_clust_value<-30
  }else{
    df_reduced<-as.data.frame((Rtsne(df_features,perplexity = 30,check_duplicates = F))$Y)
    print("perplexity: 30")
    max_clust_value<-50
  }
  row.names(df_reduced)<-row.names(df_features)
  if(nrow(df_features)<=2000 && is.null(n_batches)==T){
    start_1<-Sys.time()
    print("estimate ideal number of batches")
    dist_m<-as.matrix(dist(df_reduced))
    n_batches<-nClust(meanDist = dist_m,maxClust = max_clust_value,clusteringFunction = myKmeans)
    end_1<-Sys.time()
    time_taken<-end_1-start_1
    print("Execution time:")
    print(time_taken)
    print("Done")
  }else if(nrow(df_features)>2000 && is.null(n_batches)==T){
    ggplot_no_clusters<-ggplot(df_reduced,aes(x=V1,y=V2)) + geom_point()
    ggplot_no_clusters<- ggplot_no_clusters + theme(legend.position = "none")
    show(ggplot_no_clusters)
    n_batches<-readline(prompt = "n_batches = NULL. Please, enter number of batches: ")
    n_batches<-as.integer(n_batches)
  }else if(is.null(n_batches)==F){
    n_batches<-n_batches
  }
  out_clustering<-kmeans(x = df_reduced,centers = n_batches)
  clusters<-out_clustering$cluster
  final_df<-as.data.frame(cbind(df_reduced,clusters),stringsAsFactor=F)
  #------ visualize clustered reduced  df ------
  ggplot_clusters<-ggplot(final_df,aes(x=V1,y=V2)) + geom_point(aes(colour=factor(final_df$clusters)),size=size_points_tsne)
  ggplot_clusters<- ggplot_clusters + theme(legend.position = "none")
  ggplot_clusters<-ggplot_clusters+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  ggplot_clusters<- ggplot_clusters + theme(axis.text=element_text(size=18),axis.title.x = element_text(size = 20,face="bold"),axis.title.y = element_text(size=20,face="bold"))
  
  if(show_tsne_plot==T){
    show(ggplot_clusters)
  }
  #------ get all plots names and paths input dir -----
  print("---- get all plots names input dir ----")
  paths_all_plots<-list.files(path_dir,full.names = T,pattern = "*.csv",recursive = T)
  all_plots_names<-sapply(1:length(paths_all_plots),function(i){
    path_current_plot<-paths_all_plots[i]
    stringsplitted<-strsplit(path_current_plot,"/")[[1]]
    plot_name<-tail(stringsplitted,1)
    return(plot_name)
  })
  #----- analyze each clusters separately
  print("------- analyze each clusters separately ------")
  all_groups<-unique(final_df$clusters)
  n_groups<-length(all_groups)
  print(paste0("n_groups: ",n_groups))
  list_results<-batches_analysis(all_groups = all_groups,n_cores = n_cores,final_df = final_df,
                                 path_dir=path_dir,
                                 progress_bar = progress_bar,thr_score = thr_score,nboots=nboots)
  print("--------- extract list results ------")
  list_all_df_scores<-lapply(1:length(list_results),function(i){
    list_result_i<-list_results[[i]]
    df_scores_all_plots<-list_result_i$df_scores_all_plots
    return(df_scores_all_plots)
  })
  list_all_nodes_df<-lapply(1:length(list_results),function(i){
    list_result_i<-list_results[[i]]
    df_nodes<-list_result_i$df_nodes
    return(df_nodes)
  })
  print("------- combines all dfs -----")
  final_df_scores<-do.call(rbind,list_all_df_scores)
  final_df_nodes<-do.call(rbind,list_all_nodes_df)
  end<-Sys.time()
  time_taken<-end-start
  print("Execution time:")
  print(time_taken)
  print("Done")
  return(list(edges=final_df_scores,nodes=final_df_nodes,df_features=df_features))
}

# function to analyze each batch
batches_analysis<-function(all_groups,n_cores,final_df,path_dir,thr_score,progress_bar=T,nboots){
  splitted_path_dir<-strsplit(path_dir,"/")[[1]]
  path_dir<-paste0(splitted_path_dir,collapse = "/")
  if(progress_bar==T){
    cl <- makeSOCKcluster(n_cores)
    registerDoSNOW(cl)
    pb<-txtProgressBar(min=0,max=length(all_groups),initial = 0,style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    list_results<-foreach(i=1:length(all_groups),.options.snow=opts) %dopar% {
      library(sm)
      library(parallel)
      library(igraph)
      library(visNetwork)
      library(RColorBrewer)
      library(utils)
      library(pracma)
      library(stats)
      library(Rtsne)
      library(ggplot2)
      source("/home/rstudio/code/flowSim_similarity.R")
      source("/home/rstudio/code/flowSim_utils.R")
      source("/home/rstudio/code/flowSim_visnet_conversion.R")
      source("/home/rstudio/code/flowSim_visualization.R")
      source("/home/rstudio/code/flowSim_igraph_clustering.R")
      source("/home/rstudio/code/flowSim_exporting.R")
      source("/home/rstudio/code/flowSim_importing.R")
      current_group<-all_groups[i]
      #--- select samples name current group
      inds<-which(final_df$clusters==current_group)
      final_df_current_group<-final_df[inds,]
      samples_name_current_group<-row.names(final_df_current_group)
      # ---- get paths current samples names
      paths_plots_current_group<-sapply(1:length(samples_name_current_group),function(x){
        current_name<-samples_name_current_group[x]
        path_current_name<-paste0(path_dir,"/",current_name)
        return(path_current_name)
      })
      #----- import list dfs current samples
      list_all_input_dfs_current_group<-import_all_dfs(paths_plots = paths_plots_current_group,
                                                       n_cores = n_cores)
      #------------ Get all plot names ---------------
      # plot names will be  the nodes of our network
      plot_names<-names(list_all_input_dfs_current_group)
      #---------- get similarity scores -------------------
      names_already_analyzed<-rep("a",length(plot_names))
      list_df_scores<-list()
      for(j in 1:length(plot_names)){
        current_plot_name<-plot_names[j]
        df_scores<-get_similarity_selected_plot(list_all_input_dfs = list_all_input_dfs_current_group,name_ref_plot = current_plot_name,
                                                n_cores = n_cores, names_analyzed_vec=names_already_analyzed,n_boots = nboots)
        names_already_analyzed[j]<-current_plot_name
        list_df_scores[[j]]<-df_scores
      }
      gc()
      df_scores_all_plots<-do.call(rbind,list_df_scores)
      df_scores_all_plots$weight<-as.numeric(as.character(df_scores_all_plots$weight))
      df_scores_all_plots$from<-as.character(df_scores_all_plots$from)
      df_scores_all_plots$to<-as.character(df_scores_all_plots$to)
      # remove loops (connections to same node)
      check_loop<-df_scores_all_plots$from == df_scores_all_plots$to
      inds<-which(check_loop==T)
      if(length(inds)!=0){
        df_scores_all_plots<-df_scores_all_plots[-inds,]
      }
      # remove edges with weight < thr_score
      if(thr_score!="all"){
        inds<-which(df_scores_all_plots$weight<=thr_score)
        if(length(inds)!=0){
          df_scores_all_plots<-df_scores_all_plots[-inds,]
        }
      }
      # make nodes df
      df_nodes<-as.data.frame(plot_names)
      colnames(df_nodes)<-"id"
      # update lists
      return(list(df_scores_all_plots=df_scores_all_plots,df_nodes=df_nodes))
    }
    stopCluster(cl)
    close(pb)
  }else if(progress_bar==F){
    list_results<-mclapply(1:length(all_groups),function(i){
      current_group<-all_groups[i]
      #--- select samples name current group
      inds<-which(final_df$clusters==current_group)
      final_df_current_group<-final_df[inds,]
      samples_name_current_group<-row.names(final_df_current_group)
      # ---- get paths current samples names
      paths_plots_current_group<-sapply(1:length(samples_name_current_group),function(x){
        current_name<-samples_name_current_group[x]
        path_current_name<-paste0(path_dir,"/",current_name)
        return(path_current_name)
      })
      #----- import list dfs current samples
      list_all_input_dfs_current_group<-import_all_dfs(paths_plots = paths_plots_current_group,
                                                       n_cores = n_cores)
      #------------ Get all plot names ---------------
      # plot names will be  the nodes of our network
      plot_names<-names(list_all_input_dfs_current_group)
      #---------- get similarity scores -------------------
      names_already_analyzed<-rep("a",length(plot_names))
      list_df_scores<-list()
      for(j in 1:length(plot_names)){
        current_plot_name<-plot_names[j]
        df_scores<-get_similarity_selected_plot(list_all_input_dfs = list_all_input_dfs_current_group,name_ref_plot = current_plot_name,
                                                n_cores = n_cores, names_analyzed_vec=names_already_analyzed,n_boots = nboots)
        names_already_analyzed[j]<-current_plot_name
        list_df_scores[[j]]<-df_scores
      }
      gc()
      df_scores_all_plots<-do.call(rbind,list_df_scores)
      df_scores_all_plots$weight<-as.numeric(as.character(df_scores_all_plots$weight))
      df_scores_all_plots$from<-as.character(df_scores_all_plots$from)
      df_scores_all_plots$to<-as.character(df_scores_all_plots$to)
      # remove loops (connections to same node)
      check_loop<-df_scores_all_plots$from == df_scores_all_plots$to
      inds<-which(check_loop==T)
      if(length(inds)!=0){
        df_scores_all_plots<-df_scores_all_plots[-inds,]
      }
      # remove edges with weight < thr_score
      inds<-which(df_scores_all_plots$weight<=thr_score)
      if(length(inds)!=0){
        df_scores_all_plots<-df_scores_all_plots[-inds,]
      }
      # make nodes df
      df_nodes<-as.data.frame(plot_names)
      colnames(df_nodes)<-"id"
      # update lists
      return(list(df_scores_all_plots=df_scores_all_plots,df_nodes=df_nodes))
    },mc.cores = n_cores)
  }
  print("batches analysis Done")
  return(list_results)
}
