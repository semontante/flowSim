 

#' get distance selected file vs another local file
#' 
#' function to compare the similiraty between selected file and 
#' another file of the input directory based on the density of marker expression values. Internal function.
#' @param test_df Dataframe of selected file
#' @param loc_df Dataframe of another file of the input directory
#' @param show_plot show plot of the density comparison
#' @param nboot Number of simulations for the permuation identity test between the densities.
#' @keywords Internal
#' @return A list of two float numbers
#' $pvalue_mark_1 : pvalue first marker comparison
#' $pvalue_mark_2 : pvalue second marker comparison
get_distance_loc_vs_test<-function(test_df,loc_df,show_plot="none",nboot=50){
  set.seed(123)
  # get markers expression for test and train
  marker_1_expr_test<-test_df[,1]
  marker_2_expr_test<-test_df[,2]
  marker_1_expr_train<-loc_df[,1]
  marker_2_expr_train<-loc_df[,2]
  # get density of markers for test and train
  d1_test<-density(marker_1_expr_test,n=50)
  d2_test<-density(marker_2_expr_test,n=50)
  d1_train<-density(marker_1_expr_train,n=50)
  d2_train<-density(marker_2_expr_train,n=50)
  # get values correctly formatted
  vec_1<-rep(1,length(d1_train$y))
  vec_2<-rep(2,length(d1_test$y))
  vec_all_mark_1<-c(vec_1,vec_2)
  d_all_mark_1<-c(d1_train$y,d1_test$y)
  vec_1<-rep(1,length(d2_train$y))
  vec_2<-rep(2,length(d2_test$y))
  vec_all_mark_2<-c(vec_1,vec_2)
  d_all_mark_2<-c(d2_train$y,d2_test$y)
  # Execute significance test between the density for each marker
  output_mark_1<-sm.density.compare(d_all_mark_1,group = vec_all_mark_1,model = "equal",display=show_plot,nboot=nboot)
  output_mark_2<-sm.density.compare(d_all_mark_2,group = vec_all_mark_2,model = "equal",display=show_plot,nboot=nboot)
  pvalue_mark_1<-output_mark_1$p
  pvalue_mark_2<-output_mark_2$p
  # lower pvalue more significant differences between the two densities
  return(list(pvalue_mark_1=pvalue_mark_1,pvalue_mark_2=pvalue_mark_2))
  
}

# 
#' Get heterogeneity score
#' 
#' function to get the heterogeneity score of a dataset
#' @param igraph_network Igraph network object
#' @param partition_vec Vectors of files and clusters associations
#' @param visnetdata list of dataframes containing the edges and nodes information in VisNetwork format
#' @return A single float number indicating the final heterogeneity score
#' @export
#' @examples 
#' \donttest{gen_igraph_network(igraph_network=igraph_network,partition_vec=partition_vec,visnetdata=visnetdata)}
get_heterogeneity_score<-function(igraph_network,partition_vec,visnetdata){
  print("----- get nodes and edges info -----")
  n_edges<-nrow(visnetdata$edges)
  tot_n_clusters<-length(unique(visnetdata$nodes$group))
  n_nodes<-length(visnetdata$nodes$id)
  all_groups<-unique(visnetdata$nodes$group)
  print(paste0("tot_n_clusters: ",tot_n_clusters))
  print(paste0("n_edges: ",n_edges))
  
  print("----- get clusters info -----")
  counter_single_clusters<-0
  counter_multiple_clusters<-0
  for(i in 1:length(all_groups)){
    current_group<-all_groups[i]
    check_group<-visnetdata$nodes$group %in% current_group
    inds<-which(check_group==T)
    nodes_current_group<-visnetdata$nodes$id[inds]
    n_nodes_current_group<-length(nodes_current_group)
    if(n_nodes_current_group==1){
      counter_single_clusters<-counter_single_clusters+1
    }else{
      counter_multiple_clusters<-counter_multiple_clusters+1
    }
  }
  n_clusters_single_node<-counter_single_clusters
  n_clusters_multiple_nodes<-counter_multiple_clusters
  print(paste0("n_clusters_single_node: ",n_clusters_single_node))
  print(paste0("n_clusters_multiple_nodes: ",n_clusters_multiple_nodes))
  print("----- calculate score -----")
  edges_prop<-n_edges/n_nodes
  print(paste0("edges fraction: ",edges_prop))
  # if the number of edges is too low 
  # we calculate the heterogeneity score using fraction of single files
  score_1<-n_clusters_single_node/tot_n_clusters
  score_1<-round(score_1,2)
  print(sprintf("fraction score: %s",score_1))
  # we use modularity formula
  score_2<-modularity(igraph_network,membership = partition_vec)
  score_2<-round(score_2,2)
  print(sprintf("modularity score: %s",score_2))
  final_score<- (score_1 + score_2)/2
  if(edges_prop<=0.01){
    final_score<- score_1 
  }
  final_score<-round(final_score,2)
  return(final_score)
}


#' get_new_colors_based_on_groups
#' 
#' function to get new colors based on groups. Internal function
#' @param df_groups Dataframe containing the groups associations
#' @keywords Internal
#' @return A dataframe containing the association between groups,files and colors
get_new_colors_based_on_groups<-function(df_groups){
  numb_clusters<-length(unique(df_groups[,2]))
  # we generate vec of unique colors
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  #col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector<-brewer.pal(12, "Paired")
  col_vector <- colorRampPalette(col_vector)(numb_clusters)
  # extract random colors
  unique_groups<-unique(df_groups[,2])
  rand_int<-sample.int(length(col_vector),length(unique_groups))
  colors_selected<-col_vector[rand_int]
  df_temp<-as.data.frame(cbind(unique_groups,colors_selected))
  df_temp$unique_groups<-as.character(df_temp$unique_groups)
  final_col_vec<-rep("a",nrow(df_groups))
  df_temp$colors_selected<-as.character(df_temp$colors_selected)
  # assign colors to the clusters
  for(i in 1:nrow(df_temp)){
    current_group<-df_temp$unique_groups[i]
    current_color<-df_temp$colors_selected[i]
    inds<-which(df_groups[,2]==current_group)
    final_col_vec[inds]<-current_color
  }
  df_groups_colors<-as.data.frame(cbind(df_groups,final_col_vec),stringsAsFactors=F)
  return(df_groups_colors)
}


#' get_size_based_on_groups
#' 
#' function to get size of groups based on groups. Internal function
#' @param df_groups Dataframe containing the groups associations
#' @keywords Internal
#' @return A dataframe containing the association between groups,files and size
get_size_based_on_groups<-function(df_groups){
  n_groups<-length(unique(df_groups[,2]))
  # we generate vec of sizes
  final_size_vec<-rep(0,nrow(df_groups))
  unique_groups<-unique(df_groups[,2])
  for(i in 1:length(unique_groups)){
    current_group<-unique_groups[i]
    inds<-which(df_groups[,2]==current_group)
    df_groups_current_group<-df_groups[inds,]
    size_current_group<-nrow(df_groups_current_group)
    final_size_vec[inds]<-size_current_group
  }
  return(final_size_vec)
}


 
#' get_hull_all_gates
#' 
#' function to get the convex hull of all gates. Internal function
#' @param df_groups Dataframe containing the groups associations
#' @keywords Internal
#' @return A list containing the convex hull for all classes
get_hull_all_gates<-function(gated_df){
  colnames(gated_df)<-c("x","y","classes")
  all_classes<-unique(gated_df$classes)
  list_df_hull<-list()
  for(classes in all_classes){
    inds<-which(gated_df$classes==classes)
    df_current_classes<-gated_df[inds,]
    if(nrow(df_current_classes)>200000){ # if there are many events I use the concave function
      df_current_classes_hull_values<-as.data.frame(concaveman(as.matrix(df_current_classes[,c(1,2)]),concavity=5))
    }else{ # otherwise I use the classical convex hull
      current_hull_indices<-chull(df_current_classes[,c(1,2)]) # indices points at the border of the convex hull for the current classes
      df_current_classes_hull_values<-df_current_classes[current_hull_indices,c(1,2)]
    }
    vec_group<-rep(sprintf("%s",classes),nrow(df_current_classes_hull_values))
    df_current_hull<-cbind(df_current_classes_hull_values,vec_group)
    colnames(df_current_hull)<-c("x","y","group_gate")
    list_df_hull[[sprintf("%s",classes)]]<-df_current_hull
  }
  return(list_df_hull)
}

#' get_inds_files_selected
#' 
#' function to get indices for the paths of the files selected
#' @param path_dir path of the directory to explore
#' @param files_selected Vector containing the names of the files whose indices need to be found
#' @param return_path If True, return paths instead of indices. Default to False
#' @return Vector containing the indices of the files
#' @export
#' @examples 
#' \donttest{get_inds_files_selected(path_dir="path/to/directory",files_selected=vector_names)}

get_inds_files_selected<-function(path_dir,files_selected,return_path=F){
  paths_all_plots<-list.files(path_dir,full.names = T,pattern = "*.csv",recursive = T)
  splitted_list<-strsplit(paths_all_plots,"/")
  all_files_names<-sapply(splitted_list, function(x){
    return(tail(x,1))
  })
  inds_selected_files<-which((all_files_names %in% files_selected)==T)
  if(length(inds_selected_files) != length(files_selected)){
    warning("some files selected not found")
  }
  if(return_path==F){
    return(inds_selected_files)
  }else if(return_path==T){
    paths_selected<-paths_all_plots[inds_selected_files]
    return(paths_selected)
  }
}

# 

#' rename_files
#' 
#' Function to rename files within a directory
#' @param path_dir path of the directory with files to rename
#' @param to_replace Pattern to replace. if to_replace_fixed=F, regex expression can be used. 
#' @param string string to replace the pattern selected
#' @param to_replace_fixed if False, to_replace can be also a regex expression. Default to True.
#' @return NULL
#' @export
#' @examples 
#' \donttest{rename_files(path_dir="path/to/directory",to_replace="pattern to replace",string="string replaced",regex_expr=F)}

rename_files<-function(path_dir,to_replace,string,to_replace_fixed=T){
  old_paths<-list.files(path_dir,full.names = T,recursive = T)
  sapply(1:length(old_paths),function(i){
    old_path_i<-old_paths[i]
    new_path_i<-sub(to_replace,string,old_path_i,fixed = to_replace_fixed)
    file.rename(old_path_i,new_path_i)
  })
}

#' get_csv
#' 
#' function to convert fcs files in bivariate csv files
#' @param path path of the directory with fcs files to convert
#' @param channels Vector reporting the two channels to consider
#' @param markers Vector reporting the two markers to consider. Considered only if channels==NULL
#' @param path.output directory path to export the converted csv files
#' @return NULL
#' @export
#' @examples 
#' \donttest{get_csv(path="path/to/directory",channels=c("FSC-A","FSC-H"),path.output="path/to/directory")}

get_csv<-function(path,channels=NULL,markers=NULL,path.output){
  if(is.null(channels)==T & is.null(markers)==T){
    stop("channels and markers argument are both NULL")
  }
  FCSfiles_path_list <- list.files(path = path, recursive = F, pattern = ".fcs", full.names = T)
  fcs_files_list <-lapply(FCSfiles_path_list,read.FCS) 
  fs <- as(fcs_files_list,"flowSet")
  all_markers<-fs[[1]]@parameters@data$desc
  all_channels<-fs[[1]]@parameters@data$name
  print("All channels:")
  print(all_channels)
  print("All markers:")
  print(all_markers)
  output<-sapply(1:length(fs),function(i){
    frame<-fs[[i]]
    if(is.null(channels)==T){
      channels<-c("a","b")
      ind_1_marker<-which(as.vector(frame@parameters@data$desc)==markers[1])
      if(length(ind_1_marker)==0){
        stop("marker 1 not found")
      }
      channels[1]<-frame@parameters@data$name[ind_1_marker]
      ind_2_marker<-which(as.vector(frame@parameters@data$desc)==markers[2])
      if(length(ind_2_marker)==0){
        stop("marker 2 not found")
      }
      channels[2]<-frame@parameters@data$name[ind_2_marker]
    }
    df_data<-exprs(frame)[,channels]
    name_file<-tail(strsplit(FCSfiles_path_list[i],"/")[[1]],1)
    name_file<-gsub(name_file,pattern = ".fcs",replacement = ".csv")
    write.table(df_data,file=paste0(path.output,"/",name_file), row.names=FALSE, sep=",")
  })
}

