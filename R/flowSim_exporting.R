 
#' exports_plots
#' 
#' function to export all plots for each group
#' @param visnetdata list of dataframes containing the edges and nodes information in VisNetwork format
#' @param path_expr_data Path to directory containing the expression data of the files analyzed
#' @param path_output Path to the directory where to export the plots for each group
#' If path_output=NULL, no plots are generated. Default to NULL
#' @param n_cores Number of cores to use. Default to 1
#' @param path_gates Path of gates data, default to NULL
#' @return A list of all file names for each group.
#' @export
#' @examples 
#' \donttest{exports_plots(visnetdata=visnetdata,path_expr_data="path to input directory",path_output=NULL,n_cores=1)}
exports_plots<-function(visnetdata,path_expr_data,path_output=NULL,n_cores=1,type="dens",
                        path_gates=NULL){
  nodes_data<-visnetdata$nodes
  nodes_data$group<-as.numeric(nodes_data$group)
  groups<-unique(nodes_data$group)
  # get path expression files
  print("----- get path expression files")
  path_exprs_data_files<-list.files(path_expr_data,full.names = T,pattern = "*.csv",recursive = T)
  path_exprs_data_files_vec<-sapply(1:length(path_exprs_data_files),function(i){
    current_path<-path_exprs_data_files[i]
    splitted_path<-strsplit(current_path,"/")[[1]]
    name_path<-tail(splitted_path,1)
    return(name_path)
  })
  # path_exprs_data_files_vec contains only the names of the files without the path
  if(is.null(path_gates)==F){
    print("----- get gates data path")
    gates_files_paths<-list.files(path_gates,full.names = T,pattern = "*.csv",recursive = T)
  }
  # export plots
  print("----- export plots")
  info_exporting<-mclapply(1:length(groups),function(i){
    current_group<-groups[i]
    print(sprintf("current group exporting: %s",current_group))
    if(is.null(path_output)==F){
      path_output_group<-paste0(path_output,sprintf("/group_%s",current_group))
      suppressWarnings(dir.create(path_output_group,recursive = T))
    }
    inds<-which(nodes_data$group==current_group)
    samples_names_group_i<-nodes_data$id[inds]
    print("samples names current group:")
    print(samples_names_group_i)
    if(is.null(path_output)==F){
      for(sample_name in samples_names_group_i){
        # get writing path current file
        path_output_file<-paste0(path_output_group,sprintf("/%s.png",sample_name))
        # get expression data
        ind<-which(path_exprs_data_files_vec==sample_name)
        if(length(ind)!=1){
          stop("0 or multiple matches for sample name to export")
        }
        df_exprs<-read.csv(path_exprs_data_files[ind])
        if(is.null(path_gates)==F){
          ind<-grep(sample_name,gates_files_paths,fixed = T)
          if(length(ind)==0){
            # open file
            png(file = path_output_file,  width = 800, height = 800)
            flowSim_plot(df_exprs) # generate image
            # close file
            dev.off()
          }else{
            df_gate<-read.csv(gates_files_paths[ind],check.names = F)
            df_exprs_final<-cbind(df_exprs,df_gate[,1])
            # open file
            png(file = path_output_file,  width = 800, height = 800)
            flowSim_plot(df_exprs_final,plot_gate=T) # generate image with gates
            # close file
            dev.off()
          }
          
        }else{
          # open file
          png(file = path_output_file,  width = 800, height = 800)
          flowSim_plot(df_exprs) # generate image
          # close file
          dev.off()
        }
      }
    }
    
    return(samples_names_group_i)
  },mc.cores = n_cores)
  names(info_exporting)<-groups
  return(info_exporting)
}

#' exports_filtered_plots
#' 
#' function to export only the files to keep for each group
#' @param visnetdata list of dataframes containing the edges and nodes information in VisNetwork format
#' @param df_features A dataframe containing the original features set (before the tsne reduction) of all input files
#' @param path_expr_data Path to directory containing the expression data of the files analyzed
#' @param path_output Path to the directory where to export the plots for each group
#' If path_output=NULL, no plots are generated. Default to NULL
#' @param n_cores Number of cores to use. Default to 1
#' @param filter_thr float number indicating the cut height for the dendrogram generated 
#' by the hierarchical clustering. Default to 2.1
#' @param plot_dendrogram If True, the plot of the dendrogram is generated (It may fail for some groups). Default to False.
#' @param seed_n Set seed number. Default to 40.
#' @return A list of the selected files names to keep for each group
#' @export
#' @examples 
#' \donttest{exports_filtered_plots(visnetdata=visnetdata,df_features=df_features,
#' path_expr_data="path to input directory",path_output=NULL,n_cores=1,filter_thr=2.1,plot_dendrogram=F)}

exports_filtered_plots<-function(visnetdata,df_features,path_expr_data,
                                 path_output=NULL,n_cores=1,type="dens",filter_thr=2.1,
                                 plot_dendrogram=F,seed_n=40){
  set.seed(seed_n)
  nodes_data<-visnetdata$nodes
  nodes_data$group<-as.numeric(nodes_data$group)
  groups<-unique(nodes_data$group)

  ################ export filtered plots of subgroups 
  print("------ export filtered plots of subgroups -------")
  info_exporting<-mclapply(1:length(groups),function(i){
    current_group<-groups[i]
    #--------- get samples names current group
    print(sprintf("current group exporting: %s",current_group))
    if(is.null(path_output)==F){
      path_output_group<-paste0(path_output,sprintf("/group_%s",current_group))
      suppressWarnings(dir.create(path_output_group,recursive = T))
    }
    inds<-which(nodes_data$group==current_group)
    samples_names_group_i<-nodes_data$id[inds]
    #------- analyze features samples current group ---------
    print("analyze features samples current group")
    if(length(samples_names_group_i)==1){
      samples_names_selected<-samples_names_group_i
    }else{
      # ------- get features of the samples of the current group
      check_vec<-row.names(df_features) %in% samples_names_group_i
      inds<-which(check_vec==T)
      df_features_group<-df_features[inds,]
      if(is.null(path_output)==F){
        path_file<-paste0(path_output_group,sprintf("/df_features_%s.csv",current_group))
        write.csv(df_features_group,file = path_file,row.names = T)
      }  
      #---- scale features ---
      inds_to_scale<-c(1,2,5,6,7,8)
      df_features_group_scaled_temp<-abs(scale(df_features_group[,inds_to_scale]))
      df_features_group_scaled<-cbind(df_features_group_scaled_temp,df_features_group[,-inds_to_scale])
      if(is.null(path_output)==F){
        path_file<-paste0(path_output_group,sprintf("/df_features_scaled_%s.csv",current_group))
        write.csv(df_features_group_scaled,file = path_file,row.names = T)
      }
      #----------- get subgroups (clustering)
      dist_obj<-dist(df_features_group_scaled,method = "euclidean")
      df_dist<-as.data.frame(as.matrix(dist_obj),stringsAsFactors=F)
      row.names(df_dist)<-row.names(df_features_group_scaled)
      out_hclustering<-hclust(dist_obj,method = "complete")
      clusters<-cutree(tree = out_hclustering,h = filter_thr)
      if(is.null(path_output)==F){
        path_file<-paste0(path_output_group,sprintf("/df_dist_%s.csv",current_group))
        write.csv(df_dist,file = path_file,row.names = T)
      }
      #--------- plot clustering
      if(plot_dendrogram==T){
        path_output_file<-paste0(path_output_group,sprintf("/clustering_%s.png",current_group))
        # open file
        png(file = path_output_file,  width = 800, height = 800)
        plot(out_hclustering)
        # close file
        dev.off()
      }
      # ------- get one sample from each subgroup 
      df_features_sub_groups<-as.data.frame(cbind(row.names(df_features_group),clusters),stringsAsFactors=F)
      colnames(df_features_sub_groups)<-c("samples_name","clusters")
      if(is.null(path_output)==F){
        path_file<-paste0(path_output_group,sprintf("/clustering_%s.csv",current_group))
        write.csv(df_features_sub_groups,file = path_file,row.names = F)
      }
      unique_subgroups<-unique(df_features_sub_groups$clusters)
      samples_names_selected<-rep("a",length(unique_subgroups))
      for(g in 1:length(unique_subgroups)){
        current_subgroup<-unique_subgroups[g]
        inds<-which(df_features_sub_groups$clusters==current_subgroup)
        samples_subgroup<-df_features_sub_groups$samples_name[inds]
        set.seed(seed_n)
        rand_int<-sample(1:length(samples_subgroup),1)
        sample_select_current_subgroup<-samples_subgroup[rand_int]
        samples_names_selected[g]<-sample_select_current_subgroup
      }
    }
    #---------- export samples selected---------------
    if(is.null(path_output)==F){
      ################ get path expression files ----------
      print("----- get path expression files -----")
      path_exprs_data_files<-list.files(path_expr_data,full.names = T,pattern = "*.csv",recursive = T)
      path_exprs_data_files_vec<-sapply(1:length(path_exprs_data_files),function(i){
        current_path<-path_exprs_data_files[i]
        splitted_path<-strsplit(current_path,"/")[[1]]
        name_path<-tail(splitted_path,1)
        return(name_path)
      })
      print("export samples selected")
      #print(samples_names_selected)
      for(sample_name in samples_names_selected){
        # get writing path current file
        path_output_file<-paste0(path_output_group,sprintf("/%s.png",sample_name))
        # get expression data
        ind<-which(path_exprs_data_files_vec==sample_name)
        if(length(ind)!=1){
          stop("0 or multiple matches for sample name to export")
        }
        df_exprs<-read.csv(path_exprs_data_files[ind])
        # open file
        png(file = path_output_file,  width = 800, height = 800)
        flowSim_plot(df_exprs) # generate image
        # close file
        dev.off()
      }
    }
    samples_names_selected<-paste0(samples_names_selected,collapse = ";")
    return(c(current_group,samples_names_selected))
  },mc.cores = n_cores)
  info_exporting<-do.call(rbind,info_exporting)
  info_exporting_df<-as.data.frame(info_exporting,stringsAsFactors=F)
  colnames(info_exporting_df)<-c("group","files_filtered")
  vec<-as.character(info_exporting_df$files_filtered)
  vec_filtered<-unlist(strsplit(vec,split = ";"))
  print(sprintf("Number of files selected: %d",length(vec_filtered)))
  return(info_exporting_df)
}
