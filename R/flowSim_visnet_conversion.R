
#' gen_visnetwork_data 
#' 
#' function to convert igraph network to a list of two VisNetwork dataframes containig 
#' the edges and nodes information needed to generate a VisNetwork object
#' @param igraph_network Igraph network object
#' @param contract_net If True, it generates a contracted version of the 
#' VisNetwork dataframes (to visualize large datasets). Default to False
#' @param single_nodes_black Isolated files/nodes colored in black? Default to True.
#' @return If contracted_net=False, a list of full VisNetwork dataframes is generated, if 
#' contracted_net=True, a list of contracted VisNetwork dataframes is returned.
#' $nodes : dataframe containing the nodes information (full or contracted)
#' $edges : dataframe containing the edges information (full or contracted)
#' @export
#' @examples 
#' \donttest{gen_visnetwork_data(igraph_network=igraph_network,contract_net=F)}
gen_visnetwork_data<-function(igraph_network,contract_net=F,single_nodes_black=T){
  data<-toVisNetworkData(igraph_network)
  data$edges$weight<-as.numeric(data$edges$weight)
  data$nodes$group<-as.character(data$nodes$group)
  #------------------------- add edges tooltip ----------------------------
  print("-------generating edges tooltip-----")
  title_column<-as.character(data$edges$weight)
  data$edges$title<-factor(title_column)
  #------------------------- add width column ----------------------------
  print("-------generating edges width-----")
  width_column<-data$edges$weight*2
  data$edges$width<-width_column

  #---------- modify label column ---------------
  print("----- modify label column --------")
  new_label_column<-paste(data$nodes$group,data$nodes$label,sep="--")
  data$nodes$label<-new_label_column
  #----------------- add nodes tooltip -------------
  print("-------add nodes tooltip-----")
  data$nodes$title<-data$nodes$label
  #--------------------- hide label nodes --------------------------
  print("-------hide nodes labels-----")
  x<-rep(NA,nrow(data$nodes))
  data$nodes$label<-x
  #--------------------- set size nodes --------------------------
  print("-------set size nodes-----")
  x<-rep(10,nrow(data$nodes))
  data$nodes$size<-x
  #--------- add color for clusters composed by single nodes -----------
  print("---- add color column ----")
  df_groups<-data$nodes[,c("id","group")]
  df_colors<-get_new_colors_based_on_groups(df_groups)
  data$nodes$color<-df_colors[,3]
  data$nodes$color<-as.character(data$nodes$color)
  #--------- add size of clusters (number of nodes in each clusters) ------
  vec_sizes<-get_size_based_on_groups(df_groups)
  data$nodes$n_files<-vec_sizes
  data$nodes$n_files<-as.integer(data$nodes$n_files)
  #----- set colors single file clusters as black ------------
  if(single_nodes_black==T){
    inds<-which(data$nodes$n_files==1)
    data$nodes$color[inds]<-"black"
  }
  #---------------------- contract network if needed
  if(contract_net==T){
    print("----- contract network ------")
    data_contracted<-contract_data(data = data)
    return(data_contracted)
  }else{
    return(data)
  }
}


#' contract_data
#' 
#' function to contract the visnetdata object (containing nodes and edges info). Internal function.
#' @param data List of the full VisNetwork dataframe containing the edges and nodes information
#' @keywords Internal
#' @return A list of two VisNetwork dataframes containing the contracted information of nodes and edges
contract_data<-function(data){
  data_contracted<-list()
  data_edges<-data$edges
  data_nodes<-data$nodes
  #----------------------- make nodes data -----------------------
  print("######################## make nodes data #######################")
  unique_groups<-unique(data_nodes$group)
  vec_ids<-rep("a",length(unique_groups))
  vec_color<-rep("blue",length(unique_groups))
  vec_size<-rep(8,length(unique_groups))
  vec_label<-rep(NA,length(unique_groups))
  vec_title<-rep("a",length(unique_groups))
  vec_samples_ids<-rep("a",length(unique_groups))
  for(i in 1:length(unique_groups)){
    current_group<-unique_groups[i]
    inds<-which(data_nodes$group==current_group)
    samples_ids_current_group<-data_nodes$id[inds]
    n_files_current_group<-unique(data_nodes$n_files[inds])
    vec_ids[i]<-current_group
    vec_size[i]<-as.integer(n_files_current_group)
    string<-sprintf("group: %s, %d files",current_group,n_files_current_group)
    vec_title[i]<-string
    string<-paste0(samples_ids_current_group,collapse = "#")
    vec_samples_ids[i]<-string
  }
  data_nodes_contracted<-data.frame(id=vec_ids,color=vec_color,size=vec_size,
                                    label=vec_label,title=vec_title,
                                    files_id=vec_samples_ids,stringsAsFactors = F)
  # remove single files clusters
  inds<-which(data_nodes_contracted$size==1)
  data_nodes_contracted<-data_nodes_contracted[-inds,]
  data_contracted[["nodes"]]<-data_nodes_contracted
  
  #------------------------- make edges data ------------------------
  print("#################### make edges data #######################")
  list_info_connections<-list()
  print("--------- analysis connections all groups --------")
  for(i in 1:nrow(data_nodes_contracted)){
    current_group_id<-data_nodes_contracted$id[i]
    #print(sprintf("------- analysis connections group: %s ----------",current_group_id))
    files_id_current_group<-data_nodes_contracted$files_id[i]
    files_id_current_group<-strsplit(files_id_current_group,"#")[[1]]
    #print("files id current group:")
    #print(files_id_current_group)
    #print("----- check connections all files id current group")
    vec_groups_connected<-c()
    for(file_id in files_id_current_group){
      #print(sprintf("check connections from---to: %s",file_id))
      inds<-grep(file_id,data_edges$from,fixed=T)
      data_edges_current_file<-data_edges[inds,]
      files_id_connected<-data_edges_current_file$to
      check_id<-data$nodes$id %in% files_id_connected
      inds<-which(check_id==T)
      group_connected_to<-data$nodes$group[inds]
      #print("connected to group:")
      vec_groups_connected<-c(vec_groups_connected,group_connected_to)
    }
    # remove internal connections info
    check_internal_con<-vec_groups_connected %in% current_group_id
    inds<-which(check_internal_con==T)
    vec_groups_connected<-vec_groups_connected[-inds]
    if(length(vec_groups_connected)==0){
      vec_groups_connected<-"None"
    }
    list_info_connections[[current_group_id]]<-vec_groups_connected
  }
  check_None<-sapply(list_info_connections,function(x){
    result_check<-x=="None"
    return(unique(result_check))
    })
  inds<-which(check_None==T)
  list_info_connections<-list_info_connections[-inds]
  print(list_info_connections)
  if(length(list_info_connections)==0){
    print("All clusters are unconnected")
    data_edges_contracted<-data.frame(from=character(),to=character(),weight=integer())
    data_contracted[["edges"]]<-data_edges_contracted
    return(data_contracted)
  }
  print("----- convert list connections ------")
  all_groups_connected<-names(list_info_connections)
  list_df_edges_all_groups<-list()
  for(i in 1:length(all_groups_connected)){
    current_group_name<-all_groups_connected[i]
    info_connections_current_group<-list_info_connections[[current_group_name]]

    unique_connected_group<-unique(info_connections_current_group)
    n_groups_connected<-length(unique_connected_group)
    vec_from_current_group<-rep(current_group_name,n_groups_connected)
    vec_to_current_group<-rep("a",n_groups_connected)
    vec_weight_current_group<-rep(1,n_groups_connected)
    for(j in 1:length(unique_connected_group)){
      group<-unique_connected_group[j]
      inds<-which(info_connections_current_group==group)
      n_edges_connected_current_group<-length(inds)
      vec_to_current_group[j]<-group
      vec_weight_current_group[j]<-n_edges_connected_current_group
    }
    temp<-rep("n_connections: ",length(vec_weight_current_group))
    vec_title_current_group<-paste0(temp,vec_weight_current_group)
    m_edges<-cbind(vec_from_current_group,vec_to_current_group,vec_weight_current_group,vec_title_current_group)
    df_edges_current_group<-as.data.frame(m_edges,stringsAsFactors=F)
    colnames(df_edges_current_group)<-c("from","to","weight","title")
    list_df_edges_all_groups[[i]]<-df_edges_current_group
  }
  data_edges_contracted<-do.call(rbind,list_df_edges_all_groups)
  data_contracted[["edges"]]<-data_edges_contracted
  # -------- fix biunivocal connections ------------
  print("------------------ fix biunivocal connections -----------------")
  print(" ------ find connections to fix")
  connections_to_fix<-c()
  for(i in 1:nrow(data_edges_contracted)){
    group_i_from<-data_edges_contracted$from[i]
    group_i_to<-data_edges_contracted$to[i]
    inds<-which(data_edges_contracted$to==group_i_from)
    data_edges_contracted_temp<-data_edges_contracted[inds,]
    inds<-which(data_edges_contracted_temp$from==group_i_to)
    data_edges_contracted_temp<-data_edges_contracted_temp[inds,]
    string_1<-paste(group_i_from,group_i_to,sep="+")
    group_from_temp<-data_edges_contracted_temp$from
    group_to_temp<-data_edges_contracted_temp$to
    string_to_check<-paste(group_from_temp,group_to_temp,sep="+")
    if(nrow(data_edges_contracted_temp)!=0){
        connections_to_fix<-c(connections_to_fix,string_to_check)
    }
  }
  print("------ connections to fix:")
  print(connections_to_fix)
  if(length(connections_to_fix)==0){
    print("no connections to fix")
  }else{
    # determine couples
    connections_to_fix_inds_couple<-c()
    connections_already_analyzed<-c()
    for(i in 1:length(connections_to_fix)){
      element_i<-connections_to_fix[i]
      connection_1<-strsplit(element_i,"+",fixed=T)[[1]][1]
      connection_2<-strsplit(element_i,"+",fixed=T)[[1]][2]
      string_inverse<-paste(connection_2,connection_1,sep="+")
      if(string_inverse %in% connections_already_analyzed){
        next()
      }
      check_inverse<-connections_to_fix %in% string_inverse
      ind_inverse<-which(check_inverse==T)
      check_current<-connections_to_fix %in% element_i
      ind_current<-which(check_current==T)
      # ind_inverse<-grep(string_inverse,connections_to_fix,fixed = T)
      # ind_current<-grep(element_i,connections_to_fix,fixed = T)
      inds_couple<-c(ind_current,ind_inverse)
      if(length(inds_couple)!=2){
        print(inds_couple)
        print(connections_to_fix[inds_couple])
        stop("Error: inds_couple has more than 2 integers")
      }
      inds_couple<-paste0(inds_couple,collapse = ",")
      connections_to_fix_inds_couple<-c(connections_to_fix_inds_couple,inds_couple)
      connections_already_analyzed<-c(connections_already_analyzed,element_i,string_inverse)
    }
    connections_to_fix_inds_couple<-unique(connections_to_fix_inds_couple)
    print("------ connections to fix inds pair:")
    print(connections_to_fix_inds_couple)
    # finally fix connections
    vec_string_connections<-paste(data_edges_contracted$from,data_edges_contracted$to,sep = "+")
    vec_inds_to_remove<-c()
    for(i in 1:length(connections_to_fix_inds_couple)){
      current_couple<-connections_to_fix_inds_couple[i]
      couple_1<-strsplit(current_couple,",")[[1]][1]
      couple_2<-strsplit(current_couple,",")[[1]][2]
      couple_1<-as.integer(couple_1)
      couple_2<-as.integer(couple_2)
      connection_1<-connections_to_fix[couple_1]
      connection_2<-connections_to_fix[couple_2]
      check_connection_1<-vec_string_connections %in% connection_1
      ind_1<-which(check_connection_1==T)
      data_edges_connection_1<-data_edges_contracted[ind_1,]
      weight_connection_1<-as.integer(data_edges_connection_1$weight)
      check_connection_2<-vec_string_connections %in% connection_2
      ind_2<-which(check_connection_2==T)
      data_edges_connection_2<-data_edges_contracted[ind_2,]
      weight_connection_2<-as.integer(data_edges_connection_2$weight)
      new_weight<-weight_connection_1+weight_connection_2
      new_title<-paste0("n_connections: ",new_weight)
      data_edges_connection_1$weight<-new_weight
      data_edges_connection_1$title<-new_title
      data_edges_contracted[ind_1,]<-data_edges_connection_1
      vec_inds_to_remove<-c(vec_inds_to_remove,ind_2)
    }
    data_edges_contracted<-data_edges_contracted[-vec_inds_to_remove,]
    data_contracted[["edges"]]<-data_edges_contracted
  }
  return(data_contracted)
}



