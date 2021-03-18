##===============================================
##Some functions to analyse data
##===============================================

##===============================================
##To get the descendant edge of an internal ndoe
##===============================================

getDescendantEdges<-function(tree,node,curr=NULL, edge_list=NULL){
  edge_list <- c()
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  edge_list <- c(edge_list, which(tree$edge[,1]==node))
  curr<-c(curr,daughters)
  w<-which(daughters>length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(edge_list)
}

##===============================================
##Converting trees to nhx like labeled trees, required by bppml
##===============================================

convert_to_nhx <- function(in_tree_str){
  tree_str <- in_tree_str
  tree_vec <- strsplit(tree_str, split = "")[[1]]
  tree_vec_nhx <- c()
  
  node_id <- 0
  for(char_n in 1:length(tree_vec)){
    char_oi <- tree_vec[char_n]
    if(char_oi==":"){
      char2_n <- length(tree_vec_nhx)
      while(!(tree_vec_nhx[char2_n] %in% c(",", "(", ")"))){
        char2_n <- char2_n -1
      }
      tree_vec_nhx <- c(tree_vec_nhx[1:char2_n], c(node_id,char_oi))
      node_id <- node_id+1
    }else{
      tree_vec_nhx <- c(tree_vec_nhx, char_oi)
    }
  }
  tree_str_nhx <- paste0(tree_vec_nhx, collapse="")
  
  return(tree_str_nhx)
}

##===============================================
## A mapping of node labels required to shuffle between phylo tree objects and NHX trees
##bppml asks nodes to be numbered in the serial order in which they appear in newick string (achieved explicitly in NHX)
##but phylo has its own numbering scheme with all tips labelled first
##===============================================
newick_to_nhx_map <- function(in_tree_str){
  tree_str <- in_tree_str
  tree_vec <- strsplit(tree_str, split = "")[[1]]
  tree_vec_nhx <- c()
  
  mapped_strings <- lapply(1:i)
  
  node_id <- 0
  for(char_n in 1:length(tree_vec)){
    char_oi <- tree_vec[char_n]
    if(char_oi==":"){ #when you hit a ":"
      char2_n <- length(tree_vec_nhx)
      while(!(tree_vec_nhx[char2_n] %in% c(",", "(", ")"))){ ##as long as you dont find a separator
        char2_n <- char2_n -1 ##keep going back
      }
      tree_vec_nhx <- c(tree_vec_nhx[1:char2_n], c(node_id,char_oi)) ##add a node id and put back the ":"
      node_id <- node_id+1
      
      ##map the node id to the replaced string
      
    }else{ ##keep adding the characters unless you hit a ":"
      tree_vec_nhx <- c(tree_vec_nhx, char_oi)
    }
  }
  tree_str_nhx <- paste0(tree_vec_nhx, collapse="")
  
  return(tree_str_nhx)
}

##===============================================
##Actual analysis code
##===============================================
library(phytools)
library(cowplot)
library(tidyverse)
library(seqinr)

main_dir <- "/"
mapnh_dir <- paste0(main_dir, "Results/mapNH/")
dNdS_dir <- paste0(main_dir, "Results/dNdS/")
plots_dir <- paste0(main_dir, "Results/focal_trees/")

focal_clade_details_fn <- paste0(dNdS_dir, "focal_clades_details_nonHEG_v2.tab") ##a file containing dataset names and which model sets (parameter regimes) correspond to focal clades and whcih to control clades
focal_clade_details <- read.table(file = focal_clade_details_fn, header = T, row.names = 1, stringsAsFactors = F)
ds_vec <- rownames(focal_clade_details)
ds_parent <- focal_clade_details$parent_clade
  
set_n <- 1
set_size <- 40
gene_sets <- c("nonHEG", "HEG")

##data to be stored in a dataframe with following headings (eq=equilibrium, act=actual)
col_names_df <- c("focal_clade", "gene_set", "eqGC3_root", "eqGC3_focal", "eqGC3_sister", "actGC3_focal", "actGC3_sister", "omega_focal", "omega_sister", "dS_focal", "dS_sister", "dN_focal", "dN_sister", "dS_focal_max", "dS_other_max", "dS_max")
dNdS_gc3_data <- as.data.frame(matrix(data = NA, nrow = 2*nrow(focal_clade_details), ncol = length(col_names_df)))
colnames(dNdS_gc3_data) <- col_names_df

for(ds_n in 1:nrow(focal_clade_details)){
  ds_oi <- ds_vec[ds_n]
  ds_oi_dir <- ds_parent[ds_n]
  print(ds_oi)
  
  mapnh_sp_dir <- paste0(mapnh_dir,ds_oi_dir,"/")
  mapnh_ds_dir <- paste0(mapnh_sp_dir,ds_oi_dir,"/")
  
  #convert original tree to a node numbered tree
  newick_tree_fn <- list.files(path = mapnh_sp_dir, pattern = ".newick", full.names = F)[1]
  newick_str <- scan(file = paste0(mapnh_sp_dir,newick_tree_fn), what = "character")
  newick_tree <- read.tree(paste0(mapnh_sp_dir,newick_tree_fn))
  
  nhx_str <- convert_to_nhx(newick_str)
  writeLines(text = nhx_str, con = paste0(mapnh_sp_dir,gsub(pattern = ".newick", replacement = ".nodeName", x = newick_tree_fn)))
  start_tree_fn <- list.files(path = mapnh_sp_dir, pattern = ".nodeName", full.names = F)[1]
  in_tree <- read.tree(paste0(mapnh_sp_dir, start_tree_fn))

  ##To get the number of models and list nodes belonging to each model
  ##lets read shared details of the model from nonHEG params file
  model_params_fn <- list.files(path = mapnh_ds_dir, pattern = glob2rx(paste0("*nonHEG_orthologs",set_size,"_set",set_n,".params")), full.names = F)
  if(identical(model_params_fn, character(0))) next ##skip if no parameter file is found
  params_oi <- readLines(con = paste0(mapnh_ds_dir,model_params_fn))
  n_models_line <- params_oi[which(grepl("nonhomogeneous.number_of_models",params_oi))]
  n_models <- as.integer(strsplit(n_models_line, split = "=")[[1]][2])
  
  model_nodes <- lapply(1:n_models, function(x) NA)
  for(model_n in 1:n_models){
    model_nodes_line <- params_oi[which(grepl(paste0("model",model_n,".nodes_id="),params_oi))]
    model_nodes[[model_n]] <- strsplit(strsplit(model_nodes_line,"=")[[1]][2], ",")[[1]]
  }
  focal_model <- focal_clade_details$focal_model[ds_n]
  control_model <- focal_clade_details$control_model[ds_n]
  all_focal_model_nodes <- unlist(sapply(setdiff(1:n_models,control_model), function(x) model_nodes[[x]]))
  
  for(gene_set_n in 1:length(gene_sets)){
    gene_set <- gene_sets[gene_set_n]
    
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"focal_clade"] <- ds_vec[ds_n]
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"gene_set"] <- gene_set

    ##read model information for gene set of interest
    model_params_fn <- list.files(path = mapnh_ds_dir, pattern = glob2rx(paste0("*_",gene_set,"_orthologs",set_size,"_set", set_n,".params*")), full.names = F)
    if(identical(model_params_fn,character(0))){
      next
    }else if(length(model_params_fn)>1){ #if there more than one file search for exact file path required
        model_params_fn <- paste0("model_nonhom_ml_", ds_parent[ds_n], "_", gene_set, "_orthologs", set_size, "_set", set_n,".params")
      }
    params_oi <- readLines(con = paste0(mapnh_ds_dir,model_params_fn))

    root_params_line <- params_oi[which(grepl("nonhomogeneous.root_freq=",params_oi))]
    eqGC3_param_pos <- gregexpr(pattern = "3_Full.theta=", text = root_params_line)
    eqGC3_root_oi <- paste0(strsplit(root_params_line,"")[[1]][(eqGC3_param_pos[[1]][1]+13):(eqGC3_param_pos[[1]][1]+16)], collapse = "")
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n, "eqGC3_root"] <- as.numeric(eqGC3_root_oi)
    
    model_eqGC3 <- lapply(1:n_models, function(x) NA)
    model_omega <- lapply(1:n_models, function(x) NA)
    for(model_n in 1:n_models){
      model_params_line <- params_oi[which(grepl(paste0("model",model_n,"="),params_oi))]
      eqGC3_param_pos <- gregexpr(pattern = "3_Full.theta=", text = model_params_line)
      model_eqGC3[[model_n]] <- paste0(strsplit(model_params_line,"")[[1]][(eqGC3_param_pos[[1]][1]+13):(eqGC3_param_pos[[1]][1]+16)], collapse = "")
      
      omega_param_pos <- gregexpr(pattern = "omega=", text = model_params_line)
      model_omega[[model_n]] <- strsplit(strsplit(model_params_line,"=")[[1]][14], ")")[[1]]
    }
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n, "eqGC3_focal"] <- as.numeric(model_eqGC3[[focal_model]])
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"eqGC3_sister"] <- as.numeric(model_eqGC3[[control_model]])
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"omega_focal"] <- as.numeric(model_omega[[focal_model]])
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"omega_sister"] <- as.numeric(model_omega[[control_model]])
    
    ##identify sister tips for comparison of branch lengths (substitution counts)
    ##**this works differently if focal clade has single or multiple taxa**
    focal_nodes <- model_nodes[[focal_model]] #bppml numbering
    sister_nodes <- model_nodes[[control_model]] #bppml numbering
    
    # if(length(focal_nodes)>1 & length(sister_nodes)>1){
      
    # }
    ##if focal set has single tip and no internal nodes, then its ancestor is also MRCA of focal and sister clades
    if(length(focal_nodes)==1){
      print("case 1, focal clade has single node")
      focal_tips_phylo <- which(in_tree$tip.label %in% focal_nodes) #focal node number in phylo object
      focal_tips_labels <- in_tree$tip.label[focal_tips_phylo] #bppml numbering
      focal_mrca <- getParent(tree = in_tree, node = focal_tips_phylo) #in phylo numbering
      focal_sister_anc <- focal_mrca #in phylo numbering
      if(length(sister_nodes)==1){
        sister_mrca <- focal_sister_anc
      }else{
        sister_mrca <- setdiff(in_tree$edge[in_tree$edge[,1]==focal_sister_anc,2], focal_tips_phylo) #in phylo numbering
      }
      sister_tips_phylo <- intersect(1:length(in_tree$tip.label), getDescendants(tree = in_tree, node = sister_mrca))
      sister_tips_phylo <- setdiff(sister_tips_phylo, which(in_tree$tip.label %in% all_focal_model_nodes))
      sister_tips_labels <- in_tree$tip.label[sister_tips_phylo]
      #if all sister tips are also in another focal set, look for outgroup as control
      if(identical(setdiff(sister_tips_labels, all_focal_model_nodes),character(0))){
        print("case 1b, sister clade also has a jump on its stem branch; looking for outgroup as control")
        focal_sister_anc2 <- getParent(tree = in_tree, node = focal_sister_anc)
        if(length(sister_nodes)==1){
          sister_mrca <- focal_sister_anc2
        }else{
          sister_mrca <- setdiff(in_tree$edge[in_tree$edge[,1]==focal_sister_anc2,2], focal_sister_anc) #in phylo numbering
        }
        sister_tips_phylo <- intersect(1:length(in_tree$tip.label), getDescendants(tree = in_tree, node = sister_mrca))
        sister_tips_phylo <- setdiff(sister_tips_phylo, which(in_tree$tip.label %in% all_focal_model_nodes))
        sister_tips_labels <- in_tree$tip.label[sister_tips_phylo]
        anc_for_branch_lengths <- focal_sister_anc2
      }else{
        anc_for_branch_lengths <- focal_sister_anc
      }
    #if the focal set has single tip but also some internal nodes, implying that there is a jump on sister branch
    }else if(length(intersect(focal_nodes,in_tree$tip.label))==1){
      print("case 2, sister clade has a nested GC jump; remove the taxa affected by the nested jump from control group")
      focal_tips_phylo <- which(in_tree$tip.label %in% focal_nodes)
      focal_tips_labels <- in_tree$tip.label[focal_tips_phylo] #bppml numbering
      focal_mrca <- getParent(tree = in_tree, node = focal_tips_phylo) #phylo numbering
      focal_sister_anc <- getParent(tree = in_tree, node = focal_mrca) #phylo numbering
      if(length(sister_nodes)==1){
        sister_mrca <- focal_sister_anc
      }else{
        sister_mrca <- setdiff(in_tree$edge[in_tree$edge[,1]==focal_sister_anc,2], focal_mrca) #in phylo numbering
      }
      sister_tips_phylo <- intersect(1:length(in_tree$tip.label), getDescendants(tree = in_tree, node = sister_mrca)) #phylo numbering
      sister_tips_phylo <- setdiff(sister_tips_phylo, which(in_tree$tip.label %in% all_focal_model_nodes))
      sister_tips_labels <- in_tree$tip.label[sister_tips_phylo]
      anc_for_branch_lengths <- focal_sister_anc
    #if the focal set has multiple tips
    }else{
      print("case 3, normal sister clade with multiple taxa")
      focal_tips_phylo <- which(in_tree$tip.label %in% focal_nodes)
      focal_tips_labels <- in_tree$tip.label[focal_tips_phylo] #bppml numbering
      focal_mrca <- getMRCA(phy = in_tree, tip = focal_tips_phylo) #phylo numbering
      focal_sister_anc <- in_tree$edge[which(in_tree$edge[,2]==focal_mrca),1] #phylo numbering
      if(length(sister_nodes)==1){
        sister_mrca <- focal_sister_anc
      }else{
        sister_mrca <- setdiff(in_tree$edge[in_tree$edge[,1]==focal_sister_anc,2], focal_mrca) #in phylo numbering
      }
      sister_tips_phylo <- intersect(1:length(in_tree$tip.label), getDescendants(tree = in_tree, node = sister_mrca))
      sister_tips_phylo <- setdiff(sister_tips_phylo, which(in_tree$tip.label %in% all_focal_model_nodes))
      sister_tips_labels <- in_tree$tip.label[sister_tips_phylo]
      ##if all sister tips are also in another focal clade, look for outgroup as sister set
      if(identical(setdiff(sister_tips_labels, all_focal_model_nodes),character(0))){
        print("case 3b, sister clade also affected by jump; look for outgroup as control")
        focal_sister_anc2 <- getParent(tree = in_tree, node = focal_sister_anc)
        sister_mrca <- setdiff(in_tree$edge[which(in_tree$edge[,1]==focal_sister_anc2),2], focal_sister_anc)
        sister_tips_phylo <- intersect(1:length(in_tree$tip.label), getDescendants(tree = in_tree, node = sister_mrca))
        sister_tips_phylo <- setdiff(sister_tips_phylo, which(in_tree$tip.label %in% all_focal_model_nodes))
        sister_tips_labels <- in_tree$tip.label[sister_tips_phylo]
        anc_for_branch_lengths <- focal_sister_anc2
      }else{
        anc_for_branch_lengths <- focal_sister_anc
      }
    }
    
    short_label_tree <- newick_tree
    # GCF_short_names <- read.table(paste0(main_dir, "Input_data/phylogenies/", ds_oi_dir, "/", ds_oi_dir,"_GCF_ids_shortnames.txt"), stringsAsFactors = F, row.names = 1, header = F, sep="\t")
    GCF_short_names <- read.table(paste0("/home/saurabh/Documents/Deepas_lab_Backup/Deepas_lab/Work/GC_macro/Common_data/", ds_oi_dir, "/", ds_oi_dir,"_GCF_ids_shortnames.txt"), stringsAsFactors = F, row.names = 1, header = F, sep="\t")
    short_label_tree$tip.label <- GCF_short_names[short_label_tree$tip.label,"V2"]
    if(anc_for_branch_lengths==focal_mrca){
      focal_edges <- which(short_label_tree$edge[,1]==anc_for_branch_lengths & short_label_tree$edge[,2]==focal_tips_phylo)
    } else{
      focal_edges <- unique(c(which(short_label_tree$edge[,1]==anc_for_branch_lengths & short_label_tree$edge[,2]==focal_mrca), which(short_label_tree$edge[,2] %in% getDescendants(short_label_tree, focal_mrca)))) 
    }
    control_edges <- unique(c(which(short_label_tree$edge[,1]==anc_for_branch_lengths & short_label_tree$edge[,2]==sister_mrca), which(short_label_tree$edge[,2] %in% setdiff(getDescendants(short_label_tree, sister_mrca), c(which(in_tree$tip.label %in% all_focal_model_nodes), length(short_label_tree$tip.label) + which(in_tree$node.label %in% all_focal_model_nodes))))))
    col_edges <- rep("black", nrow(short_label_tree$edge))
    col_edges[focal_edges] <- "orange"
    col_edges[control_edges] <- "purple"
    svg(filename = paste0(plots_dir, ds_oi, "_focal_clade_tree.svg"), width = 5, height = 3)
    plot(x = short_label_tree, main = paste0(ds_oi, "\n(orange=focal; purple=control)"), edge.color = col_edges, edge.width = 2)
    dev.off()
    
    ##get actual GC content of the focal clade and sister clade sequences
    seq_fn_pattern <- paste0("*_", gene_set,"_orthologs",set_size,"_set", set_n, ".fa")
    if(identical(list.files(path = mapnh_sp_dir, pattern = glob2rx(seq_fn_pattern), full.names = F),character(0))){
      print("sequence file not found!")
    }else{
      seq_fn <- list.files(path = mapnh_sp_dir, pattern = glob2rx(seq_fn_pattern), full.names = T)
      seqs_oi <- read.fasta(file = seq_fn, seqtype = "DNA", as.string = F)
      GC3_seqs <- sapply(seqs_oi, function(x) GC3(x))
    }
      
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n, "actGC3_focal"] <- median(GC3_seqs[newick_tree$tip.label[focal_tips_phylo]])
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"actGC3_sister"] <- median(GC3_seqs[newick_tree$tip.label[sister_tips_phylo]])
    
    ##get dN, dS tree files
    dN_fn_pattern <- paste0("*_", gene_set,"_orthologs",set_size,"_set", set_n,".fa.counts_dN.dnd*")
    dS_fn_pattern <- paste0("*_", gene_set,"_orthologs",set_size,"_set", set_n,".fa.counts_dS.dnd*")
    if(identical(list.files(path = mapnh_sp_dir, pattern = glob2rx(dN_fn_pattern), full.names = F),character(0))){
      print("substitution rate estimate tree files not found")
    }else{
      dN_tree_fn <- list.files(path = mapnh_sp_dir, pattern =  glob2rx(dN_fn_pattern), full.names = F)
      dS_tree_fn <- list.files(path = mapnh_sp_dir, pattern =  glob2rx(dS_fn_pattern), full.names = F)
    }
    
    dN_tree_nodeNamed <- read.tree(text = convert_to_nhx(scan(file = paste0(mapnh_sp_dir,dN_tree_fn), what = "character")))
    dS_tree_nodeNamed <- read.tree(text = convert_to_nhx(scan(file = paste0(mapnh_sp_dir,dS_tree_fn), what = "character")))
    
    ##get maximum dS for finding saturation
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"dS_focal_max"] <- max(dS_tree_nodeNamed$edge.length[focal_edges])
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"dS_other_max"] <- max(dS_tree_nodeNamed$edge.length[dS_tree_nodeNamed$edge.length[-focal_edges]])
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"dS_max"] <- max(dS_tree_nodeNamed$edge.length)
      
    ##get cumulative branch lengths for focal and sister tips
    focal_sister_anc_dS_height <- nodeheight(tree = dS_tree_nodeNamed, node = anc_for_branch_lengths)
    focal_tips_dS_height <- nodeHeights(tree = dS_tree_nodeNamed)[which(dS_tree_nodeNamed$edge[,2] %in% focal_tips_phylo),2]
    sister_tips_dS_height <- nodeHeights(tree = dS_tree_nodeNamed)[which(dS_tree_nodeNamed$edge[,2] %in% sister_tips_phylo),2]
    focal_tips_median_dS_length <- median(focal_tips_dS_height - focal_sister_anc_dS_height)
    sister_tips_median_dS_length <- median(sister_tips_dS_height - focal_sister_anc_dS_height)
    
    focal_sister_anc_dN_height <- nodeheight(tree = dN_tree_nodeNamed, node = anc_for_branch_lengths)
    focal_tips_dN_height <- nodeHeights(tree = dN_tree_nodeNamed)[which(dN_tree_nodeNamed$edge[,2] %in% focal_tips_phylo),2]
    sister_tips_dN_height <- nodeHeights(tree = dN_tree_nodeNamed)[which(dN_tree_nodeNamed$edge[,2] %in% sister_tips_phylo),2]
    focal_tips_median_dN_length <- median(focal_tips_dN_height - focal_sister_anc_dN_height)
    sister_tips_median_dN_length <- median(sister_tips_dN_height - focal_sister_anc_dN_height)
    
    ##MOST IMPORTANT: take the ratio of median cumulative lengths of focal vs control (sister) clade
    median_dS_ratio <- focal_tips_median_dS_length/sister_tips_median_dS_length
    median_dN_ratio <- focal_tips_median_dN_length/sister_tips_median_dN_length
    
    ##Save substitution rate data to the dataframe
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"dS_focal"] <- focal_tips_median_dS_length
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"dS_sister"] <- sister_tips_median_dS_length
    
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"dN_focal"] <- focal_tips_median_dN_length
    dNdS_gc3_data[(ds_n-1)*2+gene_set_n,"dN_sister"] <- sister_tips_median_dN_length
  }
}

dNdS_gc3_data$eqDirectionA <- (dNdS_gc3_data$eqGC3_focal - dNdS_gc3_data$eqGC3_root)/abs(dNdS_gc3_data$eqGC3_focal - dNdS_gc3_data$eqGC3_root)
dNdS_gc3_data$eqDirectionB <- (dNdS_gc3_data$eqGC3_focal - dNdS_gc3_data$eqGC3_sister)/abs(dNdS_gc3_data$eqGC3_focal - dNdS_gc3_data$eqGC3_sister)
dNdS_gc3_data$eqDirectionC <- (dNdS_gc3_data$actGC3_focal - dNdS_gc3_data$actGC3_sister)/abs(dNdS_gc3_data$actGC3_focal - dNdS_gc3_data$actGC3_sister)
write.table(x = dNdS_gc3_data, file = paste0(dNdS_dir,"focal_clades_dNdS_gc3_metadata_orthologs", set_size,"_set",set_n,"_revised.tab"), quote = F, row.names = F, col.names = T)
