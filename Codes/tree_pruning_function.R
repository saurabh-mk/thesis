prune_tree <- function(tree_oi, dist_cutoff){
  if(dist_cutoff==0) return(tree_oi)
  require(phytools)
  dist_cutoff <- dist_cutoff
  
  tree_pruned <- tree_oi
  
  #check for duplicate tip labels, this could be a problem in subsequent analysis
  if(length(unique(tree_pruned$tip.label))<length(tree_pruned$tip.label)){
    usr_input <- readline(prompt = "Tree has duplicate tip labels. Should we proceed? (y/n): ")
    if(usr_input=="y"){
      print("You have chosen to ignore duplicates!")
    }
    else{
      print("Do something about it!")
    }
  }
  
  pruning_done <- F
  iter_oi <- 1
  
  while(pruning_done==F){
    # print(paste0("beginning iteration number:     ", as.character(iter_oi)))
    iter_oi <- iter_oi + 1
    
    #get parameters for remaining tree
    n_tips <- length(tree_pruned$tip.label)
    #order internal nodes accoding to height, start from nodes that are closest to root
    # int_nodes <- n_tips+(1:tree_pruned$Nnode)
    # int_heights <- sapply(int_nodes, function(x) nodeheight(tree_pruned, x))
    edges_oi <- cbind(tree_pruned$edge,tree_pruned$edge.length) 
    ##next we determine the edges that are shorter than dist cutoff
    edges_oi <- subset(edges_oi,edges_oi[,3]<dist_cutoff)
    ##here we discard edges that are deep in the tree
    ##i.e only deal with edges whose derived nodes are the tips of the tree;
    ##this way we cut down on the number of nodes the program analyses in the next step
    edges_oi <- subset(edges_oi,edges_oi[,2]<length(tree_pruned$tip.label))
    int_nodes <- edges_oi[,1]
    if(identical(int_nodes,numeric(0))){
      print("cannot prune tree with this threshold")
      return(tree_pruned)
    }
      
    for(int_node in int_nodes){
      node_height_oi <- nodeheight(tree_pruned, int_node)
      all_desc <- getDescendants(tree_pruned, int_node)
      leaves <- all_desc[which(all_desc %in% 1:n_tips)]
      dist_leaves <- sapply(leaves, function(x) nodeheight(tree_pruned,x)-node_height_oi)
      node_oi <- int_node
      if(length(which(dist_leaves<dist_cutoff))>1) break  #break if you find a node with > 1 descendant worth pruning
    }
    #print(paste0("removing descendants of ", as.character(node_oi)))
    
    #   node_height_oi <- nodeheight(tree_pruned, node_oi)
    #   all_desc <- getDescendants(tree_pruned, node_oi)
    #   leaves <- all_desc[which(all_desc %in% 1:n_tips)]
    #   dist_leaves <- sapply(leaves, function(x) nodeheight(tree_pruned,x)-node_height_oi)
    
    keep_leaves <- leaves[which(dist_leaves>dist_cutoff)]
    keep_leaves <- c(keep_leaves, leaves[sample(which(dist_leaves<dist_cutoff), 1)]) #sample one from within cutoff leaves
    drop_leaves <- setdiff(leaves, keep_leaves)
    #print("dropped leaves...")
    #print(as.character(tree_pruned$tip.label[drop_leaves]))
    tree_pruned <- drop.tip(tree_pruned, drop_leaves)
    # plot(tree_pruned)
    if(int_node==tail(int_nodes,1)) break #break if you have scanned all nodes and did not find a node with dscendants within cutoff
  }
  return(tree_pruned)
}