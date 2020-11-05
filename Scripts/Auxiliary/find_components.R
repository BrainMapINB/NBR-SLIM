find_components <- function(all_edges, thr, tri_pos, nnodes){
  
  # Find significant edges
  edges <- which(all_edges > thr)
  if(length(edges) > 0){
    # Initializate objects
    node_comp <- 1:nnodes
    lista <- vector("list",1)
    component <- vector("integer", length(edges))
    # Store edge strength above threshold
    strength <- all_edges[edges] - thr
    # Find components
    for(jj in 1:length(edges)) component[jj] <- node_comp[tri_pos[edges[jj],]] <- min(node_comp[tri_pos[edges[jj],]])
    # Store results
    lista <- cbind(edges, tri_pos[edges,1], tri_pos[edges, 2], component, strength)
    colnames(lista) <- c("2Dcol", "3Drow","3Dcol", "comp", "strn")
  } else lista <- NULL
  return(lista)
}