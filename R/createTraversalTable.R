#' Returns the leaves of a binary tree and a post-order traversal
#'
#' @param tree a phylo object from the ape package
#' @return a list containing the tree tips and a postorder traversal
#' @export
createTraversalTable <- function(tree){
  tree = ape::reorder.phylo(tree, "postorder")
  traversal=as.matrix(data.table::data.table(id=1:nrow(tree$edge),tree$edge)[,.(child.id=V2,edge.id=id),by="V1"])
  tips=setdiff(tree$edge[,2],tree$edge[,1])  
  return(list(tips=tips,traversal=traversal))  
}
