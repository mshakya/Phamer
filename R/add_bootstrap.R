#' @importFrom ape nodelabels

NULL

#' An extension of nodelabels to add symbols indicating the bootstrap support pass certain threshold

#' @param phy A tree object with bootstrap support values in nodes.
#' @param threshold A threshold value. Only the nodes that have support higher than this number will be plotted.
#' @param pch Shape to add
#'
#'
#' @return plots shapes indicating bootstrap support in the current plot window
#'
#' @export

add_bootstrap <- function(phy, pch, threshold){

    interpret_bs <- unlist(lapply(phy$node.label, function(x) if (x == "Root") {x} else if (x == "") {""} else if (as.integer(x) > threshold)  {"bs"} else {""}))

    node_pos <- which(interpret_bs %in% c("bs"))

    print(node_pos)
    ape::nodelabels(node = node_pos, pch = pch)
}

