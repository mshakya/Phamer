#' @importFrom ape drop.tip

NULL

#' Drops, all but one taxa from a monophyletic clade, keeps the info on total number of taxa

#' @param df A dataframe with two columns, first column has taxa name matching phy, and second column has its group affiliation.
#' @param phy A tree object.
#' @param backbone If set to TRUE return backbone tree, else returns a regular tree with dropped branches. (default=TRUE)
#' @param size If set to TRUE draws clades proportional to number of taxa. If FALSE, draws clades that are of same size. (default=TRUE)
#'
#'
#' @return A backbone tree with
#'
#' @export

backbone_tree <- function(df, phy, size=TRUE, backbone=TRUE){

    #############################################################################
    # empty list
    tips_to_drop <- list()
    #############################################################################


    groupings <- as.vector(unique(df$clade.name))
    trans <- data.frame(taxa = c(), size = c(), clade_name = c())

    # loop through the groupings
    for (i in 1:length(groupings)) {

        # get all the taxa name associated with that grouping
        taxa_in_group <- as.vector(subset(df, clade.name == groupings[i] )$tip)

        if (length(taxa_in_group) > 1) {

        tips_to_drop_from_group <- taxa_in_group[1:length(taxa_in_group) - 1]

        tips_to_drop[[i]] <- tips_to_drop_from_group

        tips_to_keep <- setdiff(taxa_in_group, tips_to_drop_from_group)

        meandistnode <- setNames(phy$edge.length[sapply(1:length(phy$tip.label), function(x, y) which(y == x), y = phy$edge[,2])], phy$tip.label)[tips_to_keep][[1]]

        if (size) {

        trans <- rbind(trans, data.frame(tip.label = tips_to_keep, clade.label = groupings[i], N = length(taxa_in_group), depth = meandistnode))
        }
        else {
        trans <- rbind(trans, data.frame(tip.label = tips_to_keep, clade.label = groupings[i], N = 10, depth = meandistnode))
        }
        }

    }

    trimmed_t <- ape::drop.tip(phy, tip = unlist(tips_to_drop ))

    backboneoftree <- phytools::phylo.toBackbone(trimmed_t, trans)
    if (backbone) {
        return(backboneoftree)
    }
    else{
        return(trimmed_t)
    }

}

