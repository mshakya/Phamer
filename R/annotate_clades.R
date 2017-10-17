#' @importFrom ape nodelabels

NULL

#' Plot and label clades of a tree

#' @param phy A tree object with bootstrap support values as node labels.
#' @param annot_table A data frame that has each taxa in monophyletic clades assigned to a group. The first column name must be taxa name and the second column must be clade names and name clade.name.
#'
#'
#' @return A plot object
#'
#' @export

annotate_clades <- function(phy, annot_table){

    # check if taxa name of the tree matches with the first column
    if ( !all(phy$tip.label %in% annot_table[, 1])) {
        stop(paste("Specified tree taxa and table's first column do not match, fix, and try again!"))
    }

    # check if second column of the data frame is named clade.name
    if ( colnames(annot_table)[2] != "clade.name" ) {
        stop(paste("Specified data frame second column is not named clade.name, fix, and try again!"))
    }


    # convert to a ggtree object
    ggtree_phy <- ggtree::ggtree(phy)

    phylo_groups <- as.vector(unique(annot_table$clade.name))

    for (group in phylo_groups) {
        #  get names of taxa that are within the group
        taxa_in_group <- as.vector(subset(annot_table, clade.name == group  )[, 1])
        if (length(taxa_in_group) > 1) {
        # get the most recent common ancestor of those taxa
        MRCA <- phytools::findMRCA(phy, tips = taxa_in_group)
        # if the taxa is not null

        ggtree_phy <- ggtree_phy + ggtree::geom_cladelabel(node = MRCA, label = group, align = T)
}
}

    return(ggtree_phy + ggtree::geom_tiplab())
    }

