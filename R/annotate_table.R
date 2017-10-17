NULL

#' Assigns unassigned taxa in a clade to a group based on assignment of other taxa from the table

#' @param phy A rooted tree object.
#' @param annot_table A data frame where each row has a taxa and corresponding groups. Groups must be in column labeled `Phylo.group` and the first column must be taxa name in `phy`.
#'
#' @return A dataframe where all the taxa are now assigned to a group based on monophylyl.
#'
#' @export


annotate_table <- function(phy, annot_table) {
    #check if tree is rooted.
    if ( !ape::is.rooted(phy) ) {
        stop(paste("Specified tree is not rooted, please root the tree and try again!"))
    }

    if ( !all(phy$tip.label %in% annot_table[, 1])) {
        stop(paste("Specified tree taxa and table's first column do not match, fix, and try again!"))
    }

  # get list of groups from the annotation table
  phylo_groups <- as.vector(unique(annot_table$Phylo.group))
  final_table <- data.frame(c(tip = c(), clade.name = c()))
  for (group in phylo_groups) {
    #  get names of taxa that are within the group
    taxa_in_group <- as.vector(subset(annot_table, Phylo.group == group  )$Organism.name.strain)
    # get the most recent common ancestor of those taxa
    MRCA = phytools::findMRCA(phy, tips = taxa_in_group)
    # if the taxa is not null
    if (!is.null(MRCA)) {
      # get all the tip label associated with it
      all_tips = geiger::tips(phy, MRCA)
      clade_name <- data.frame(tip = all_tips, clade.name = rep(group, length(all_tips)))
      final_table <- rbind(clade_name, final_table)
    }
    else if (is.null(MRCA)) {
      clade_name <- data.frame(tip = taxa_in_group, clade.name = taxa_in_group)
      final_table <- rbind(clade_name, final_table)
    }
  }

  df1 <- final_table %>% dplyr::group_by(tip) %>% dplyr::arrange(dplyr::desc(clade.name)) %>% dplyr::filter(row_number() == 1) %>% data.frame()

  # taxa that were not assigned to monophyletic groups are assigned to group corresponding to their name
  df1$clade.name <- ifelse(df1$clade.name == "Unknowns", as.character(df1$tip), as.character(df1$clade.name))
  return(df1)
}
