# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


makebackbone <- function(df, phy, plot_type){
  # input: A dataframe with two columns, first column has taxa name matching
  # trees, and second column has its group affiliation
  # input: A tree file

  # collapses monophyletic clade if they belong to same group and draws a backbone tree

  # ouputs a tree object that can be plotted to a tree

  #############################################################################
  # get list of taxa name
  listofspecies <- phy$tip.label
  # empty list
  listtopreserve <- list()
  # empty list
  lengthofclades <- list()
  # empty list
  meandistnode <- list()
  # empty list
  newedgelengths <- list()
  #############################################################################

  # extract all groups as vectors from the table
  groupings <- as.vector(unique(df$clade.name))
  # loop through the groupings
  for (i in 1:length(groupings)){

    # get all the taxa name associated with that grouping
    taxa_in_group <- as.vector(subset(df, clade.name == groupings[i] )$tip)

    # find the MRCA node number of that group
    bestmrca <- getMRCA(phy, taxa_in_group )

    # if that taxa has more than 1 taxa (collapse)
    if (length(taxa_in_group) > 1){

      # get all descendants of MRCA
      mrcatips <- phy$tip.label[unlist(phangorn::Descendants(phy,bestmrca, type="tips") )]

      # keep the first taxa in the list by assigning it to list to preserve
      listtopreserve[i] <- mrcatips[1]

      # get the mean distance to MRCA from all the taxa in this group
      meandistnode[i] <- mean(dist.nodes(phy)[unlist(lapply(mrcatips,
                                                            function(x) grep(x, phy$tip.label) ) ),bestmrca] )
    }
    # if there is only one taxa that belongs to the group
    else if (length(taxa_in_group) == 1) {
      mrcatips <- taxa_in_group
      listtopreserve[i] <- mrcatips
      # just get distance from that taxa from its parent
      meandistnode[i] <- setNames(phy$edge.length[sapply(1:n,function(x,y)   which(y == x),y = phy$edge[,2])], mrcatips)
    }


    # number of genome in clades, this will determine the size
    lengthofclades[i] <- length(mrcatips)

    # remove, all but one branch that belongs to that group
    provtree <- drop.tip(phy, mrcatips, trim.internal = F, subtree = T)

    # number of taxa in tree after dropping branches
    n3 <- length(provtree$tip.label)

    # assign the new edge length
    newedgelengths  <- setNames(provtree$edge.length[sapply(1:n3,function(x,y)
      which(y == x),y = provtree$edge[,2])],
      provtree$tip.label)[provtree$tip.label[grep("tips",provtree$tip.label)] ]

  }

  # drop tips
  newtree <- drop.tip(phy, setdiff(listofspecies,unlist(listtopreserve)),
                      trim.internal = T)
  n <- length(newtree$tip.label)
  newtree$edge.length[sapply(1:n,function(x,y)
    which(y == x),y = newtree$edge[,2])] <- unlist(newedgelengths) + unlist(meandistnode)


  if (plot_type == "size") {
    trans <- data.frame(tip.label = newtree$tip.label,clade.label = groupings,
                        N = unlist(lengthofclades), depth = unlist(meandistnode) )
  }

  else if (plot_type == "same") {
    trans <- data.frame(tip.label = newtree$tip.label,clade.label = groupings,
                        N = rep (10, length(groupings)), depth = unlist(meandistnode) )
  }


  rownames(trans) <- NULL
  backboneoftree <- phytools::phylo.toBackbone(newtree,trans)
  return(backboneoftree)
}

