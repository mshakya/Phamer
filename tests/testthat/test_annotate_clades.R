
#
# read the tree file
t <- ape::read.tree("../test_data/EschShigGroup.tree")

# read the table with metadata
a <- read.csv("../test_data/genome_phylogroups.csv", header = T)

# root the tree
troot <- ape::root(t, outgroup = "E__fergusonii_ATCC_35469", resolve.root = T)

# annotate clades of trees
annotated_ggtree <- annotate_clades(troot, a)

expect_equal(dim(annotated_ggtree$data[1]), 69)
