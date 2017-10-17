
#
# read the tree file
t <- ape::read.tree("../test_data/EschShigSalmRef.fasttree")

# read the table with metadata
a <- read.csv("../test_data/Eschericia_phylotypes.csv", header = T)

# root the tree
troot <- ape::root(t, node = 697, resolve.root = T)

# annotate the taxa in the table
df <- annotate_table(troot, a)

subset_t <- backbone_tree(df, troot, size=F)
expect_equal(length(subset_t$tip.label), 59)
