# read the tree
t <- ape::read.tree("../test_data/EschShigSalmRef.fasttree")
# read the csv
a <- read.csv("../test_data/Eschericia_phylotypes.csv", header = T, sep = ",")
# check for error
expect_error(annotate_table(t, a), "Specified tree is not rooted, please root the tree and try again!")
# root the tree
troot <- ape::root(t, node=697, resolve.root = T)
# create annotation table
df <- annotate_table(troot, a)
#see if the annotation table has correct dimension
expect_that(dim(df)[1], equals(629))
