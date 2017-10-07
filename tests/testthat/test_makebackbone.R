
#
# read the tree file
t <- ape::read.tree("../test_data/EschShigSalmRef.fasttree")

# read the table with metadata
a <- read.csv("../test_data/Eschericia_phylotypes.csv", header = T)

# root the tree
troot <- ape::root(t, node = 697, resolve.root = T)

# annotate the taxa in the table
df <- annotate_table(troot, a)

backbone_tree <- makebackbone(df, t, "same")
pdf("Rplot.pdf")
plot(backbone_tree)
dev.off()
expect_equal(backbone_tree$Nnode, 58)
