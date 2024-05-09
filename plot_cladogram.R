
library(data.tree)
library(treemap)
library(ape)
library(ggtree)
library(tidytree)
library(ggplot2)
library(gridExtra)

rm(list = ls())
gc()

tab <- read.table("Data/rel_abundance.tab", sep = "\t", comment.char = "#")
colnames(tab) <- c("Taxonomy0", "C. suffusus", "C. vittatus")
tab$Taxonomy <- tab$Taxonomy0
head(tab)
tab <- tab[grepl("k__", tab$Taxonomy),]
head(tab)
dim(tab)

# Filter out taxa
tab <- tab[rowSums(tab[, c("C. suffusus", "C. vittatus")] > 0.005) >= 1,]
dim(tab)

# Process taxonomy string
tab$Taxonomy <- gsub("_+", "_", tab$Taxonomy)
head(tab)

tab$Taxonomy <- gsub(";s_", "-", tab$Taxonomy)
head(tab)

tab$Taxonomy <- gsub(";g_", ";", tab$Taxonomy)
head(tab)

# Get pathString to build the tree
tab$pathString <- gsub(";", "/", tab$Taxonomy)
head(tab)

tab$phylum <- gsub("^.+\\/p_|\\/c_.+", "", tab$pathString) #gsub("^.+\\/p_|\\/c_.+", "", tab$pathString)
tab$tip <- gsub("^.+\\/", "", tab$pathString)
head(tab)

# tab <- tab[!grepl("s_NA", tab$V1),]

# Get the top N 
topN <- 5

tab <-  tab[order(apply( tab[, c("C. suffusus", "C. vittatus")], 1, function(x) diff(range(x))), decreasing = TRUE),]
head(tab)
selIds_diffRange <- tab$tip[seq(topN)]

tab <- tab[order(tab[, "C. suffusus"], decreasing = TRUE),]
head(tab)
selIds_CS <- tab$tip[seq(topN)]

tab <- tab[order(tab[, "C. vittatus"], decreasing = TRUE),]
head(tab)
selIds_CV <- tab$tip[seq(topN)]

# ids_iram_patt <- "monas|putida|pitt|oxytoca|nematodiphila"
# ids_spec_iram <- tab$tip[grepl(ids_iram_patt, tab$tip)]

# ids_spec_iram <- c(paste(, ))


selIds_list <- list(diffRange = selIds_diffRange, 
                    "C. suffusus" = selIds_CS, 
                    "C. vittatus" = selIds_CV,
                    ids_spec_iram = list("C. suffusus" = c("Pseudomonas-NA", "Klebsiella-oxytoca", "Acinetobacter-calcoaceticus-pittii"), 
                                         "C. vittatus" = c("Pseudomonas-flavescens-putida", "Serratia-marcescens-nematodiphila")))


# tab$top_size <- ifelse(tab$tip %in% selIds, 5, 1)
# tab$shape <- ifelse(tab$V2 > tab$V3, 19, 15)

# cs_ids <- tab$tip[tab$shape == 19]
row.names(tab) <- tab$tip

taxonomic_tree <- as.Node(tab)
taxonomic_tree_phylo <- as.phylo.Node(taxonomic_tree)


cls <- split(tab$tip, f = tab$phylum)

tree <- groupOTU(taxonomic_tree_phylo, cls)



# tree$edge.length <- rep(10, length(tree$edge.length ))


pdf(file = "test.pdf", width = 14, height = 12)

g_list <- list()
for (n in names(selIds_list)[4]) {
  
  ids <- selIds_list[[n]]
  
  g <- ggtree(tree, layout = "circular", size = 1, ladderize = FALSE) + 
    geom_point(alpha = 0) + 
    geom_tiplab(hjust = -0.15, size = 2.85, aes(color = group), geom = "text", fontface = 3) + 
    xlim(0, 320) + #ylim(0, 100) +
    guides(color = guide_legend(title = "Phylum")) +
    geom_tippoint(size = ifelse(tip.label(tree) %in% unlist(ids), 3.5, 1),
                  colour = ifelse(tip.label(tree) %in% unlist(ids), 
                                  ifelse(tip.label(tree) %in% ids[["CS"]], "yellow2", "red2"), "black")) +
    theme(axis.text = element_text(face = "italic")) 
    # geom_text(show.legend = FALSE)
    # geom_tiplab("a") +
    # ggtitle(n)
  
  print(g)
  
  g_list[[n]] <- g
  
}

# library(reshape2)
tab2heatmap <- reshape2::melt(data = tab[, c("tip", "C. suffusus", "C. vittatus")])#, id = c("CS", "CV"))
tab2heatmap$value <- signif(tab2heatmap$value, 4)
tab2heatmap$percentage <- tab2heatmap$value * 100

g_heatmap <- ggplot(data = tab2heatmap, aes(x = variable, y = tip, fill = percentage))  + 
  xlab("") +
  ylab("") +
  geom_tile(colour = "gray45")  +
  geom_text(aes(x = variable, y = (tip), label = percentage)) +
  scale_fill_gradient(low = "white", high = "steelblue4")  + theme_minimal() +
  theme(axis.text = element_text(size = 10), 
        axis.text.x = element_text(face = "italic", angle = 45, vjust = 0.6),
        axis.text.y = element_text(face = "italic")) +
  labs(fill = "Percentage (%)")

print(g_heatmap)


# library(gridExtra)

gridExtra::grid.arrange(arrangeGrob(g_list[["C. suffusus"]], g_list[["C. vittatus"]], ncol = 1), 
                        g_heatmap,
                        ncol = 2) 


gridExtra::grid.arrange(g_list[["ids_spec_iram"]], 
                        g_heatmap,
                        ncol = 2) 

dev.off()



# References ----
# https://yulab-smu.top/treedata-book/chapter4.html
# https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html
# https://cran.r-project.org/web/packages/data.tree/vignettes/data.tree.html




















