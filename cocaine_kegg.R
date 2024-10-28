# install.packages(c("openxlsx", "ggplot2", "stringr"))
# BiocManager::install(c("clusterProfiler", "enrichplot", "GOplot", "DOSE","org.Hs.eg.db"))
library(openxlsx)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(GOplot)
library(DOSE)
library(stringr)
library("org.Hs.eg.db")
getwd()
setwd("C:/Users/86188/Desktop/geo2")
data <- read.csv('./deg_ppi_nodes_large.csv')

# Convert gene SYMBOLS to ENTREZID using 'bitr' function
gene<-bitr(data$Gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
# write.table(gene,file = './symbol_entrezid.csv',sep=",",row.names = FALSE)
# Perform KEGG pathway enrichment analysis
KEGG<-enrichKEGG(
  gene$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)

write.table(KEGG@result,file = './kegg_pathways.csv',sep=",",row.names = FALSE)

p <- dotplot(KEGG,color='pvalue')
p <- p + theme(axis.text.x = element_text(size = 32),   # Size for x tick labels
               axis.text.y = element_text(size = 32),
               axis.title.x = element_text(size =32),
               legend.text = element_text(size=32),
               legend.title = element_text(size=32),
) + scale_x_continuous(breaks = c(0.075, 0.175), labels = c("0.075", "0.175"))
 

ggsave(
  filename = "dotplot_tmp.png", 
  width = 10,             
  height = 10,            
  units = "in",         
  dpi = 600              
)
# enrichplot::cnetplot(KEGG,circular=TRUE,colorEdge = TRUE)

# Create a heatmap of the KEGG enrichment results
# enrichplot::heatplot(KEGG,showCategory = 10)

# Preprocess data for GOplot library
KEGGplotIn<-KEGG[1:10,c(1,2,5,8)]
KEGGplotIn$geneID <-str_replace_all(KEGGplotIn$geneID,'/',',')
names(KEGGplotIn)<-c('ID','Term','adj_pval','Genes')
KEGGplotIn$Category = "KEGG"
rownames(gene) <- gene$SYMBOL

# Read differential expressed genes (DEGs) data from CSV
genedata<- read.csv('DEGs.csv')
genedata <- genedata[,c(7,1)]
colnames(genedata) <- c('ID','logFC')

# Match DEGs SYMBOLS to ENTREZID
ENTREZID <- as.character(gene[genedata$ID,'ENTREZID'])
ENTREZID[is.na(ENTREZID)]=""
genedata$ID2 <- ENTREZID
genedata <- genedata[genedata['ID2'] != '',]
genedata <- genedata[,c(3,2)]
colnames(genedata) <- c('ID','logFC')

# Prepare data for GOplot circle plot
circ<-GOplot::circle_dat(KEGGplotIn,genedata)

# # Filter genes of interest for chord plot
# genedata2 <- genedata[genedata$ID %in% c('2353','1958','3725','3569'),]
chord<-chord_dat(data = circ,genes = genedata)
symbol <- gene[match(rownames(chord),gene$ENTREZID),'SYMBOL']
rownames(chord) <- symbol
# 
# # Create a GO chord plot
# p <- GOChord(
#   data = chord,
#   title = '',
#   space = 0,
#   limit = c(1,1),
#   gene.order = 'logFC',
#   gene.space = 0.25,
#   gene.size = 20,
#   lfc.col = c('red','white','blue'), #Colors for Up and down genes
#   #ribbon.col = brewer.pal(length(GOplotIn$Term)),#GO Term colors
#   process.label = 15 #GO Term Font Size
# )
# 
# 
# # Save the GO chord plot as a high-resolution PNG image
# ggsave(
#   filename = "chord.png", 
#   width = 30,             
#   height = 30,            
#   units = "in",         
#   dpi = 300              
# )


GOCircle(circ,table.legend='True')

# GOCircle(circ, lfc.col = c('purple', 'orange'))
# GOCircle(circ, zsc.col = c('yellow', 'black', 'cyan'))
# 
# GOCircle(circ, label.size = 5, label.fontface = 'italic')
write.csv(KEGGplotIn,'table_circle.csv')
ggsave(
  filename = "circle_tmp.png", 
  width = 16,             
  height = 10,            
  units = "in",         
  dpi = 600              
)


data(EC)

GOCluster(circ, KEGGplotIn$Term)
p<-GOCluster(circ,KEGGplotIn$Term,clust.by='logFC',lfc.col=c('darkgoldenrod1','black','cyan1'))
p <- p + theme(
  text = element_text(size = 12),
  # legend.position = "Bottom", 
  legend.direction = "vertical",
  legend.text = element_text(size = 12),
  # legend.title = element_text(size = 14)
)
ggsave(
  filename = "gocluster_tmp.png", 
  width = 15 ,             
  height = 15,            
  units = "in",         
  dpi = 600              
)
