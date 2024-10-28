# Ensure BiocManager is installed
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Uncomment to install necessary packages
# BiocManager::install("limma")
# BiocManager::install("GEOquery")
# install.packages("devtools")
# devtools::install_github("BioSenior/ggVolcano")
getwd()
setwd("C:/Users/86188/Desktop/geo2")

# Load required libraries
library(GEOquery)
library(limma)
library(ggVolcano)

# Download data from GEO
eSet <- getGEO("GSE54839", destdir = ".", getGPL = T)
exprSet = exprs(eSet[[1]])
fdata = fData(eSet[[1]])
pdata = pData(eSet[[1]])
dim(exprSet)  
dim(fdata)
dim(pdata)

num_samples <- dim(pdata)[1]

# Create sample grouping based on condition
sml <- c()
for (element in pdata$title) {
  if (grepl('con',element)){
    sml <- c(sml,'0')
  }else {
    sml <- c(sml,'1')
  }
}

# Adjust group names based on sample order
if (sml[1] == '0'){
  groups = make.names(c('control','addiction'))
}else {
  groups = make.names(c('addiction','control'))
}

fac <- factor(sml)
levels(fac) <- groups
design <- model.matrix(~fac+0)
colnames(design) <- levels(fac)
rownames(design) <- colnames(exprSet)

# Create contrast matrix for differential expression analysis
ctr <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=ctr, levels=design)

# Merge gene and expression data
genes <- fdata[,c("ID","Symbol")]
exprSet <- data.frame(exprSet)
exprSet$ID <- rownames(exprSet)
exprSet <- merge(exprSet,genes,by.x = "ID",by.y = 'ID')
exprSet <- exprSet[,-1]

# Handle duplicate genes
rowMeans = apply(exprSet[,c(1:num_samples)],1,function(x) mean(as.numeric(x), na.rm = T))
exprSet = exprSet[order(rowMeans, decreasing = T),]
exprSet_2 = exprSet[!duplicated(exprSet[, dim(exprSet)[2]]),] 
exprSet_na = na.omit(exprSet_2)   
explan_final = exprSet_na[exprSet_na$'Symbol' != "",]

# Handle cases where one probe maps to multiple genes
explan_final <- data.frame(explan_final[-grep("/",explan_final$"Symbol"),]) 
rownames(explan_final) <- explan_final$"Symbol"
explan_final <- explan_final[,c(1:num_samples)]

# Convert to log2 scale if necessary
qx <- as.numeric(quantile(explan_final,c(0.,0.25,0.5,0.75,0.99,1.0),na.rm=T))
LogC <- (qx[5]>100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) {
  explan_final[which(explan_final <= 0)] <- NaN
  explan_final <- log2(explan_final)
}
explan_final <- explan_final[complete.cases(explan_final),]

# Differential expression analysis
fit <- lmFit(explan_final, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tT2 <- topTable(fit2, sort.by="B", adjust.method="fdr", number=Inf)
tT2$geneName <- rownames(tT2)
data <- add_regulate(tT2, log2FC_name = "logFC", fdr_name = "P.Value",log2FC = 0.2630344, fdr = 0.05)

## test
##write.table(tT2,file = './DEGs.csv',sep=",",row.names = FALSE)

p <- ggvolcano(data, x = "log2FoldChange", y = "padj", 
               fills = c("#e94234","#b4b4d8","#269846"),
               colors = c("#e94234","#b4b4d8","#269846"),
               label = "row", label_number= 0, 
               output = FALSE, log2FC_cut = 0.2630344, FDR_cut = 0.05,legend_position = "UR",
               x_lab = '',y_lab = ''
)
p <- p + theme(
  legend.title = element_blank(),
  axis.text.x = element_text(size = 14,colour = 'black'),
  legend.text = element_text(size=14),# Size for x tick labels
               axis.text.y = element_text(size = 14,colour = 'black'),
               
               legend.position = 'top')
# p <- p + theme(
#   axis.line = element_line(colour = "black", size = 2)  # Adjust thickness for both x and y axes
# )
p + theme(plot.background = element_rect(fill = "transparent", 
))

ggsave( filename = "Cocaine_volcano_tmp.png", width = 4, height = 4, units = "in", dpi =
          600 ,bg = "transparent")




# # Create a volcano plot
# ggvolcano(data, x = "log2FoldChange", y = "padj", label = "row", label_number
# = 10, output = FALSE, log2FC_cut = 0.2630344, FDR_cut = 0.05, y_lab='Log10
# P-Value')
# #




ggsave( filename = "volcano.png", width = 15, height = 10, units = "in", dpi =
300 ) 
# Write differentially expressed genes (DEGs) to CSV 
regulate_data <- data[data$regulate != 'Normal',]
colnames(regulate_data)[colnames(regulate_data) == "padj"] <- "PVal"
write.table(regulate_data,file = './DEGs.csv',sep=",",row.names = FALSE)
