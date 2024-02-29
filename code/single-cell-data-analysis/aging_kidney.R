#--------------------------------------------------------------------------------------------------------
### aging kidney analysis R script
#--------------------------------------------------------------------------------------------------------

library(data.table)
setwd("/aging/")
file_list <- list.files()
print(file_list)

sampleNames = mapply(function(x) x[1], strsplit(file_list, split = '__'))
sampleNames
sampleNames <- gsub('_','_',sampleNames)
meta.info = data.table(file_Path = file_list, sampleID = sampleNames)
meta.info[, SubType:=mapply(function(x) x[1], strsplit(file_list, '-'))]
DT::datatable(meta.info)

norm_method = 'standard'
treat = FALSE
show = FALSE
### set the threshold of MT genome percentage: 20%
MT_cut = 100 

batch_list = list()

for(i in 1:nrow(meta.info)){
  dir_of_10X = meta.info$file_Path[i]
  sampleID = meta.info$sampleID[i]
  subtype = meta.info$SubType[i]
  seq_batch = meta.info$Batch[i]
  print(paste("Starting Processing", sampleID, 'at', Sys.time()))
  sc_mat = Read10X(dir_of_10X)
  #sc_mat<-sc_mat[rownames(sc_mat)%in%lnc_list2$lnc_list2,]
  # colnames(sc_mat) = make.names(colnames(sc_mat), unique = T)
  
  sc_tumor <- CreateSeuratObject(counts = sc_mat, project = sampleID, min.cells = 3, min.features = 0)
  sc_tumor
  
  sc_tumor[["percent.mt"]] <- PercentageFeatureSet(sc_tumor, pattern = "^mt-")
  
  #qc = VlnPlot(sc_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
#  sc_tumor <- subset(sc_tumor, subset = nFeature_RNA < 4000  & percent.mt <= 30 & nCount_RNA < 20000 )
  
  if(norm_method == 'SCT'){
    sc_tumor <- SCTransform(sc_tumor, vars.to.regress = 'percent.mt')
  } else{
    sc_tumor <- NormalizeData(sc_tumor, verbose = F)
    sc_tumor <- FindVariableFeatures(sc_tumor, selection.method = "vst", 
                                     nfeatures = 2000, verbose = FALSE)
    
    sc_tumor <- ScaleData(sc_tumor)
  } 
  
  
  sc_tumor <- RenameCells(object = sc_tumor, add.cell.id = sampleID)  # add sample name as prefix
  if(treat == TRUE){
    sc_tumor <- RunPCA(sc_tumor, verbose = T)
    sc_tumor <- RunUMAP(sc_tumor, dims = 1:30, verbose = FALSE)
    sc_tumor <- RunTSNE(sc_tumor, dims = 1:30, verbose = FALSE)
    
    sc_tumor <- FindNeighbors(sc_tumor, dims = 1:30, verbose = FALSE)
    sc_tumor <- FindClusters(sc_tumor, verbose = FALSE)
    
    umap = DimPlot(sc_tumor, label = TRUE) + NoLegend()
    tsne = TSNEPlot(sc_tumor, label = TRUE) + NoLegend()
    
    if(show == TRUE){
      #pdf(paste0('~/Project/sc_ATC/figure/QC/',sampleID,'_Seurat_QC.pdf'),12,7.5)
      print(qc)
      plot_grid(plot1,plot2)
      dev.off()
      #pdf(paste0('~/Project/sc_ATC/figure/Clustering/',sampleID,'_Seurat_Dimension_reduction_and_Clustering.pdf'))
      print(umap)
      print(tsne)
    }
  }
  
  #dev.off()
  ### add information
  sc_tumor@meta.data[,'SubType'] = subtype
  sc_tumor@meta.data[,'SampleID'] = sampleID
  sc_tumor@meta.data[,'Batch'] = seq_batch
  ### add batch data into a list to aging_kidney
  #assign(paste0(sampleID,'.Seurat'),sc_tumor)
  batch_list <- append(batch_list, sc_tumor)
  print(paste("Finishing Proliferationcessing", sampleID, 'at', Sys.time()))
}


aging_kidney <- Reduce(function(x,y){merge(x,y,merge.data = TRUE)},batch_list)
aging_kidney <- subset(aging_kidney, subset = nFeature_RNA < 4000  & percent.mt <= 30 & nCount_RNA < 20000 )
VlnPlot(aging_kidney, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0, group.by = "orig.ident")


#------------------------------------------------------------------------------------------------------------------------------
### identical cells isolation from pre-capture and post-capture
#------------------------------------------------------------------------------------------------------------------------------

a<-subset(aging_kidney, cells=names(aging_kidney@active.ident[aging_kidney$orig.ident%in%orig[1]]))
a2<-subset(aging_kidney, cells=names(aging_kidney@active.ident[aging_kidney$orig.ident%in%orig[1 + 8]]))

a_col = mapply(function(x) x[3], strsplit(colnames(a), split = '_'))
a2_col = mapply(function(x) x[3], strsplit(colnames(a2), split = '_'))

a3<-a_col[a_col%in%a2_col]
length(a3)

a3_col <- subset(a, cells = names(a@active.ident[mapply(function(x) x[3], strsplit(colnames(a), split = '_')) %in% a3]))
a4_col <- subset(a2, cells = names(a2@active.ident[mapply(function(x) x[3], strsplit(colnames(a2), split = '_')) %in% a3]))
dim(a3_col)
dim(a4_col)

a4<-rbind(as.matrix(colnames(a3_col)),as.matrix(colnames(a4_col)))

for (i in c(2:8)) {
  
  a<-subset(aging_kidney, cells=names(aging_kidney@active.ident[aging_kidney$orig.ident%in%orig[i]]))
  a2<-subset(aging_kidney, cells=names(aging_kidney@active.ident[aging_kidney$orig.ident%in%orig[i + 8]]))
  
  a_col = mapply(function(x) x[3], strsplit(colnames(a), split = '_'))
  a2_col = mapply(function(x) x[3], strsplit(colnames(a2), split = '_'))
  
  a3<-a_col[a_col%in%a2_col]
  print(length(a3))
  
  a3_col <- subset(a, cells = names(a@active.ident[mapply(function(x) x[3], strsplit(colnames(a), split = '_')) %in% a3]))
  a4_col <- subset(a2, cells = names(a2@active.ident[mapply(function(x) x[3], strsplit(colnames(a2), split = '_')) %in% a3]))
  dim(a3_col)
  dim(a4_col)
  
  a4<-rbind(as.matrix(colnames(a3_col)),as.matrix(colnames(a4_col)),a4)
  
}

aging_kidney<-subset(aging_kidney, cells=names(aging_kidney@active.ident[colnames(aging_kidney) %in% a4]))


#---------------------------------------------------------------------------------------------------------------------------------------------------------
### doublelet romove
#---------------------------------------------------------------------------------------------------------------------------------------------------------

library(DoubletFinder)

orig <- names(table(aging_kidney$orig.ident))
double<-list("a1",	"a2",	"a3",	"a4",	"a5",	"a6",	"a7",	"a8")
table(aging_kidney$orig.ident)


# chagne the argument value according to the sample

a<-subset(aging_kidney,cells = names(aging_kidney@active.ident[aging_kidney$orig.ident%in%orig[16]]))
a <- NormalizeData(a, verbose = F)
a <- FindVariableFeatures(a, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
a <- ScaleData(a, vars.to.regress = c("nCount_RNA","nFeature_RNA", "percent.mt"), verbose = TRUE)
a <- RunPCA(a, npcs = 50, verbose = TRUE)
ElbowPlot(object = a, ndims = 50)
a <- RunUMAP(a, dims = 1:30, verbose = TRUE)
sweep.res.list_kidney <- paramSweep_v3(a, PCs = 1:30, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)
nExp <- round(ncol(a) * 0.011)  # expect 4% doublets
a <- doubletFinder_v3(a, pN = 0.02, pK = 0.18, nExp = nExp, PCs = 1:30)
DF.name = colnames(a@meta.data)[grepl("DF.classification", colnames(a@meta.data))]

double[[8]]<-a

DimPlot(a,label = T)+ NoLegend()

a <- double[[1]] 
b<-a@meta.data[a@meta.data[,8]%in%"Doublet",]
d<-as.matrix(b)

for(i in 2:8){
   a <- double[[i]] 
   b<-a@meta.data[a@meta.data[,8]%in%"Doublet",]
   c<-as.matrix(b)
   d<-rbind(c,d)
   table(d[,1])
}


aging_kidney<-subset(aging_kidney,cells = names(aging_kidney@active.ident[!colnames(aging_kidney)%in%rownames(d)]))

#-------------------------------------------------------------------------------------------------------- 
### sum pre-capture and post-capture
#--------------------------------------------------------------------------------------------------------

orig<-names(table(aging_kidney$orig.ident))
b2<-list("a1",	"a2",	"a3",	"a4",	"a5",	"a6",	"a7",	"a8")

for (i in 1:8) {  
  a<-subset(aging_kidney, cells=names(aging_kidney@active.ident[aging_kidney$orig.ident%in%orig[i]]))
  a2<-subset(aging_kidney, cells=names(aging_kidney@active.ident[aging_kidney$orig.ident%in%orig[i + 8]]))
  
  b <- (as.matrix(a2@assays$RNA@counts)+ as.matrix(a@assays$RNA@counts))
  b2[i] <-  CreateSeuratObject(b)  
}

aging_kidney<-merge( x = b2[[1]], y = c(b2[[2]],b2[[3]],b2[[4]],b2[[5]],b2[[6]],b2[[7]],b2[[8]]))
dim(aging_kidney)
table(aging_kidney$orig.ident)


#-------------------------------------------------------------------------------------------------------- 
### single cell anlaysis workflow (Seurat)
#-------------------------------------------------------------------------------------------------------- 

library(Seurat)
aging_kidney <- NormalizeData(aging_kidney)
aging_kidney <- FindVariableFeatures(aging_kidney, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
aging_kidney[["percent.mt"]] <- PercentageFeatureSet(aging_kidney, pattern = "^mt-")
VlnPlot(aging_kidney, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, pt.size = 0, group.by = "information")

aging_kidney <- ScaleData(aging_kidney, vars.to.regress = c('nCount_RNA', 'percent.mt', 'nFeature_RNA'), verbose = TRUE)
aging_kidney <- RunPCA(aging_kidney, npcs=50, verbose=T)
ElbowPlot(object = aging_kidney, ndims = 50)
num_PC = 33
aging_kidney <- RunUMAP(aging_kidney, reduction = 'pca', dims=1:num_PC)
aging_kidney <- FindNeighbors(aging_kidney, reduction = 'pca', dims=1:num_PC)
aging_kidney <- FindClusters(aging_kidney, resolution = 1)

DimPlot(aging_kidney, label = T)


#-------------------------------------------------------------------------------------------------------- 
### batch correction - fastMNN
#--------------------------------------------------------------------------------------------------------

library(SeuratWrappers)

head(aging_kidney@meta.data)
table(aging_kidney$information)
aging_kidney.list <- SplitObject(aging_kidney, split.by="information")

for (i in 1:length(aging_kidney.list)) {
  aging_kidney.list[[i]] <- NormalizeData(aging_kidney.list[[i]],assay="RNA",normalization.method="LogNormalize")
  aging_kidney.list[[i]] <- FindVariableFeatures(aging_kidney.list[[i]], selection.method="vst", nfeatures=2000)
}

aging_kidney <- RunFastMNN(object.list = aging_kidney.list)
aging_kidney <- RunUMAP(aging_kidney, reduction = "mnn", dims = 1:11)
aging_kidney <- FindNeighbors(aging_kidney, reduction = "mnn", dims = 1:11)
aging_kidney <- FindClusters(aging_kidney, resolution = 1)
DimPlot(aging_kidney, label = TRUE, label.size = 5)

musGenes <- c(
  "Nphs1","Nphs2", # podocyte
  "Cldn1","Ncam1", # PEC
  "Gata3","Pdgfrb", # mesangial cell
  "Ren1", "Akr1b7", # Juxtaglomerular cells
  "Lrp2", # PT
  "Pax8", "Ppara", # STC
  "Slc5a2", "Slc13a3", # PCT
  "Cyp7b1","Serpina1a", # PST
  "Clcnka", # asc LOH
  "Umod", "Slc12a1", # ALOH
  "Nos1",  # macula densa
  "Slc12a3", "Pvalb", # DCT
  "Hsd11b2", "Aqp2", #PC
  "Atp6v1g3",  "Slc26a4", # IC
  "Pecam1", "Ehd3", "Kdr",  # GEC
  "Alox12", "Ptprr",	# aterioles
  "Mki67", "Top2a", # stem cell
  "Pax2", "Prom1", # stem cell
  "Cd79a", # B cell
  "Cd3d", # lymphocyte
  "Cd4", # Cd4
  "Cd8a", # Cd8
  "Klrb1c", # NK
  "C1qa", "C1qb",   # M1
  "Ccr2", "F13a1",  # M2
  "Cd68", "Lst1", # Macrophage
  "Hba-a2"
)

DotPlot(aging_kidney, col.min = 0, features = musGenes, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

aging_kidney <- RenameIdents(object = aging_kidney, "0"="lymphocytes",	"1"="podocyte",	"2"="GEC",	"3"="PT1",	"4"="PT2",	"5"="PEC",	"6"="JGC",	"7"="endo_podo",	"8"="mesangial",	"9"="myeloid",	"10"="aterioles",	"11"="ALOH",	"12"="CD PC",	"13"="DCT",	"14"="endo_immune",	"15"="endo_ALOH",	"16"="RBC")
aging_kidney@active.ident = factor(aging_kidney@active.ident, levels = c("podocyte",	"mesangial", "JGC",		"PEC",	"PT1",	"PT2",	"ALOH",	"CD PC",	"GEC",	"endo_podo",	"endo_ALOH", "endo_immune",	"lymphocytes",	"myeloid",	"RBC"))

DimPlot(aging_kidney, label = TRUE, label.size = 3)

lnc_marker = c( "2810433D01Rik",	"Wt1os",	"C130021I20Rik",	"Gm15866",	"4921504A21Rik",	"Gm13814",	"Gm43462",	"Dubr",	"Gm26771",	"Carmn",	"Gm12840",	"E030013I19Rik",	"Gm14964",	"Mannr",	"A830082K12Rik",	"Gm42679",	"Gm26632",	"Gm13861",	"Tbx3os1",	"Gm12122",	"9230112J17Rik",	"Gm49542",	"Gm42397",	"Gm43190",	"Gm43581",	"C330002G04Rik",	"Gm19950",	"0610005C13Rik",	"BC049987",	"Gm27216",	"Gm4208",	"Gm44785",	"Pantr1",	"Gm47708",	"Gm17750",	"Gm11789",	"1110019B22Rik",	"Gm12927",	"Gm17546",	"Bvht",	"Gm9917",	"3110099E03Rik",	"4631405J19Rik",	"Plet1os",	"Gm31814",	"AW112010",	"Mir142hg",	"Gm5547",	"Mirt1",	"Gm36161",	"Gm26740",	"Gm43065",	"Gm19585",	"Gm31243",	"AI839979",	"Tnfsf13os",	"Slc36a3os" )

DotPlot(aging_kidney, col.min = 0, features = unique(lnc_marker), cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
  

#-------------------------------------------------------------------------------------------------------- 
### stacked barplot
#--------------------------------------------------------------------------------------------------------

table(aging_kidney)
name <- names(table(aging_kidney@active.ident))

a2 <- subset(aging_kidney, idents = name[1])
a3 <- as.data.frame(table(a2$information))
a3$ident <- name[1]
a4 <- a3

for (i in 2:length(name)) {
  a2 <- subset(aging_kidney, idents = name[i])
  a3 <- as.data.frame(table(a2$information))
  a3$ident <- name[i]
  a4 <- rbind(a4,a3)
}

a4$ident = factor(a4$ident, levels = names(table(aging_kidney@active.ident)))
a4$Var1 = factor(a4$Var1, levels = names(table(aging_kidney$age)))
a4$Var1 = factor(a4$Var1, levels = names(table(aging_kidney$tissue)))
a4$Var1 = factor(a4$Var1, levels = c( "pre-rep1-young_glom", "pre-rep2-young_glom", "pre-rep1-old_glom", "pre-rep2-old_glom", "pre-rep1-young_tubule", "pre-rep2-young_tubule", "pre-rep1-old_tubule", "pre-rep2-old_tubule" ))

cols = c( "#80B1D3", "#FB8072" )
cols = c("red", "navy") 
cols = c("yellowgreen", "red","blue", "gold4", "coral",	"chartreuse",	"yellow2","forestgreen",	"olivedrab",	"palegreen",	"royalblue",	"dodgerblue",	"skyblue",	"violet",	"peachpuff",	"pink2",	"darkgoldenrod1",	"slateblue2",	"lemonchiffon3",	"burlywood3",	"aquamarine3",	"lightblue3",	"deeppink3",	"salmon1",	"darkorange3", "lightgoldenrod2",	"darkorange",	"lightpink1",	"lightcoral",	"hotpink1",	"mediumpurple1",	"slateblue2",	"lemonchiffon3",	"burlywood3",	"aquamarine3",	"lightblue3",	"deeppink3",	"seashell2")


# Stacked + percent
ggplot(a4, aes(fill=Var1, y=Freq, x=ident)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual( values = cols ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 12))


#--------------------------------------------------------------------------------------------------------
### idents barplot
#--------------------------------------------------------------------------------------------------------

a<-subset(aging_kidney, cells = names(aging_kidney@active.ident[aging_kidney$age%in%"old"]))
kidney_data<-as.data.frame(a@active.ident)
colnames(kidney_data) <- "ident"

ggplot(kidney_data, aes(x=as.factor(ident), fill=as.factor(ident) )) +  
  geom_bar( ) + theme_classic() +
  scale_fill_manual(values = c( "#A63603",     "#FDBE85",  "#E6550D",  "#FD8D3C"     ,"#A1D99B","#238B45",    "#4292C6" ,"#9ECAE1", "#FEB24C", "#08306B",  "#4292C6", "#2171B5"    ,"#54278F", "#9E9AC8"     ,"#FCCDE5", "#DE2D26", "#FB6A4A", "#FCAE91"     ,"#543005", "#8C510A", "#BF812D","#01665E")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) 

#-----------------------------------------------------------------------------------------------------------------------------
### dotplot
#-----------------------------------------------------------------------------------------------------------------------------

name = names(table(aging_kidney@active.ident))

table(aging_kidney@active.ident)
a = subset(aging_kidney, idents = name[1])
bim <- FindMarkers(a, ident.1 =  "old",  ident.2 = "young", only.pos = FALSE, verbose = T, group.by = "age", test.use = "wilcox", min.pct = 0)
bim = bim[order(bim$avg_log2FC, decreasing = T),]
bim$ident<-name[1]
bim2<-bim
head(bim2[order(bim3$avg_log2FC, decreasing = T),],20)
bim3 = bim2[rownames(bim2)%in%lncrna$gene_list,]
bim3$gene = rownames(bim3)
bim3


for ( i in c(2:9,11:15)) {
  a = subset(aging_kidney, idents = name[i])
  bim <- FindMarkers(a, ident.1 =  "old",  ident.2 = "young", only.pos = FALSE, verbose = T, group.by = "age", test.use = "wilcox", logfc.threshold = 0.15,  min.pct = 0)
  bim$ident<-name[i]
  bim2<-bim
  #  bim2<-bim2[bim2$p_val<0.05,]  
  bim2 = bim2[rownames(bim2)%in%lnc$gene_list,]
  bim2$gene = rownames(bim2)
  bim3 = rbind(bim3, bim2)
  
}

genes_to_plot = c(  "D730003I15Rik",	"Dubr",	"C130021I20Rik",	"Gm15866",	"Sp3os",	"2810433D01Rik",	"Gm26802",	"4631405J19Rik",	"Gm13861",	"Gm12840",	"Gm42679",	"Gm43847",	"Gm14964",	"Gm26771",	"Gm2415",	"1200007C13Rik",	"E030013I19Rik",	"Gm43112",	"Gm26652",	"Gm42397",	"BC049987",	"Gm15638",	"Gm14764",	"0610040F04Rik",	"Gm15563",	"Gm48281",	"3110045C21Rik",	"Gm12326",	"Gm48287",	"1110019B22Rik",	"Gm16310",	"Gm15848",	"AI838599",	"Hoxd3os1",	"Gm15441",	"Gm20485",	"Bvht",	"Meg3",	"Plet1os",	"3110099E03Rik",	"Gm13470",	"Gm37350",	"Gm15987",	"Tnfsf13os",	"Gm31243",	"Gm2682",	"1810041H14Rik",	"Gm16152",	"Gm30054",	"A930024E05Rik",	"AW112010",	"AI839979",	"Gm46224" )

DotPlot(aging_kidney, col.min = 0, features = unique(genes_to_plot), cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

exp_mat <- as.matrix(aging_kidney[["RNA"]]@data[genes_to_plot,])
aging_kidney_aging$ident = aging_kidney_aging@active.ident

meta <- aging_kidney_aging@meta.data %>% 
  select(ident)
meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
meta <- pivot_longer(meta, -ident, names_to="Gene", values_to="Expression")

meta_summary <- meta %>%
  group_by(ident, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100)

meta_summary$Gene = factor(meta_summary$Gene, levels = genes_to_plot)
meta_summary$fc = dot$avg_log2FC
meta_summary = as.data.frame(meta_summary)

lnc = c(  "D730003I15Rik",	"Dubr",	"C130021I20Rik",	"Gm15866",	"Sp3os",	"2810433D01Rik",	"Gm26802",	"4631405J19Rik",	"Gm13861",	"Gm12840",	"Gm42679",	"Gm43847",	"Gm14964",	"Gm26771",	"Gm2415",	"1200007C13Rik",	"E030013I19Rik",	"Gm43112",	"Gm26652",	"Gm42397",	"BC049987",	"Gm15638",	"Gm14764",	"0610040F04Rik",	"Gm15563",	"Gm48281",	"3110045C21Rik",	"Gm12326",	"Gm48287",	"1110019B22Rik",	"Gm16310",	"Gm15848",	"AI838599",	"Hoxd3os1",	"Gm15441",	"Gm20485",	"Bvht",	"Meg3",	"Plet1os",	"3110099E03Rik",	"Gm13470",	"Gm37350",	"Gm15987",	"Tnfsf13os",	"Gm31243",	"Gm2682",	"1810041H14Rik",	"Gm16152",	"Gm30054",	"A930024E05Rik",	"AW112010",	"AI839979",	"Gm46224" )

meta_summary = meta_summary[meta_summary$Gene%in%lnc , ]
meta_summary[meta_summary$fc < -0.5,]$fc = -0.5
meta_summary[meta_summary$fc > 0.5,]$fc = 0.5
    
ggplot(meta_summary, aes(x=Gene, y=ident)) +
  geom_point(aes(size = Pct, fill = fc), color="white", shape=21) +
  scale_size("% detected", range = c(0,11)) +
  scale_fill_gradientn(colours = c("blue", "grey", "Red"), limits = c(-0.5, 0.5)) +
    theme_classic() +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14)) 
                        
#-----------------------------------------------------------------------------------------------------------------------------------------------
### correlation
#-----------------------------------------------------------------------------------------------------------------------------------------------

kidney_subset<-subset(aging_kidney, idents = c("mesangial"))

matrix1<-kidney_subset@assays$RNA@counts
dim(matrix1)
matrix1<-as.matrix(matrix1)
matrix1<-t(matrix1)
matrix1<-as.data.frame(matrix1)
#options(warn=0)
matrix1<-cor(matrix1)


data = matrix1["Gm12840",]
data = data.frame(data)
dim(data)
data = na.omit(data)

data$gene = data$data
data = data[order(data$data, decreasing = T), ]
dim(data)
data = data[ 2:17994, ]
head(data)
data$gene = 1:17993
head(data)

ggplot(data, aes(x=gene, y=data)) +
  geom_line( color="navy", size=1, alpha=0.9, linetype=1) +
  theme_classic()


#--------------------------------------------------------------------------------------------------------
### senescence and inflammation gene set score
#--------------------------------------------------------------------------------------------------------

#GOBP senescence
senescece_gene <- list(c("Nsmce2",	"Hras",	"Prmt6",	"Hmga1b",	"Hmga1",	"Smc5",	"Nek4",	"Zmpste24",	"Zfp277",	"Nuak1",	"Arg2",	"Opa1",	"Pnpt1",	"Id2",	"Kras",	"Bmpr1a",	"Pla2r1",	"Mapkapk5",	"Trp53",	"Pawr",	"Smc6",	"Spi1",	"Arntl",	"Calr",	"B2m",	"Map2k1",	"Cav1",	"Vash1",	"Kat6a",	"Cgas",	"Abi3",	"Rsl1d1",	"Icmt",	"Eef1e1",	"Ypel3",	"Suv39h1",	"Mapk14",	"Srf",	"Ecrg4",	"Ing2",	"Prkcd",	"Cdkn2a",	"Cdkn2b",	"Pml",	"Wnt16",	"Cdkn1a"))

#hallmarke inflammatory
inflammatory_gene <- list(c("Abca1",	"Abi1",	"Acvr1b",	"Acvr2a",	"Adm",	"Adora2b",	"Adrm1",	"Ahr",	"Aplnr",	"Aqp9",	"Atp2a2",	"Atp2b1",	"Atp2c1",	"Axl",	"Bdkrb1",	"Best1",	"Bst2",	"Btg2",	"C3ar1",	"C5ar1",	"Calcrl",	"Ccl17",	"Ccl12",	"Ccl20",	"Ccl22",	"Ccl24",	"Ccl5",	"Ccl7",	"Ccr7",	"Ccrl2",	"Cd14",	"Cd40",	"Cd48",	"Cd55",	"Cd69",	"Cd70",	"Cd82",	"Cdkn1a",	"Chst2",	"Clec5a",	"Cmklr1",	"Csf1",	"Csf3",	"Csf3r",	"Cx3cl1",	"Cxcl10",	"Cxcl11",	"Cxcl5",	"Cxcl9",	"Cxcr6",	"Cybb",	"Dcbld2",	"Ebi3",	"Edn1",	"Eif2ak2",	"Emp3",	"Adgre1",	"Ereg",	"F3",	"Ffar2",	"Fpr1",	"Fzd5",	"Gabbr1",	"Gch1",	"Gna15",	"Gnai3",	"Gp1ba",	"Gpc3",	"Gpr132",	"Gpr183",	"Has2",	"Hbegf",	"Hif1a",	"Hpn",	"Hrh1",	"Icam1",	"Icam4",	"Icosl",	"Ifitm1",	"Ifnar1",	"Ifngr2",	"Il10",	"Il10ra",	"Il12b",	"Il15",	"Il15ra",	"Il18",	"Il18r1",	"Il18rap",	"Il1a",	"Il1b",	"Il1r1",	"Il2rb",	"Il4ra",	"Il6",	"Il7r",	"Inhba",	"Irak2",	"Irf1",	"Irf7",	"Itga5",	"Itgb3",	"Itgb8",	"Kcna3",	"Kcnj2",	"Kcnmb2",	"Kif1b",	"Klf6",	"Lamp3",	"Lck",	"Lcp2",	"Ldlr",	"Lif",	"Lpar1",	"Lta",	"Ly6e",	"Lyn",	"Marco",	"Mefv",	"Mep1a",	"Met",	"Mmp14",	"Msr1",	"Mxd1",	"Myc",	"Nampt",	"Ndp",	"Nfkb1",	"Nfkbia",	"Nlrp3",	"Nmi",	"Nmur1",	"Nod2",	"Npffr2",	"Olr1",	"Oprk1",	"Osm",	"Osmr",	"P2rx4",	"P2rx7",	"P2ry2",	"Pcdh7",	"Pde4b",	"Pdpn",	"Pik3r5",	"Plaur",	"Prok2",	"Psen1",	"Ptafr",	"Ptger2",	"Ptger4",	"Ptgir",	"Ptpre",	"Pvr",	"Raf1",	"Rasgrp1",	"Rela",	"Rgs1",	"Rgs16",	"Rhog",	"Ripk2",	"Rnf144b",	"Ros1",	"Rtp4",	"Scarf1",	"Scn1b",	"Sele",	"Sell",	"Selenos",	"Sema4d",	"Serpine1",	"Sgms2",	"Slamf1",	"Slc11a2",	"Slc1a2",	"Slc28a2",	"Slc31a1",	"Slc31a2",	"Slc4a4",	"Slc7a1",	"Slc7a2",	"Sphk1",	"Sri",	"Stab1",	"Tacr1",	"Tacr3",	"Tapbp",	"Timp1",	"Tlr1",	"Tlr2",	"Tlr3",	"Tnfaip6",	"Tnfrsf1b",	"Tnfrsf9",	"Tnfsf10",	"Tnfsf15",	"Tnfsf9",	"Tpbg",	"Vip"))


aging_kidney <- AddModuleScore(
  object = aging_kidney,
  features = senescece_gene,
  ctrl = 46,
  name = 'senescece_gene'
)

VlnPlot(a, features = c("inflammatory_gene1", "senescece_gene1"), group.by = "age")


#--------------------------------------------------------------------------------------------------------
### score gene exprssion
#--------------------------------------------------------------------------------------------------------

a = subset(aging_kidney, idents  = "mesangial")

table(a@active.ident)
summary(a$senescece_gene1)
summary(a$inflammatory_gene1)

a2 = a@meta.data
dim(a2)
a3 = a2[a2$senescece_gene1 > mean(a2$senescece_gene1),]
dim(a3)

Idents(a, cells = colnames(a)) = "low"
Idents(a, cells = rownames(a3)) = "high"
table(a@active.ident)

a@active.ident = factor(a@active.ident, levels = c("low", "high"))
VlnPlot(a, features = "Gm14964", cols = c("navy",  "red"), pt.size = 1) + stat_compare_means()

  
#--------------------------------------------------------------------------------------------------------
### beeswarm plot
#--------------------------------------------------------------------------------------------------------

ident_names<-names(table(aging_kidney@active.ident))
a<-subset(aging_kidney, idents = ident_names[1])
bim <- FindMarkers(a, ident.1 =  "young",  ident.2 = "old", only.pos = FALSE  , verbose = T, test.use = "wilcox", group.by = "age", logfc.threshold = 0.1)
bim$ident<-ident_names[1]
bim2<-bim
bim2<-bim2[bim2$p_val<0.05,]
#bim2<-bim2[abs(bim2$avg_logFC)>1,]
bim2<-bim2

for (i in 2:length(ident_names)) {
  
  a<-subset(aging_kidney, idents = ident_names[i])
  bim <- FindMarkers(a, ident.1 =  "young",  ident.2 = "old", only.pos = FALSE  , verbose = T, test.use = "wilcox", group.by = "age", logfc.threshold = 0.1)
  bim$ident<-ident_names[i]
  bim<-bim[bim$p_val<0.05,]
  bim2<-rbind(bim,bim2)
  
}

bim3<-bim2[rownames(bim2)%in%lnc$gene_list,]

library(beeswarm)
a = c("#A63603",     "#FDBE85",  "#E6550D",  "#FD8D3C"     ,"#A1D99B","#41AB5D","#238B45"      ,"#9ECAE1", "#FEB24C", "#08306B",  "#4292C6", "#2171B5"    ,"#54278F", "#9E9AC8"     ,"#FCCDE5", "#DE2D26", "#FB6A4A", "#FCAE91"     ,"#543005", "#8C510A", "#BF812D","#01665E")

beeswarm(avg_log2FC ~ ident, data=bim3, pch=19, method="swarm", col=a, cex=0.5, corral = c("wrap"))


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### volcano plot
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

name = names(table(aging_kidney@active.ident))

table(aging_kidney@active.ident)
a = subset(aging_kidney, idents = name[1])
bim <- FindMarkers(a, ident.1 =  "old",  ident.2 = "diet", only.pos = FALSE, verbose = T, group.by = "age", test.use = "wilcox", logfc.threshold = 0, min.pct = 0, min.diff.pct = 0)
bim$ident<-name[1]
bim2<-bim
bim3 = bim2[order(bim2$avg_log2FC, decreasing = F),]

bim3 = bim3[rownames(bim3)%in%lnc$gene_list,]
bim3$Gene = rownames(bim3)

for ( i in c(2:15, 17:19)) {
  a = subset(aging_kidney, idents = name[i])
  bim <- FindMarkers(a, ident.1 =  "old",  ident.2 = "diet", only.pos = FALSE, verbose = T, group.by = "age", test.use = "wilcox", logfc.threshold = 0,  min.pct = 0, min.diff.pct = 0)
  bim$ident<-name[i]
  bim2<-bim
  #  bim2<-bim2[bim2$p_val<0.05,]  
  bim2 = bim2[rownames(bim2)%in%lnc$gene_list,]
  bim2$Gene = rownames(bim2)
  bim3 = rbind(bim3, bim2)  
}

bim2 = bim

library(dplyr)
bim3 <- bim2 %>% 
  mutate(
    Expression = case_when(avg_log2FC >= 0.2 & p_val <= 0.05 ~ "Up-regulated",
                           avg_log2FC <= -0.2 & p_val <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

head(bim3) 
bim3$gene = NA
bim3$gene[ rownames(bim3) %in% rownames(bim2)] = rownames(bim2[ rownames(bim2) %in% rownames(bim3), ]) 


library(ggplot2)
ggplot(bim3, aes(avg_log2FC, -log(  p_val,10), colour=Expression)) +
  geom_point(size=1) +
  scale_color_manual(values=c("blue", "grey","red"))+
  # xlim(c(-4.5, 4.5)) +
  geom_vline(xintercept=c(-0.2,0.2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 2,lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="Differential expression")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())


