#--------------------------------------------------------------------------------------------------------
### lncRNA atlas analysis R script
#--------------------------------------------------------------------------------------------------------

setwd('/data/')
file_list <- list.files()
print(file_list)

sampleNames = mapply(function(x) x[length(x)], strsplit(file_list, split = '__'))
sampleNames <- gsub('-','_',sampleNames)
meta.info = data.table(file_Path = file_list, sampleID = sampleNames)

meta.info[, SubType:= mapply(function(x) x[1], strsplit(sampleID, '-'))]
DT::datatable(meta.info)

norm_method = 'standard'
treat = FALSE
show = FALSE
MT_cut = 100

batch_list = list()

for(i in 1:nrow(meta.info)){
  dir_of_10X = meta.info$file_Path[i]
  sampleID = meta.info$sampleID[i]
  subtype = meta.info$SubType[i]
  seq_batch = meta.info$Batch[i]
  print(paste("Starting processing", sampleID, 'at', Sys.time()))
  sc_mat = Read10X(dir_of_10X)
  
  lncrna <- CreateSeuratObject(counts = sc_mat, project = sampleID, min.cells = 1)
  lncrna
  
  lncrna[["percent.mt"]] <- PercentageFeatureSet(lncrna, pattern = "^mt-")
  
  qc = VlnPlot(lncrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  plot1 <- FeatureScatter(lncrna, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(lncrna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  
  lncrna <- subset(lncrna, percent.mt <= MT_cut & nFeature_RNA < 4000 )
  
  if(norm_method == 'SCT'){
    lncrna <- SCTransform(lncrna, vars.to.regress = 'percent.mt')
  } else{
    lncrna <- NormalizeData(lncrna, verbose = F)
    lncrna <- FindVariableFeatures(lncrna, selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)
    
    lncrna <- ScaleData(lncrna, verbose = F)
  }
  
  
  lncrna <- RenameCells(object = lncrna, add.cell.id = sampleID)  # add sample name as prefix
  if(treat == TRUE){
    lncrna <- RunPCA(lncrna, verbose = T)
    lncrna <- RunUMAP(lncrna, dims = 1:30, verbose = FALSE)
    lncrna <- RunTSNE(lncrna, dims = 1:30, verbose = FALSE)
    
    lncrna <- FindNeighbors(lncrna, dims = 1:30, verbose = FALSE)
    lncrna <- FindClusters(lncrna, verbose = FALSE)
    
    umap = DimPlot(lncrna, label = TRUE) + NoLegend()
    tsne = TSNEPlot(lncrna, label = TRUE) + NoLegend()
    
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
  lncrna@meta.data[,'SubType'] = subtype
  lncrna@meta.data[,'SampleID'] = sampleID
  lncrna@meta.data[,'Batch'] = seq_batch
  ### add batch data into a list to tissue
  #assign(paste0(sampleID,'.Seurat'),lncrna)
  batch_list <- append(batch_list, lncrna)
  print(paste("Finishing processing", sampleID, 'at', Sys.time()))
}

tissue <- Reduce(function(x,y){tissue(x,y,tissue.data = TRUE)}, batch_list)

VlnPlot(tissue, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

tissue <- subset(tissue, subset =  nFeature_RNA < 4000  & percent.mt <= 30 & nCount_RNA < 20000 )

VlnPlot(tissue, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

tissue <- FindVariableFeatures(tissue, selection.method = "vst", nfeatures = 2000, verbose = FALSE)



#--------------------------------------------------------------------------------------------------------
### identical cells isolation from pre-capture and post-capture
#--------------------------------------------------------------------------------------------------------

name = names(table(tissue$orig.ident))
name

for (i in 1:11) {
  
  a = subset(tissue, idents = name[1])
  head(colnames(a))
  a2 = subset(tissue, idents = name[2])
  head(colnames(a2))

  a3 = mapply(function(x) x[4], strsplit(colnames(a), "_"))
  head(a3)
  a4 = mapply(function(x) x[4], strsplit(colnames(a2), "_"))
  head(a4)

  print(table(a3%in%a4))

}

a = subset(tissue, idents = name[1])
head(colnames(a))
a2 = subset(tissue, idents = name[9])
head(colnames(a2))

a3 = mapply(function(x) x[4], strsplit(colnames(a), "_"))
head(a3)
a4 = mapply(function(x) x[4], strsplit(colnames(a2), "_"))
head(a4)

print(table(a3%in%a4))
a5 = a3[a3%in%a4]
a5 = as.data.frame(a5)
a7 = rbind(as.matrix(paste0(name[9], "_", a5$a5)), as.matrix(paste0(name[16], "_", a5$a5  )))


for (i in c(1,3,19,21,23)) {
  
  a = subset(tissue, idents = name[i])
  head(colnames(a))
  a2 = subset(tissue, idents = name[i+1])
  head(colnames(a2))
  
  a3 = mapply(function(x) x[4], strsplit(colnames(a), "_"))
  head(a3)
  a4 = mapply(function(x) x[4], strsplit(colnames(a2), "_"))
  head(a4)
  
  print(table(a3%in%a4))
  a5 = a3[a3%in%a4]
  a5 = as.data.frame(a5)
  a6 = rbind(as.matrix(paste0(name[i], "_", a5$a5)), as.matrix(paste0(name[i+1], "_", a5$a5  )))
  a7 = rbind(a6, a7)  
}

for (i in c(5,6,7,8,10,11)) {
  
  a = subset(tissue, idents = name[i])
  head(colnames(a))
  a2 = subset(tissue, idents = name[i+7])
  head(colnames(a2))
  
  a3 = mapply(function(x) x[5], strsplit(colnames(a), "_"))
  head(a3)
  a4 = mapply(function(x) x[5], strsplit(colnames(a2), "_"))
  head(a4)
  
  print(table(a3%in%a4))
  a5 = a3[a3%in%a4]
  a5 = as.data.frame(a5)
  a6 = rbind(as.matrix(paste0(name[i], "_", a5$a5)), as.matrix(paste0(name[i+7], "_", a5$a5  )))
  a7 = rbind(a6, a7)  
}

a7 = as.data.frame(a7)
head(a7)

tissue<-subset(tissue, cells = a7$V1)


#---------------------------------------------------------------------------------------------------------------------------------------------------------
### doublelet remove
#---------------------------------------------------------------------------------------------------------------------------------------------------------

orig.ident <- names(table(tissue$orig.ident))

library(DoubletFinder)

double<-list("a1",	"a2",	"a3",	"a4",	"a5",	"a6",	"a7",	"a8")
table(tissue$orig.ident)

# chagne the argument value according to the sample

a <- subset(tissue, cells = names(tissue@active.ident[tissue$orig.ident%in%orig[1]]))
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
a <- doubletFinder_v3(a, pN = 0.3, pK = 0.18, nExp = nExp, PCs = 1:30)
DF.name = colnames(a@meta.data)[grepl("DF.classification", colnames(a@meta.data))]



##############################################################################################
### Each tissue subclustering
##############################################################################################

##############################################################################################
### kidney
##############################################################################################

name = names(table(tissue@active.ident))
name

kidney_lnc<-subset(tissue, idents = name[5:11] )
kidney_mrna<-subset(tissue, idents = name[12:18] )

kidney_mat<-(kidney_mrna@assays$RNA@counts+kidney_lnc@assays$RNA@counts)

kidney <- CollapseSpeciesExpressionMatrix(kidney_mat)
kidney <- CreateSeuratObject(counts = kidney, min.cells = 3)
kidney[["percent.mt"]] <- PercentageFeatureSet(kidney, pattern = "^mt-")

VlnPlot(kidney, group.by = "orig.ident" , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

kidney <- NormalizeData(kidney, normalization.method = "LogNormalize", scale.factor = 10000)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 2000)
kidney <- ScaleData(kidney, vars.to.regress = c('nCount_RNA', 'percent.mt', 'nFeatures'), verbose = TRUE)
kidney <- RunPCA(kidney, npcs = 50, verbose = FALSE)

ElbowPlot(object = kidney, ndims = 50)

kidney <- FindNeighbors(kidney, dims = 1:30)
kidney <- RunUMAP(kidney, dims = 1:30)
kidney <- FindClusters(kidney, resolution = 2)


a = kidney_mrna@meta.data
a2 = kidney@meta.data

head(a)
head(a2)

a$cell = rownames(a)
a2$cell = rownames(a2)

a2 = tissue(a, a2, "cell")
head(a2)
rownames(a2) = a2$cell
a2 = a2[rownames(kidney@meta.data), ]
head(kidney@meta.data)
head(a2)

kidney@meta.data = a2


library(SeuratWrappers)

head(kidney)
table(kidney$orig.ident.x)
kidney.list <- SplitObject(kidney, split.by = "orig.ident.x")

for (i in 1:length(lung.list)) {
  kidney.list[[i]] <- NormalizeData(kidney.list[[i]],assay="RNA",normalization.method="LogNormalize")
  kidney.list[[i]] <- FindVariableFeatures(kidney.list[[i]], selection.method="vst", nfeatures=2000)
}

kidney <- RunFastMNN(object.list = kidney.list)

kidney5 <- RunUMAP(kidney, reduction = "mnn", dims = 1:30)
kidney5 <- FindNeighbors(kidney5, reduction = "mnn", dims = 1:30)
kidney5 <- FindClusters(kidney5, resolution = 0.5)

DimPlot(kidney5, label = TRUE, label.size = 5) 

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


DotPlot(kidney, col.min = 0, features = musGenes, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


kidney<- RenameIdents(object = kidney, "0"="PT1",	"1"="PT2",	"2"="GEC",	"3"="PT2",	"4"="PT1",	"5"="aterioles",	"6"="mesangial",	"7"="DCT",	"8"="ALOH",	"9"="Macro",	"10"="PC",	"11"="PEC",	"12"="GEC",	"13"="podocyte",	"14"="CD8",	"15"="PST",	"16"="IC",	"17"="JGC",	"18"="B cell",	"19"="RBC" )
ident = c( "podocyte",	"mesangial",	"JGC",	"PEC",	"PT1",	"PT2",	"PST",	"ALOH",	"DCT",	"PC",	"IC",	"GEC",	"aterioles",	"B cell",	"CD8",	"Macro",	"RBC" )
kidney@active.ident = factor( kidney@active.ident, levels = ident)


lnc = c("2810433D01Rik",	"Wt1os",	"C130021I20Rik",	"Gm15866",	"Gm13814",	"Dubr", "Gm26771",		"E030013I19Rik",	"A830082K12Rik", "Carmn",	"Gm12840",	"Gm13861",	"9530034E10Rik",	"Tbx3os1",	"4631405J19Rik",	"Gm42679",	"Dhx58os",	"Gm26802", "Gm19950",	"AI314278",	"D630029K05Rik",	"0610005C13Rik",	"Gm42937",	"Gm42397",	"C730027H18Rik",	"Gm49542",	"Pantr1", "Eprn",	"Gm17750",	"1110019B22Rik",	"Gm11789",	"Gm12927", "Gm47708",	"Hoxb5os", "Lhx1os",	"Gm15848",	"4930461G14Rik",	"D630024D03Rik",	"Hoxd3os1",	"AI838599",	"Gm20485",	"Gm16136",	"Gm13594",	"Gm53",	"2610016A17Rik",	"Gm17634",	"Gm12121",	"Bvht",	"3110099E03Rik",	"Gm35167",	"Gm31243",	"Mir142hg",	"Gm26740",	"Gm30211",	"AW112010",	"Gm15472",	"Gm19585",	"Gm2682",	"Gm43065",	"Gm34084",	"F630028O10Rik" )

DotPlot(kidney, col.min = 0, features = unique(lnc), cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


lnc = c("2810433D01Rik",	"Gm26771",	"Gm13861",	"Gm26802",  "Gm19950", "Gm26917",	 "Gm12927",  "Gm15848",  "AI838599",	"2610016A17Rik", "Bvht",	"3110099E03Rik",	"Gm31243",	"Gm15472",		"Gm34084" )

VlnPlot(kidney, features = unique(lnc), stack = T, flip = T ) + NoLegend()

#--------------------------------------------------------------------------------------------------------------------------------------
#correlation heatmap
#--------------------------------------------------------------------------------------------------------------------------------------

list = c("2810433D01Rik",	"Gm26771",	"Gm13861",	"Gm26802",  "Gm19950", "Gm26917",	 "Gm12927",  "Gm15848",  "AI838599",	"2610016A17Rik", "Bvht",	"3110099E03Rik",	"Gm31243",	"Gm15472",		"Gm34084")
VlnPlot(kidney, features = list , stack = T, flip = T) + NoLegend()

list<-list[list%in%rownames(kidney)]

matrix1<-kidney@assays$RNA@counts
matrix1<-matrix1[rownames(matrix1)%in%list,]
dim(matrix1)
matrix1<-as.matrix(matrix1)
matrix1<-t(matrix1)
matrix1<-as.data.frame(matrix1)
#options(warn=0)
matrix1<-cor(matrix1)

matrix2<-matrix1[list,]
matrix2<-matrix2[,list]
dim(matrix2)

rownames(matrix2)<-factor(rownames(matrix2), levels = list) 
colnames(matrix2)<-factor(colnames(matrix2), levels = list) 

library(corrplot)
corrplot(matrix2, method = "color", col = colorRampPalette(c("blue","white","red"))(200), tl.col = c("black"))



##############################################################################################
### lung
##############################################################################################

lung_mrna<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("lung_mrna_lung")]))
lung_lnc<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("lung_lnc_lung")]))

lung_mat<-(lung_mrna@assays$RNA@counts+lung_lnc@assays$RNA@counts)

lung   <- CreateSeuratObject(counts = lung_mat, min.cells = 3)
lung[["percent.mt"]] <- PercentageFeatureSet(lung, pattern = "^mt-")

VlnPlot(lung, group.by = "orig.ident" , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

lung   <- NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
lung   <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 2000)
lung   <- ScaleData(lung, vars.to.regress = c('nCount_RNA', 'percent.mt', 'nFeatures'), verbose = TRUE)
lung   <- RunPCA(lung, npcs = 50, verbose = FALSE)

ElbowPlot(object = lung, ndims = 50)

lung   <- FindNeighbors(lung, dims = 1:30)
lung   <- RunUMAP(lung, dims = 1:30)
lung   <- FindClusters(lung, resolution = 9)

DimPlot(lung, label = TRUE) + NoLegend()

gene = c("Ager",	"Egfr",	"Pdpn",	"Clic5",	"Sftpb",	"Sftpa1",	"Abca3",	"Sftpc",	"Ccdc113",	"Cyp2f2",	"Scgb3a2",	"Scgb1a1",	"Ptprb",	"Cd93",	"Tmem100",	"Dcn",	"Lum",	"Tgfbi",	"Col1a2",	"Upk3b",	" Wt1",	" Lrrn4",	" Msln",	"Mustn1",	"Acta2",	"Tagln",	"Notch3",	"Pdgfrb",	"Mki67",	"Top2a",	"Ube2c",	"Mmrn1",	"Prox1",	"Tbx1",	"Sema3d",	"Cd79a",	"Cd79b",	"Cd3d",	"Lef1",	"Cxcr6",	"Icos",	"Nkg7",	"Gzma",	"Lpl",	"Abcg1",	"Cd209a",	"Mgl2",	"Ccl3",	"Ccl9",	"Ccl4",	"S100a4",	"Gngt2",	"Lst1",	"S100a9",	"S100a8",	"Hba-a1",	"Hba-a2")

DotPlot(lung, col.min = 0, features = gene[gene%in%rownames(lung)], cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

lung   <- RenameIdents(object = lung, "0"="Naïve T",	"1"="B cell",	"2"="Fibroblast1",	"3"="NK",	"4"="Fibroblast1",	"5"="Gcap",	"6"="AT2",	"7"="AM",	"8"="Monocyte2",	"9"="Monocyte3",	"10"="Fibroblast2",	"11"="Gcap",	"12"="DC",	"13"="proliferating",	"14"="B cell2",	"15"="Gcap",	"16"="Gcap",	"17"="RBC",	"18"="Pericyte",	"19"="AT1",	"20"="Fibroblast3",	"21"="CD8",	"22"="Myofibro",	"23"="Neutrophils",	"24"="Gcap",	"25"="Basophils",	"26"="Mesothelial",	"27"="SMC2",	"28"="Immune stromal cell",	"29"="Ciliated",	"30"="Myeloid" )

ident = c("AT1",	"AT2",	"Ciliated",	"Gcap",	"Fibroblast1",	"Fibroblast2",	"Fibroblast3",	"Myofibro",	"Mesothelial",	"SMC2",	"Immune stromal cell",	"Pericyte",	"Naïve T",	"CD8",	"NK",	"B cell",	"B cell2",	"AM",	"Basophils",	"DC",	"Monocyte2",	"Monocyte3",	"Neutrophils",	"proliferating",	"Myeloid",	"RBC")
lung@active.ident = factor(lung@active.ident, levels = ident)


gene = c( "2410004I01Rik",	"Gm16212",	"Gm42702",	"2810433D01Rik",	"Gm42696",	"Gm43560",	"Gm42701",	"Gm43259",	"Platr22",	"Gm42699",	"Gm43707",	"Gm19990",	"Gm44751",	"9230117E06Rik",	"Gm28856",	"A530020G20Rik",	"Gm16136",	"3300002A11Rik",	"Gm29538",	"C230072F16Rik",	"D630024D03Rik",	"C330002G04Rik",	"D430036J16Rik",	"Smkr-ps",	"Gm16006",	"AU040972",	"Gm13580",	"1810008I18Rik",	"Gm44000",	"4933406B17Rik",	"Uckl1os",	"A930001A20Rik",	"Med9os",	"1810041H14Rik",	"Gm15441",	"D830026I12Rik",	"Gm15889",	"AU022754",	"Gm12984",	"Gm32950",	"Rptoros",	"Gm15587",	"1700095J03Rik",	"Gm20755",	"9130015G15Rik",	"Gm15728",	"9930004E17Rik",	"A230057D06Rik",	"Gm30025",	"Bvht",	"Gm9917",	"Gm26632",	"1200007C13Rik",	"Meg3",	"6330403K07Rik",	"Tnfsf13os",	"Gm26873",	"2310015D24Rik",	"4833412C05Rik",	"Hoxd3os1",	"4833422C13Rik",	"Gm2682",	"Gm19585",	"Gm33104",	"Gm31243",	"Gm30211",	"Bach2it1",	"AI504432",	"Gm41307",	"Gm45237",	"Gm15657",	"Gm36723",	"Gm16120",	"Gm36161",	"AI839979",	"Mirt2",	"Ptgs2os2",	"Gm43445" )

DotPlot(lung, col.min = 0, features = gene, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#--------------------------------------------------------------------------------------------------------------------------------------
#correlation heatmap
#--------------------------------------------------------------------------------------------------------------------------------------

list = c( "2410004I01Rik",	"Gm44751",	"C330002G04Rik",	"Bvht",	"1200007C13Rik", "Meg3",	"6330403K07Rik",	"Gm26873",	"4833412C05Rik",	"Gm2682",	"Gm31243",	"Gm41307",	"Gm45237",	"Gm36723",	"Gm36161",	"AI839979",	"Mirt2" )
list<-list[list%in%rownames(lung)]
VlnPlot(lung, features = list , stack = T, flip = T) + NoLegend()

list = c( "2410004I01Rik",	"Gm44751",	"C330002G04Rik",	"Bvht",	"1200007C13Rik",	"Meg3",	"6330403K07Rik",	"Gm26873",	"4833412C05Rik",	"Gm2682",	"Gm31243",	"Gm41307",	"Gm45237",	"Gm36723",	"Gm36161",	"AI839979",	"Mirt2",	"Ager", "Sftpb", "Ccdc113", "Ptprb", "Col1a1", "Col1a2", "Tgfbi", "Upk3b", "Mustn1", "Cd3d", "Cd79a", "Lpl", "Ccl3", "Cd209a", "S100a4", "Gngt2", "S100a9" )

matrix1<-lung@assays$RNA@counts
matrix1<-matrix1[rownames(matrix1)%in%list,]
dim(matrix1)
matrix1<-as.matrix(matrix1)
matrix1<-t(matrix1)
matrix1<-as.data.frame(matrix1)
#options(warn=0)
matrix1<-cor(matrix1)

matrix2<-matrix1[list,]
matrix2<-matrix2[,list]
dim(matrix2)

rownames(matrix2)<-factor(rownames(matrix2), levels = list) 
colnames(matrix2)<-factor(colnames(matrix2), levels = list) 

library(corrplot)
corrplot(matrix2, method = "color", col = colorRampPalette(c("blue","white","red"))(200), tl.col = c("black"))



##############################################################################################
### liver
##############################################################################################

liver_mrna<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("liver_mrna_liver")]))
liver_lnc<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("liver_lnc_liver")]))

liver_mat<-(liver_mrna@assays$RNA@counts+liver_lnc@assays$RNA@counts)

liver  <- CreateSeuratObject(counts = liver_mat, min.cells = 3)
liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^mt-")

VlnPlot(liver, group.by = "orig.ident" , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

liver  <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = 10000)
liver  <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 2000)
liver  <- ScaleData(liver, vars.to.regress = c('nCount_RNA', 'percent.mt', 'nFeatures'), verbose = TRUE)
liver  <- RunPCA(liver, npcs = 50, verbose = FALSE)

ElbowPlot(object = liver, ndims = 50)

liver  <- FindNeighbors(liver, dims = 1:21)
liver  <- RunUMAP(liver, dims = 1:21)
liver  <- FindClusters(liver, resolution = 2)

DimPlot(liver, label = TRUE)

gene = c( "Cd8a", "Cd4", "Sell", "Ccr7", "Lef1", "Cyp2e1",	"Apob",	"Pck1",	"Dcn",	"Cygb",	"Kdr",	"Eng",	"Gpihbp1",	"Pecam1",	"Cd79b",	"Cd79a",	"Trbc2",	"Cd3d",	"Thy1",	"Gzma",	"Nkg7",	"Top2a",	"Mki67",	"C1qa",	"Csf1r",	"Clec4f",	"Adgre1",	"Fut4",	"Itgam",	"Cd24a",	"Ly6g",	"S100a8",	"S100a9",	"Ly6c2",	"Clec10a",	"Fcgr1",	"Ccr2",	"Spn",	"Nr4a1",	"Hbb-bt",	"Hba-a2" )

DotPlot(liver, col.min = 0, features = gene[gene%in%rownames(liver)], cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

liver  <- RenameIdents(object = liver, "0"="B cell",	"1"="Endo1",	"2"="Kupffer",	"3"="T cell",	"4"="T cell2",	"5"="Mono",	"6"="Endo2",	"7"="Mono2",	"8"="NK",	"9"="Endo3",	"10"="RBC",	"11"="DC",	"12"="NK/T",	"13"="pat.Mono",	"14"="Endo4",	"15"="Neutro",	"16"="Cycling",	"17"="hepatocyte1",	"18"="T cell3",	"19"="Endo5",	"20"="hepatocyte2",	"21"="Stellate cell",	"22"="Granulocyte",	"23"="Granulocyte" )

ident = c( "hepatocyte1",	"hepatocyte2",	"Stellate cell",	"Endo1",	"Endo2",	"Endo3",	"Endo4",	"Endo5",	"B cell",	"T cell",	"T cell2",	"T cell3",	"NK/T",	"NK",	"Cycling",	"Kupffer",	"Mono",	"Mono2",	"pat.Mono",	"Neutro",	"Granulocyte",	"DC",	"RBC" )
liver@active.ident = factor(liver@active.ident, levels = ident)

DimPlot(liver, label = TRUE)


gene = c( "0610005C13Rik",	"1810008I18Rik",	"4732465J04Rik",	"Gm12909",	"B630019A10Rik",	"C730027H18Rik",	"Egfros",	"Gm43189",	"Gm36041",	"Gm30262",	"Gm49417",	"0610031O16Rik",	"Gm32063",	"9030616G12Rik",	"Gm31583",	"4930404H11Rik",	"AI480526",	"Gm40770",	"Gm42067",	"2210039B01Rik",	"Gm45792",	"Gm48302",	"Gm17491",	"Gm15998",	"Gm43258",	"E430024I08Rik",	"AI463229",	"C920006O11Rik",	"Gm48086",	"Gm17196",	"Hnf4aos",	"Gm34667",	"9030622O22Rik",	"Meg3",	"Gm42679",	"Carmn",	"Snhg18",	"Gm15222",	"Airn",	"4833422C13Rik",	"Gm14964",	"5033428I22Rik",	"9530082P21Rik",	"Dubr",	"Plet1os",	"Bvht",	"Gm16014",	"Gm31243",	"Gm30211",	"A630023P12Rik",	"Gm2682",	 "Gm20069",	"Gm19585",	"Gm28875",	"Lockd",	"Gm43914",	"4930455G09Rik",	"AI839979",	"Gm36161",	"Gm13373",	"Gm31814",	"1600010M07Rik",	"Mirt2",	"Ino80dos",	"Gm5547",	"Gm45715",	"Arhgap27os2" )

DotPlot(liver, col.min = 0, features = gene, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


#--------------------------------------------------------------------------------------------------------------------------------------
#correlation heatmap
#--------------------------------------------------------------------------------------------------------------------------------------

liver2 =  liver
liver2 <- RenameIdents(object = liver2, "Endo2" = "Endo1", "Endo3" = "Endo1", "Endo4" = "Endo1",  "Endo5" = "Endo1", "T cell2" = "T cell", "T cell3" = "T cell", "Mono2" = "Mono" )
table(liver2@active.ident)

ident = c( "hepatocyte1",	"hepatocyte2",	"Stellate cell",	"Endo1",	"B cell",	"T cell",		"NK/T",	"NK",	"Cycling",	"Kupffer",	"Mono",	"pat.Mono",	"Neutro",	"Granulocyte",	"DC",	"RBC" )
liver2@active.ident = factor(liver2@active.ident, levels = ident)

list = c(  "0610005C13Rik",	"Egfros",	"Gm42679",	"Plet1os",	"Gm31243",	"Gm2682",	"Gm28875",	"Lockd",	"Gm43914",	"AI839979",	"Gm31814", "1600010M07Rik",	"Gm45715" )
VlnPlot(a, features = list , stack = T, flip = T) + NoLegend()

list = c( "0610005C13Rik",	"Egfros",	"Gm42679",	"Plet1os",	"Gm31243",	"Gm2682",	"Gm28875",	"Lockd",	"Gm43914",	"AI839979",	"Gm31814", "1600010M07Rik",	"Gm45715",	"Cyp2e1",	"Pck1",	"Dcn",	"Kdr",	"Cd79a",	"Cd3d",	"Klrc1",	"Mki67",	"C1qa",	"Itgam", "Nr4a1",	"S100a8",	"Ly6c2" )
list<-list[list%in%rownames(a)]

matrix1<-liver2@assays$RNA@counts
matrix1<-matrix1[rownames(matrix1)%in%list,]
dim(matrix1)
matrix1<-as.matrix(matrix1)
matrix1<-t(matrix1)
matrix1<-as.data.frame(matrix1)
#options(warn=0)
matrix1<-cor(matrix1)

matrix2<-matrix1[list,]
matrix2<-matrix2[,list]
dim(matrix2)

rownames(matrix2)<-factor(rownames(matrix2), levels = list) 
colnames(matrix2)<-factor(colnames(matrix2), levels = list) 

library(corrplot)
corrplot(matrix2, method = "color", col = colorRampPalette(c("blue","white","red"))(200), tl.col = c("black"))



##############################################################################################
### intestine
##############################################################################################


intestine_mrna<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("intestine_mrna_intestine")]))
intestine_lnc<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("intestine_lnc_intestine")]))

intestine_mat<-(intestine_mrna@assays$RNA@counts+intestine_lnc@assays$RNA@counts)

intestine     <- CreateSeuratObject(counts = intestine_mat, min.cells = 3)
intestine[["percent.mt"]] <- PercentageFeatureSet(intestine, pattern = "^mt-")

VlnPlot(intestine, group.by = "orig.ident" , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

intestine     <- NormalizeData(intestine, normalization.method = "LogNormalize", scale.factor = 10000)
intestine     <- FindVariableFeatures(intestine, selection.method = "vst", nfeatures = 2000)
intestine     <- ScaleData(intestine, vars.to.regress = c('nCount_RNA', 'percent.mt', 'nFeatures'), verbose = TRUE)
intestine     <- RunPCA(intestine, npcs = 50, verbose = FALSE)

ElbowPlot(object = intestine, ndims = 50)

intestine     <- FindNeighbors(intestine, dims = 1:28)
intestine     <- RunUMAP(intestine, dims = 1:28)
intestine     <- FindClusters(intestine, resolution = 2)

DimPlot(intestine, label = TRUE)

gene = c( "Lgr5",	"Lrig1",	"Ascl2",	"Olfm4",	"Smoc2",	"Mki67",	"Ube2c",	"Col1a1",	"Wnt5a",	"Bmp2",	"Foxl1",	"Dcn",	"Vim",	"Col1a2",	"Pdgfrb",	"Notch3",	"Acta2",	"Mustn1",	"Myl9",	"Kdr",	"Pecam1",	"Lyve1",	"Prox1",	"Cd79a",	"Cd79b",	"Sdc1",	"Tnfrsf17",	"Jchain",	"Cd3d",	"Cd3e",	"Gzmk",	"Nkg7",	"Il7r",	"Rora",	"C1qa",	"C1qb",	"Siglech",	"Cd300c", "Cd74", "Cd68" )

DotPlot(intestine, col.min = 0, features = gene[gene%in%rownames(intestine)], cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

intestine     <- RenameIdents(object = intestine, "0"="Endo1",	"1"="B cell",	"2"="Myofibro1",	"3"="Cd4 T",	"4"="Cd8 T",	"5"="Endo2",	"6"="plasma",	"7"="Proli_ISC",	"8"="Endo3",	"9"="Fibro",	"10"="Naïve_Cd8",	"11"="proli_B cell",	"12"="VTT",	"13"="Endo4",	"14"="NK",	"15"="Naïve_Cd4",	"16"="Myofibro2",	"17"="Proli",	"18"="Myeloid",	"19"="Macro",	"20"="plasma2",	"21"="VTT2",	"22"="DC",	"23"="Pdgfrb+Fibro",	"24"="LV" )

ident = c( "Proli_ISC",	"Proli",	"VTT",	"VTT2",	"Endo1",	"Endo2",	"Endo3",	"Endo4",	"LV",	"Pdgfrb+Fibro",	"Fibro",	"Myofibro1",	"Myofibro2",	"proli_B cell",	"B cell",	"plasma",	"plasma2",	"Naïve_Cd4",	"Cd4 T",	"Naïve_Cd8",	"Cd8 T",	"NK",	"Myeloid", "Macro",	"DC" )
intestine@active.ident = factor(intestine@active.ident, levels = ident)

gene = c( "Lockd",	"Gm26682",	"Gm26903",	"Gm15222",	"Gm2670",	"Gm28875",	"Tbx3os1",	"Bvht",	"Gm12002",	"4833403J16Rik",	"Gm32688",	"9530082P21Rik",	"Hoxd3os1",	"Gm38910",	"Gm15902",	"1810008I18Rik",	"D5Ertd615e",	"D830026I12Rik",	"Gm11837",	"4833422C13Rik",	"Gm14964",	"Gm13861",	"B430010I23Rik",	"Adap2os",	"Gm13470",	"Gm26802",	"D630033O11Rik",	"Dnm3os",	"Gm7932",	"Meg3",	"Gm5084",	"4833412C05Rik",	"6330403K07Rik",	"Gm39459",	"2610035D17Rik",	"Gm2447",	"Gm31814",	"2010007H06Rik",	"3300005D01Rik",	"Gm15991",	"Gm4610",	"4930426D05Rik",	"5830418P13Rik",	"Gm31243",	"Gm19951",	"Gm34215",	"Kcnmb4os2",	"4833407H14Rik",	"Gm15472",	"A630023P12Rik",	"Gm20069",	"Gm29243",	"Gm2682",	"Bach2it1",	"1810041H14Rik",	"Gm14029",	"Gm36723",	"Gm43445",	"Slc36a3os",	"Gm34084",	"2900052N01Rik",	"Gm13391",	"Gm43914",	"Arhgap27os2", "Gm31718" )

DotPlot(intestine, col.min = 0, features = gene, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))



#--------------------------------------------------------------------------------------------------------------------------------------
#correlation heatmap
#--------------------------------------------------------------------------------------------------------------------------------------

intestine2 =  intestine

intestine2 <- RenameIdents(object = intestine2, "VTT2" = "VTT", "Endo2" = "Endo1", "Endo3" = "Endo1", "Endo4" = "Endo1",  "Myofibro2" = "Myofibro1", "plasma2" = "plasma" )
table(intestine2@active.ident)

ident = c( 	"Proli",	"VTT",	"Endo1",	"LV",	"Pdgfrb+Fibro",	"Fibro",	"Myofibro1",		"proli_B cell",	"B cell",	"plasma",	"Naïve_Cd4",	"Cd4 T",	"Naïve_Cd8",	"Cd8 T",	"NK",	"Myeloid", "Macro",	"DC" )
intestine2@active.ident = factor(intestine2@active.ident, levels = ident)
ident = names(table(intestine2@active.ident))

list = c(  "Lockd",	"Gm26903",	"Bvht",	"Gm15902",	"Gm13861",	"D630033O11Rik",	"4833412C05Rik",	"2610035D17Rik",	"Gm31243",	"Gm34215",	"A630023P12Rik",	"Gm20069",	"Gm2682",	"1810041H14Rik",	"Gm36723",	"Slc36a3os",	"2900052N01Rik",	"Gm43914" )
VlnPlot(a, features = list , stack = T, flip = T) + NoLegend()

list = c( "Ube2c",	"Wnt5a",	"Kdr",	"Lyve1",	"Pdgfrb",	"Dcn",	"Acta2",	"Top2a",	"Cd79a",	"Sdc1",	"Lef1", "Cd4",	"Sell",	"Cd8a",	"Nkg7",	"Cd74",	"C1qa",	"Siglech", "Lockd",	"Gm26903",	"Bvht",	"Gm15902",	"Gm13861",	"D630033O11Rik",	"4833412C05Rik",	"2610035D17Rik",	"Gm31243",	"Gm34215",	"A630023P12Rik", 	"Gm20069",	"Gm2682",	"1810041H14Rik",	"Gm36723", "Slc36a3os",	"2900052N01Rik",	"Gm43914"   )
list<-list[list%in%rownames(a)]
length(list)
list[1:5]


matrix1<-intestine2@assays$RNA@counts
matrix1<-matrix1[rownames(matrix1)%in%list,]
dim(matrix1)
matrix1<-as.matrix(matrix1)
matrix1<-t(matrix1)
matrix1<-as.data.frame(matrix1)
#options(warn=0)
matrix1<-cor(matrix1)

matrix2<-matrix1[list,]
matrix2<-matrix2[,list]
dim(matrix2)

rownames(matrix2)<-factor(rownames(matrix2), levels = list) 
colnames(matrix2)<-factor(colnames(matrix2), levels = list) 

library(corrplot)
corrplot(matrix2, method = "color", col = colorRampPalette(c("blue","white","red"))(200), tl.col = c("black"))



##############################################################################################
### heart
##############################################################################################


heart_mrna<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("heart_mrna_heart")]))
heart_lnc<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("heart_lnc_heart")]))

heart_mat<-(heart_mrna@assays$RNA@counts+heart_lnc@assays$RNA@counts)

heart   <- CreateSeuratObject(counts = heart_mat, min.cells = 3)
heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "^mt-")
VlnPlot(heart, group.by = "orig.ident" , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

heart   <- NormalizeData(heart, normalization.method = "LogNormalize", scale.factor = 10000)
heart   <- FindVariableFeatures(heart, selection.method = "vst", nfeatures = 2000)
heart   <- ScaleData(heart, vars.to.regress = c('nCount_RNA', 'percent.mt', 'nFeatures'), verbose = TRUE)
heart   <- RunPCA(heart, npcs = 50, verbose = FALSE)

ElbowPlot(object = heart, ndims = 50)

heart   <- FindNeighbors(heart, dims = 1:32)
heart   <- RunUMAP(heart, dims = 1:32)
heart   <- FindClusters(heart, resolution = 2)

DimPlot(heart, label = TRUE) + NoLegend()

musGenes <- c( "Myh6",	"Tnnt2",	"Tnni3",	"Kdr",	"Emcn",	"Pecam1",	"Lyve1",	"Cldn5",	"Dcn",	"Col1a1",	"Col1a2",	"Tagln",	"Mustn1",	"Myh11",	"Kcnj8",	"Vtn",	"Colec11",	"Cd79a",	"Cd79b",	"Nkg7",	"Gzma",	"Adgre1",	"Cd14",	"Napsa",	"Lgals3",	"Plbd1",	"S100a8",	"S100a9",	"Gp9",	"Itgb3",	"Itga2b",	"Hbb-bt",	"Hba-a1" )
DotPlot(heart, col.min = 0, features = musGenes, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

heart   <- RenameIdents(object = heart, "0"="Fibro1",	"1"="Fibro1",	"2"="Fibro3",	"3"="RM",	"4"="Fibro1",	"5"="Endo",	"6"="Fibro3",	"7"="RM2",	"8"="SMC",	"9"="B cell",	"10"="Fibro1",	"11"="Endo2",	"12"="LE",	"13"="Myofibro",	"14"="T cell",	"15"="Pericyte",	"16"="Endo3",	"17"="DC",	"18"="Granulocyte",	"19"="RBC",	"20"="Platelet",	"21"="Fibro7",	"22"="Cardiomyocyte" )

ident = c( "Proli_ISC",	"Proli",	"VTT",	"VTT2",	"Endo1",	"Endo2",	"Endo3",	"Endo4",	"LV",	"Pdgfrb+Fibro",	"Fibro",	"Myofibro1",	"Myofibro2",	"proli_B cell",	"B cell",	"plasma",	"plasma2",	"Naïve_Cd4",	"Cd4 T",	"Naïve_Cd8",	"Cd8 T",	"NK",	"Myeloid", "Macro",	"DC" )
heart@active.ident = factor(heart@active.ident, levels = ident)

DimPlot(heart, label = TRUE)

gene = c( "Lockd",	"Gm26682",	"Gm26903",	"Gm15222",	"Gm2670",	"Gm28875",	"Tbx3os1",	"Bvht",	"Gm12002",	"4833403J16Rik",	"Gm32688",	"9530082P21Rik",	"Hoxd3os1",	"Gm38910",	"Gm15902",	"1810008I18Rik",	"D5Ertd615e",	"D830026I12Rik",	"Gm11837",	"4833422C13Rik",	"Gm14964",	"Gm13861",	"B430010I23Rik",	"Adap2os",	"Gm13470",	"Gm26802",	"D630033O11Rik",	"Dnm3os",	"Gm7932",	"Meg3",	"Gm5084",	"4833412C05Rik",	"6330403K07Rik",	"Gm39459",	"2610035D17Rik",	"Gm2447",	"Gm31814",	"2010007H06Rik",	"3300005D01Rik",	"Gm15991",	"Gm4610",	"4930426D05Rik",	"5830418P13Rik",	"Gm31243",	"Gm19951",	"Gm34215",	"Kcnmb4os2",	"4833407H14Rik",	"Gm15472",	"A630023P12Rik",	"Gm20069",	"Gm29243",	"Gm2682",	"Bach2it1",	"1810041H14Rik",	"Gm14029",	"Gm36723",	"Gm43445",	"Slc36a3os",	"Gm34084",	"2900052N01Rik",	"Gm13391",	"Gm43914",	"Arhgap27os2", "Gm31718" )

DotPlot(heart, col.min = 0, features = gene, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))



#--------------------------------------------------------------------------------------------------------------------------------------
#correlation heatmap
#--------------------------------------------------------------------------------------------------------------------------------------

heart2 =  heart
heart2 <- RenameIdents(object = heart2, "VTT2" = "VTT", "Endo2" = "Endo1", "Endo3" = "Endo1", "Endo4" = "Endo1",  "Myofibro2" = "Myofibro1", "plasma2" = "plasma" )
table(heart2@active.ident)

ident = c( "Proli",	"VTT",	"Endo1",	"LV",	"Pdgfrb+Fibro",	"Fibro",	"Myofibro1", "proli_B cell",	"B cell",	"plasma",	"Naïve_Cd4",	"Cd4 T",	"Naïve_Cd8",	"Cd8 T",	"NK",	"Myeloid", "Macro",	"DC" )
heart2@active.ident = factor(heart2@active.ident, levels = ident)

list = c(  "Lockd",	"Gm26903",	"Bvht",	"Gm15902",	"Gm13861",	"D630033O11Rik",	"4833412C05Rik",	"2610035D17Rik",	"Gm31243",	"Gm34215",	"A630023P12Rik",	"Gm20069",	"Gm2682",	"1810041H14Rik",	"Gm36723",	"Slc36a3os",	"2900052N01Rik",	"Gm43914" )
VlnPlot(a, features = list , stack = T, flip = T) + NoLegend()

list = c( "Ube2c",	"Wnt5a",	"Kdr",	"Lyve1",	"Pdgfrb",	"Dcn",	"Acta2",	"Top2a",	"Cd79a",	"Sdc1",	"Lef1", "Cd4",	"Sell",	"Cd8a",	"Nkg7",	"Cd74",	"C1qa",	"Siglech", "Lockd",	"Gm26903",	"Bvht",	"Gm15902",	"Gm13861",	"D630033O11Rik",	"4833412C05Rik",	"2610035D17Rik",	"Gm31243",	"Gm34215",	"A630023P12Rik", 	"Gm20069",	"Gm2682",	"1810041H14Rik",	"Gm36723", "Slc36a3os",	"2900052N01Rik",	"Gm43914"   )
list<-list[list%in%rownames(a)]
length(list)
list[1:5]

matrix1<-heart2@assays$RNA@counts
matrix1<-matrix1[rownames(matrix1)%in%list,]
dim(matrix1)
matrix1<-as.matrix(matrix1)
matrix1<-t(matrix1)
matrix1<-as.data.frame(matrix1)
#options(warn=0)
matrix1<-cor(matrix1)

matrix2<-matrix1[list,]
matrix2<-matrix2[,list]
dim(matrix2)

rownames(matrix2)<-factor(rownames(matrix2), levels = list) 
colnames(matrix2)<-factor(colnames(matrix2), levels = list) 

library(corrplot)
corrplot(matrix2, method = "color", col = colorRampPalette(c("blue","white","red"))(200), tl.col = c("black"))



##############################################################################################
### thymus
##############################################################################################

thymus_mrna<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("thymus_mrna_thymus")]))
thymus_lnc<-subset(tissue, cells=names(tissue@active.ident[tissue$orig.ident %in% c("thymus_lnc_thymus")]))

thymus_mat<-(thymus_mrna@assays$RNA@counts+thymus_lnc@assays$RNA@counts)

thymus    <- CreateSeuratObject(counts = thymus_mat, min.cells = 3)
thymus[["percent.mt"]] <- PercentageFeatureSet(thymus, pattern = "^mt-")

VlnPlot(thymus, group.by = "orig.ident" , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

thymus    <- NormalizeData(thymus, normalization.method = "LogNormalize", scale.factor = 10000)
thymus    <- FindVariableFeatures(thymus, selection.method = "vst", nfeatures = 2000)
thymus    <- ScaleData(thymus, vars.to.regress = c('nCount_RNA', 'percent.mt', 'nFeatures'), verbose = TRUE)
thymus    <- RunPCA(thymus, npcs = 50, verbose = FALSE)

ElbowPlot(object = thymus, ndims = 50)

thymus    <- FindNeighbors(thymus, dims = 1:21)
thymus    <- RunUMAP(thymus, dims = 1:21)
thymus    <- FindClusters(thymus, resolution = 1.7)

DimPlot(thymus, label = TRUE) + NoLegend()

musGenes <- c( "Cd4",	"Cd8a",	"Rorc",	"Cd40lg",	"Sell",	"S1pr1",	"Itm2a",	"Id2",	"Ccr9",	"Tox",	"Ikzf2",	"Id3",	"Cldn10",	"Nkg7",	"Klrd1",	"Pcna",	"Mki67",	"Trdc",	"Runx3",	"Birc5",	"Hes1",	"Trgc2",	"Il2ra",	"Adgrg1",	"Notch1",	"Runx1",	"Cd8",	"H2-Aa",	"Napsa",	"Cd74",	"Tyrobp",	"Trgc4",	"Trgc1" )

DotPlot(thymus, col.min = 0, features = musGenes, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

thymus    <- RenameIdents(object = thymus, "0"="DP",	"1"="DP",	"2"="DP",	"3"="Cd4 SP(M)",	"4"="DP",	"5"="Cd8 SP",	"6"="Cd4 SP",	"7"="αβ T",	"8"="CD8 αα",	"9"="NK",	"10"="DP",	"11"="DP(P)",	"12"="DP",	"13"="DP2(P)",	"14"="DN",	"15"="ISP",	"16"="DC",	"17"="DP" )
ident = c( "DN",	"ISP",	"DP(P)",	"DP2(P)",	"DP",	"αβ T",	"CD8 αα",	"Cd8 SP",	"Cd4 SP",	"Cd4 SP(M)",	"NK", "DC" )
thymus@active.ident = factor(thymus@active.ident, levels = ident)

DimPlot(thymus, label = TRUE)

gene = c("Gm26930",	"Gm37065",	"Gm12840",	"A630023P12Rik",	"Gm44175",	"Lockd",	"H19",	"5830468F06Rik",	"Gm42047",	"5830418P13Rik",	"Gm12305",	"Brip1os",	"Dgkeos",	"2610020C07Rik",	"B630019A10Rik",	"Gm15340",	"AW112010",	"Gm38405",	"Gm37552")

DotPlot(thymus, col.min = 0, features = gene, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))



#--------------------------------------------------------------------------------------------------------------------------------------
# thymus trajectory heatmap
#--------------------------------------------------------------------------------------------------------------------------------------

table(thymus@active.ident)

thymus_subset<-subset(thymus, cells=names(thymus@active.ident[thymus@active.ident %in% c( "DN", "ISP",	"DP(P)",	"DP2(P)",	 "DP", "CD8 αα", "Cd8 SP", "Cd4 SP",	"Cd4 SP(M)" )]))

thymus_subset$nCount_RNA=NULL

thymus_subset$nFeature_RNA=NULL

thymus_subset$percent.mt=NULL

thymus_subset$orig.ident<-thymus_subset@active.ident

thymus_subset@meta.data=thymus_subset@meta.data[,1:3]

x=GetAssayData(object = thymus_subset, slot = "counts")[rownames(GetAssayData(object = thymus_subset, slot = "counts")) %in% rownames(GetAssayData(object = thymus_subset)), colnames(GetAssayData(object = thymus_subset, slot = "counts")) %in% colnames(GetAssayData(object = thymus_subset))]

target=x[,colnames(x) %in% rownames(thymus_subset@meta.data)]

target<-as.matrix(target)

target=as.matrix(t(target))

name=tissue(thymus_subset@meta.data, target, by="row.names")

rownames(name)=name$Row.names

name=name[,-1]

name[1:5,1:4]

thymus_subset@meta.data=name[,1:4]

name=name[,-1:-3]

target=t(name)

target[1:5,1:4]

g=as.data.frame(rownames(x))

rownames(g)=rownames(x)

colnames(g)=c("gene_short_name")


library(monocle)

pd <- new("AnnotatedDataFrame", data = thymus_subset@meta.data)

fd <- new("AnnotatedDataFrame", data = g)

HSMM <- newCellDataSet(as.matrix(target), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)

HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 0.1)

print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

disp_table <- dispersionTable(HSMM)

ordering_genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 2 * dispersion_fit)$gene_id

HSMM <- setOrderingFilter(HSMM, ordering_genes)

plot_ordering_genes(HSMM)

HSMM <- reduceDimension(HSMM, max_components=2)

HSMM <- orderCells(HSMM, reverse=FALSE)

plot_cell_trajectory(HSMM, show_tree = FALSE ,show_backbone = FALSE, color_by="orig.ident", show_branch_points = FALSE, cell_size=1) + scale_color_manual(values= c("darkorchid","RED1", "royalblue4", "goldenrod2","turquoise3","palegreen2","orangered2","navy","maroon4","deepskyblue","ivory4","royalblue4","magenta2","orangered2","bisque3","gold3","lemonchiffon3","ivory4","royalblue4","lightsalmon4","indianred3","mediumpurple4","mintcream","olivedrab3","powderblue","plum4","sienna3","thistle4"))

#------------------------------------------------------------------------------------------------------------------------
###coordination change
#------------------------------------------------------------------------------------------------------------------------

HSMM2 = HSMM

set.seed(1)

coor<-runif(table(HSMM2$orig.ident)[1],1,2)
coor <- as.data.frame(coor)
coor$V2 = names(table(HSMM2$orig.ident))[1]
coor2 = coor

for (i in 2:9) {
  coor = runif(table(HSMM2$orig.ident)[i],i,i+1)
  coor <- as.data.frame(coor)
  coor$V2 = names(table(HSMM2$orig.ident))[i]
  coor2 = rbind(as.matrix(coor2), as.matrix(coor)) 
}

coor2 = as.data.frame(coor2)
dim(HSMM2)
dim(coor2)

table(HSMM2$orig.ident)
coor2$V2 = factor(coor2$V2, levels = names(table(HSMM2$orig.ident)))
table(coor2$V2)

head(t(HSMM2@reducedDimS))


a = as.data.frame(t(HSMM2@reducedDimS))
head(a)

a2 = as.data.frame(HSMM2@phenoData@data)
head(a2)

a2$coor = "a"

a2[a2$orig.ident%in%names(table(HSMM2$orig.ident))[1],]$coor  = coor2[coor2$V2%in%names(table(HSMM2$orig.ident))[1],]$coor

for ( i in 2:9) {
  a2[a2$orig.ident%in%names(table(HSMM2$orig.ident))[i],]$coor  = coor2[coor2$V2%in%names(table(HSMM2$orig.ident))[i],]$coor
}


a2$coor2 = 1

a2$coor = as.numeric(a2$coor)
class(a2$coor)

head(a2)

HSMM2@reducedDimS[1,] = as.numeric(a2[,9])
HSMM2@reducedDimS[2,] = as.numeric(a2[,10])

HSMM2$Pseudotime = as.numeric(a2[,9])

plot_cell_trajectory(HSMM2, show_tree = FALSE ,show_backbone = FALSE, color_by="orig.ident", show_branch_points = FALSE, cell_size=1) 
plot_cell_trajectory(HSMM2, show_tree = FALSE ,show_backbone = FALSE, color_by="Pseudotime", show_branch_points = FALSE, cell_size=1) 



HSMM_expressed_genes <- row.names(subset(fData(HSMM2), num_cells_expressed >= 10))

HSMM_filtered <- HSMM2[HSMM_expressed_genes,]

diff_test_res <- differentialGeneTest(HSMM_filtered, fullModelFormulaStr="~sm.ns(Pseudotime)")

diff=diff_test_res[,c("gene_short_name", "pval", "qval","use_for_ordering")]

diff=diff[diff$qval < 0.01, ]

diff=diff[order(diff$qval),]


genelist <- row.names(diff_test_res[order(diff_test_res$qval)[1:50],])

genelist <- c("Hes1", "Trdc", "Ptcra", "Cd4", "Cd8a", "Pcna", "Mki67" ,"Ccnd3", "Trac", "Ccr9", "Satb1", "Notch1", "Rag1", "Rag2", "Rorc", "Ets2", "Lef1", "Gata3", "Runx3", "Zbtb7b",  "Id2")

genelist <- c("Runx2", "Hes1",	"Runx1", "Notch1",	"Rag1",	"Rag2",	"Ets2", "Rorc",	"Gata3", "Tox", "Lef1",	"Zbtb7b", "Trgc2", "Trdc", "Ptcra", "Mki67", "Ccnd3", "Cd4", "Cd8a",  "Trac", "Ccr7","Gm4632",	"2700038G22Rik",	"Lockd",	"Gm15340", "Gm28112",	"Gm14718",	"A630023P12Rik")


library(pals)
plot_pseudotime_heatmap(HSMM2[genelist,], num_clusters = 1, cores =1, show_rownames = T,  hmcols = coolwarm(100), cluster_rows = FALSE)





##############################################################################################
### combined tissues
##############################################################################################

library(Seurat)

tissue_tissue<- tissue(heart, y = c(intestine, kidney, lung, liver, thymus), project = "lncRNA")
tissue_tissue<- FindVariableFeatures(tissue5, selection.method = "vst", nfeatures = 2000)
tissue_tissue<- ScaleData(tissue5, vars.to.regress = c('nCount_RNA', 'percent.mt', 'nFeature_RNA'), verbose = TRUE )
tissue_tissue<- RunPCA(tissue5, npcs = 50, verbose = FALSE)

ElbowPlot(object = tissue5, ndims = 50 )

tissue_tissue<- FindNeighbors(tissue5, dims = 1:40)
tissue_tissue<- RunUMAP(tissue5, dims = 1:40)

#---------------------------------------------------------------------------------------------------------------------------------------------------------
### batch correction
#---------------------------------------------------------------------------------------------------------------------------------------------------------

library(SeuratWrappers)

lncrna.list)<- SplitObject(tissue5, split.by="orig.ident")

for (i in 1:length(lung.list)) {
  lncrna.list)[[i]] <- NormalizeDatalncrna.list)[[i]],assay="RNA",normalization.method="LogNormalize")
  lncrna.list)[[i]] <- FindVariableFeatureslncrna.list)[[i]], selection.method="vst", nfeatures=2000)
}

tissue_tissue<- RunFastMNN(object.list = lncrna.list)
tissue_tissue<- RunUMAP(tissue5, reduction = "mnn", dims = 1:40)
tissue_tissue<- FindNeighbors(tissue5, reduction = "mnn", dims = 1:40)
tissue_tissue<- FindClusters(tissue5, resolution = 1)

DimPlot(tissue5, label = TRUE, label.size = 5, group.by = "orig.ident")  


tissue5.marker <- FindAllMarkers(tissue5, only.pos = TRUE, min.diff.pct = 0.2)
tissue5.marker2 = tissue5.marker[tissue5.marker$p_val_adj< 0.01, ]
tissue5.marker2 = tissue5.marker2 tissue5.marker2 avg_log2FC > 0.5, ]


markers = c("Lrp2", "Slc5a2", "Cyp2e1",	"Apob", "Slc12a1", "Umod",	"Myh6",	"Tnnt2",	"Mki67", "Olfm4", "Scgb3a2", "Ccdc113",  "Ager",	"Pdpn",	"Sftpb",	"Sftpa1", "Dcn", "Lum", "Col1a1",  "Nphs1", "Nphs2", 
             "Gata3","Pdgfrb", "Ren1",	"Akr1b7", "Kcnj8",	"Colec11", 	"Mustn1",	"Myh11", "Emcn",	"Kdr", "Pecam1",	"Cldn5", "Lyve1", "Cd4", "Cd8a", "Cd3d", "Cd3e", "Klrb1c", "Nkg7",
             "Cd79a", "Cd79b", "Tnfrsf17",	"Jchain",	"Csf1r", "C1qa", "Adgre1", "Cd68", "Abcg1", "Cd14", "Itgam", "Hbb-bt",	"Hba-a1" )

DotPlot(tissue5, col.min = 0, features = markers, cols = c("grey", "red"), scale = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


ident = c("PT",	"hepatocyte1",	"hepatocyte2",	"kidney_other_tubule",	"Cardiomyocyte",	"Proli_ISC",	"Ciliated",	"AT1",	"AT2",	"Fibroblast2",	"Fibroblast1",	"Fibro3",	"podocyte",	"mesangial",	"JGC",	"Pericyte",	"Myofibro",	"GEC",	"Endo1",	"Endo2",	"LE",	"DN_DP",	"T cell",	"NK",	"B cell",	"B cell2",	"plasma1",	"Kupffer",	"Myeloid1",	"AM",	"myeloid",	"RBC")
tissue5@active.ident = factor(tissue5@active.ident, levels = ident)

col = c( "yellowgreen",	"darkorchid",	"RED1",	"#FDB462", "royalblue4",	"indianred3",	"goldenrod2",	"firebrick3",	"palegreen2",	"orangered2",	"#8DD3C7", 	"maroon4","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#BC80BD", "#B3CDE3", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC" , "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2")

DimPlot(tissue5, label = TRUE, cols = col)

