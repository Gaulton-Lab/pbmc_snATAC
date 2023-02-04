suppressPackageStartupMessages(library(cicero))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(tidyr))

clus <- makeCluster(32)

wd = '/nfs/lab/projects/pbmc_snATAC/analysis_v2/run_cicero_paola/'
setwd(wd)
args<-commandArgs(TRUE)

input_mat  = args[1]
out_prefix = args[2]
pbmcid     = args[3]

data           = read.table( input_mat ,header=T,  sep="\t")
colnames(data) = c("peak", "barcode", "value")

sc.umap = read.table( "../peaks/pbmc1-15_clusterLabels.txt", header=T, row.names=1,  sep="\t") ### map barcode to cluster 
rownames(sc.umap) = sc.umap$X



### subset to one sample
sc.umap = droplevels(subset(sc.umap, sample == pbmcid))
data    = droplevels(subset(data, barcode %in% rownames(sc.umap)))


### convert to sparse matrix peak x barcode
sc.data <- with(data, sparseMatrix(i=as.numeric(as.factor(peak)), j=as.numeric(as.factor(barcode)), 
                                   x=value, dimnames=list(levels(as.factor(peak)), levels(as.factor(barcode)))))
rownames(sc.data) <- paste0('chr', gsub('-','_', gsub(':','_',rownames(sc.data))))


sc.data.subset <- sc.data


cellinfo <-data.frame(cells=colnames(sc.data.subset))
row.names(cellinfo) <- cellinfo$cells
dhsinfo <- data.frame(site_name=rownames(sc.data.subset))
row.names(dhsinfo) <- dhsinfo$site_name
dhsinfo <- cbind(dhsinfo, stringr::str_split_fixed(dhsinfo$site_name, "_", 3))
names(dhsinfo) <- c('site_name','chr','bp1','bp2')
dhsinfo$chr <- gsub('chr','', dhsinfo$chr)
dhsinfo$bp1 <- as.numeric(as.character(dhsinfo$bp1))
dhsinfo$bp2 <- as.numeric(as.character(dhsinfo$bp2))

input_cds <- suppressWarnings(newCellDataSet(as(sc.data.subset, 'dgCMatrix'),
                                             phenoData = methods::new('AnnotatedDataFrame', data = cellinfo),
                                             featureData = methods::new('AnnotatedDataFrame', data = dhsinfo),
                                             expressionFamily=negbinomial.size(),
                                             lowerDetectionLimit=0))
input_cds@expressionFamily <- binomialff()
input_cds@expressionFamily@vfamily <- 'binomialff'
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

input_cds <- input_cds[fData(input_cds)$num_cells_expressed > 0,]
umap_coords <- sc.umap[colnames(sc.data.subset), c('UMAP1','UMAP2')]
colnames(umap_coords) <- NULL

cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords, k=30)
window <- 1e6
data('human.hg19.genome')
distance_parameters <- estimate_distance_parameter(cicero_cds, window=window, maxit=100, sample_num=100, distance_constraint=500000, genomic_coords=human.hg19.genome)
mean_distance_parameter <- mean(unlist(distance_parameters))
cicero_out <- generate_cicero_models(cicero_cds, distance_parameter=mean_distance_parameter, window=window, genomic_coords=human.hg19.genome)
conns <- assemble_connections(cicero_out, silent=FALSE)
#saveRDS(conns, file.path(wd, paste0('cicero/', cluster, '.1MB_cicero_conns.rds')))
#write.table(conns, file.path(wd, paste0('output/', celltype, '.cicero_conns_dedup.txt')), sep='\t', quote=FALSE, row.names=FALSE)

## this step is to remove duplicated connections
conns = conns[order(-conns$coaccess),]
bed = cbind(str_split_fixed(conns[,1], "\\_", 3 ), str_split_fixed(conns[,2], "\\_", 3 ))

ord = matrix(parRapply(clus, bed, function(x) x[order(as.numeric(x))] ), ncol=6, byrow=T)

ord = cbind(ord[, c(5,1:2,5,3:4)], conns$coaccess)
dedup = ord[!duplicated(ord[,1:6]),]                      

dist = as.numeric(dedup[,6])-as.numeric(dedup[,2])
dedup  = subset(dedup, dist >10000)
                       
dedup = data.frame( Peak1 = paste(dedup[,1], dedup[,2], dedup[,3], sep="_")  , 
                    Peak2 = paste(dedup[,4], dedup[,5], dedup[,6], sep="_") , coaccess = dedup[,7]  )


write.table(dedup, file.path(wd, paste0('output_per_sample/', out_prefix, '.cicero_conns_dedup.txt')), sep='\t', quote=FALSE, row.names=FALSE)

writeLines(c("tot_pairs", nrow(dedup),"coaccess", sum(as.numeric(dedup[,3])>0.05, na.rm=T)), paste0('output_per_sample/', out_prefix,"summary"))


###################################
####### Compare PCHIC #######
###################################
setwd(file.path(wd, 'output_per_sample'))

#### Promoter Cature HiC Primary cells
dir = "/nfs/lab/projects/pbmc_snATAC/data/publicdata/"
pc = read.table(paste0( dir, 'PCHiC_peak_matrix_cutoff5.tsv'), header=T)

cells = colnames(pc)[14:ncol(pc)-2]
wcel = data.frame(Peak1 = paste0("chr", pc$baitChr,"_" , pc$baitStart, "_" , pc$baitEnd),
                  Peak2 = paste0("chr", pc$oeChr,"_" , pc$oeStart, "_" , pc$oeEnd),
                 pc[,cells],
                 dist = abs(pc$dist))
wcel = subset(wcel, dist <  1050000)
wcel = subset(wcel, dist >  10000)




load_filtered_connections = function( cluster) {

infile   = paste0(cluster,   '.cicero_conns_dedup.txt')
outfile  = paste0(cluster, '.cicero_conns_dedup.bedpe')
baits    = '/nfs/lab/projects/pbmc_snATAC/data/publicdata/Digest_Human_HindIII_baits.bed'
filtfile = paste0(cluster, '.cicero_conns_baits')

    if(!file.exists(filtfile)){
#sed = paste("tail -n +2", infile, "| sed 's/_/\\t/g' | pgltools sort >" , outfile) ## double escape
#sed = paste( "sed 's/_/\\t/g'", infile, "| pgltools sort >" , outfile) ## double escape
sed = paste("awk -F \"\\t\" '{print $1,$2,$3}'", infile, "| tail -n +2 | sed 's/_/\\t/g' | pgltools sort >" , outfile) 
        
pgl = paste('pgltools intersect1D -a', outfile ,'-b' , baits,  '-wa -d 1000 >', filtfile)

system(sed)
system(pgl)
}
conns = fread(filtfile, data.table=F, header=F, sep="\t")
conns = conns[!duplicated(conns[,1:6]),]
    
cic = data.frame(Peak1 = paste(conns[,1], conns[,2], conns[,3], sep="_"), 
                 Peak2 = paste(conns[,4], conns[,5], conns[,6], sep="_"), coaccess = conns[,8])
cic$dist = conns[,6]- conns[,2]
return(cic)
}


cic = load_filtered_connections( cluster=out_prefix )
cic = subset(cic, dist>10000)

cic$signif <- cic$coaccess > 0.05 

test = matrix(unlist(mclapply(cells, function(x) compare_connections(cic, data.frame(wcel[wcel[,x]>=5,1:2]), maxgap=1000),
                  mc.cores =18)), ncol=length(cells))
colnames(test) = cells   
    
tab = matrix(unlist(mclapply(cells, function(x)fisher.test(rbind(table(test[cic$signif==FALSE, x]), 
                                              table(test[cic$signif, x])))[c("p.value", "estimate")] )), nrow=2)
colnames(tab)=cells
rownames(tab) = c("pv.pos", "or.pos")
write.table (tab, paste0(out_prefix, ".comparePCHIC.fisher"),  sep="\t", quote=F)

