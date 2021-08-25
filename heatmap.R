#!/share/nas1/lihz/software/Miniconda3/envs/R35/bin/R

library('getopt');
library('ComplexHeatmap')
library('psych')
library('circlize')

##########
spec = matrix(c(
       'help' , 'h', 0, "logical",
       'micro' , 'm', 1, "character",
       'metab' , 'b', 1, "character",
       'samples' , 's', 1, "character",
       'outfile' , 'o', 1, "character",
       'top' , 't', 1, "integer"
        ), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
        cat(getopt(spec, usage=TRUE));
        cat("
Usage example:
        Rscript  heatmap.r  -m micro.txt -b metabolite.txt -s AA1,AA2,AA3 -o heatmapFile
Options:
-h              NULL            get this help
-m      character       the file for microbita abudance [forced:firstLine is header,each colum separate by tab]
-b      character       the file for metabolite abudance [forced:firstLine is header,each colum separate by tab]
-s      character       the sample for display [forced:sample1,sample2,sample3...]
-t      integer         show Top microbita or metabolite [defaults:15]
-o      character       the filename for output graph [forced]
\n")
        q(status=1);
}

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }
# check non-null args
if ( is.null(opt$micro) )      {print_usage(spec) }
if ( is.null(opt$metab) )      {print_usage(spec) }
if ( is.null(opt$samples) )    {print_usage(spec) }
if ( is.null(opt$outfile) )    {print_usage(spec) }
if ( is.null(opt$top) )        {opt$top = 15}

#reading data
samples <- strsplit(opt$samples,',')[[1]]
microData <-read.table(file=opt$micro, header = T,comment.char="", sep="\t",row.names = 1)
metabData <- read.table(file=opt$metab, header = T,comment.char="", sep="\t",row.names = 1)

# choose sample data
if ((length(intersect(samples,colnames(microData))) != length(samples)) || (length(intersect(samples,colnames(metabData))) != length(samples))){
  print("sample not in microbita or metabolite file ,please check...")
  q(status=1)
}

## get top abudance microbita or metabolite
top <- function(data,topValue=15){
  data <- data[order(rowSums(data),decreasing=TRUE),]
  if (topValue > nrow(data)) {topValue <- nrow(data)}
  data <- data[1:topValue,]
  return(data)
}
microAnalysisData <- top(microData[,samples], topValue=opt$top)
metabAnalysisData <- top(metabData[,samples], topValue=opt$top)

###### calculate correction
corr <- corr.test(t(microAnalysisData), t(metabAnalysisData),method="pearson",adjust="fdr")
significant <- function(x) cut(x,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
pvalueLabel <- t(apply(corr$p,1,significant))
correction <- corr$r

########### z-score ###########
microZscoreData <- t(apply(microAnalysisData, 1, scale))
colnames(microZscoreData) <- samples
metabZscoreData <- t(apply(metabAnalysisData, 1, scale))
colnames(metabZscoreData) <- samples


### draw
col_fun1 = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun2 = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))

p1 <- Heatmap(correction,col=col_fun1,show_heatmap_legend=F,
              bottom_annotation = HeatmapAnnotation(foo=metabZscoreData,col = list(foo = col_fun2),
                                                    annotation_name_gp = gpar(fontsize=8),
                                                    show_legend = F),
              cell_fun = function(j, i, x, y, width, height, fill) {grid.text(pvalueLabel[i, j], x, y, gp = gpar(fontsize = 10))},
              column_title = 'Correction Coefficient',
              column_title_gp = gpar(fontsize = 10),
              width = unit(6, "cm"), height = unit(8, "cm"),
              row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8))
p2 <- Heatmap(microZscoreData,col=col_fun2,show_heatmap_legend=F,
              column_title = 'Microbiome',column_title_gp = gpar(fontsize = 10),
              width = unit(6, "cm"), height = unit(8, "cm"),
              row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 8))

lgd1 = Legend(title_gp = gpar(fontsize = 10),title = "Correction", col_fun = col_fun1, at = c(1, 0, -1), labels = c("1", "0", "-1"))
lgd2 = Legend(title_gp = gpar(fontsize = 10),title = "Microbiome", col_fun = col_fun2, at = c(max(microZscoreData), 0, min(microZscoreData)), labels = c(ceiling(max(microZscoreData)), "0", floor(min(microZscoreData))))
lgd3 = Legend(title_gp = gpar(fontsize = 10),title = "Metabolome", col_fun = col_fun2, at = c(max(metabZscoreData), 0, min(metabZscoreData)), labels = c(ceiling(max(metabZscoreData)), "0", floor(min(metabZscoreData))))
pdf(file = paste(opt$outfile,".pdf",sep = ""),width = 16, height = 8)
draw(p1+p2, ht_gap = unit(1, "cm"))
draw(lgd1,x = unit(0.9, "npc"), y = unit(0.6, "npc"))
draw(lgd2,x = unit(0.9, "npc"), y = unit(0.5, "npc"))
draw(lgd3,x = unit(0.9, "npc"), y = unit(0.4, "npc"))
dev.off()
