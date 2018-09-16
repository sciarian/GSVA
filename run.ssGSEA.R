load('server-side//RData//lincs_signatures_cmpd.RData')
load('server-side//RData//lincs_signatures_shrna.RData')
require(GSVA)
require(GSA)

#gurantee that rows are in the same order
sum(rownames(lincs_signatures_cmpd) == rownames(lincs_signatures_shrna))

msigdb          <-  GSA.read.gmt("client-side/meta.data/h.all.v6.1.entrez.gmt")
genesets        <-  msigdb$genesets
names(genesets) <-  msigdb$geneset.names 
combined.expr.matrix    <- cbind(lincs_signatures_cmpd,lincs_signatures_shrna)
intersected.gene        <- intersect(genesets[['HALLMARK_INTERFERON_GAMMA_RESPONSE']],rownames(lincs_signatures_cmpd)) # Great! 33 genes are covered by L1000
intersected.gene        <- intersect(genesets[['HALLMARK_INTERFERON_ALPHA_RESPONSE']],rownames(lincs_signatures_cmpd)) # Oh, only 8 genes are covered by L1000

oncogenic.geneset.gsea.results    <-  gsva(combined.expr.matrix, genesets[c('HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_INTERFERON_ALPHA_RESPONSE')], method = 'ssgsea') #ggsea
hist(oncogenic.geneset.gsea.results[1,],breaks=200)
boxplot(oncogenic.geneset.gsea.results[1,])
IFN.score                         <- oncogenic.geneset.gsea.results[1,]

outlier.upper <- quantile(IFN.score)[4] + 1.5 * IQR(IFN.score)
outlier.lower <- quantile(IFN.score)[2] - 1.5 * IQR(IFN.score)
IFN.score[IFN.score > outlier.upper ]
IFN.score[IFN.score < outlier.lower ]

load('server-side/RData/signature_meta.RData')

flag     <- match(x=as.character(names(IFN.score)),table = signature_meta$id)
new.meta <- signature_meta[flag,]

df <- data.frame(score=IFN.score,cell.line=new.meta$cell_id)
require(ggplot2)
ggplot(df) + geom_boxplot(aes(x=cell.line,y=score))

## pick out perturbations
x    <- IFN.score[IFN.score > outlier.upper ] %>% sort
id   <- names(x)
flag <- match(x=as.character(id),table = signature_meta$id)
tmp  <- signature_meta[flag,]

## should correct base-line level accroding to cell line
