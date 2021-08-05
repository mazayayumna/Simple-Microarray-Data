library(GEOquery)
library(affy)
library(hgu133a.db)
library(hgu133acdf)

setwd("D:/KULIAHAN/Belajar/4. Apply/Nalagenetics/microarray")
f <- list.files(pattern = ".CEL")
j=0
i=0
for(i=1, i<=j, j++ in i:9)
{
	raw.data <- ReadAffy(verbose=FALSE, filenames=f[i], cdfname="hgu133acdf")
	data.rma.norm = rma(raw.data)
}
rma = exprs(data.rma.norm)
write.table(rma, file="rma.txt", quote=FALSE, sep="\t")

probes=row.names(rma)
ls("package:hgu133a.db")
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound=NA))
 rma=cbind(probes,Symbols,Entrez_IDs,rma)
write.table(rma, file = "annotation.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)