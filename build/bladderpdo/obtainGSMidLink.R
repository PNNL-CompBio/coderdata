gse <- GEOquery::getGEO("GSE103990",GSEMatrix=FALSE)
GEOquery::Meta(gse)
GSMs <- GEOquery::Meta(gse)$sample_id
gsmlinkDf <- data.frame(GSM = GSMs, sampleid = NA, RNAid = NA)
for (i in 1:length(GSMs)){
  GSMinfo <- GEOquery::getGEO(gsmlinkDf[i,"GSM"])
  #description + title
  gsmlinkDf[i,'sampleid'] <- GEOquery::Meta(GSMinfo)$title
  gsmlinkDf[i, 'RNAid'] <- GEOquery::Meta(GSMinfo)$description[2]
}

write.csv(gsmlinkDf, file="./gsmlinkDf.csv", row.names = F)
