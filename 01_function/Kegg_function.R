library(GOstats)
library(KEGG.db)
library(org.Mm.eg.db)
library(Category)
library(pathview)

# --- GO AND KEGG ENRICHMENT ---
myGK <- function (geneId) {
  mygk <- list()
  
  entrezId <- mget(geneId, org.Mm.egSYMBOL2EG, ifnotfound = NA)
  entrezId <- entrezId[! is.na(entrezId)]
  entrezId <- as.character(entrezId)
  
  goAnn <- get("org.Mm.egGO")
  universe <- Lkeys(goAnn)
  for (category in c("BP", "MF", "CC")) {
    params <- new("GOHyperGParams", geneIds = entrezId, universeGeneIds = universe, annotation = "org.Mm.eg.db", 
                  ontology = category, pvalueCutoff = 0.001, testDirection = "over")  
    over = hyperGTest(params)
    go <- summary(over)
    glist <- geneIdsByCategory(over)
    glist <- sapply(glist, function(x) {y <- mget(x, envir=org.Mm.egSYMBOL); paste(y, collapse=";")})
    
    go$Symbols <- glist[as.character(go[, 1])]
    
#   if (nrow(go) > 10) go <- go[1:10, ]
    mygk[[category]] <- go
  }
 
  keggAnn <- get("org.Mm.egPATH")
  universe <- Lkeys(keggAnn)
  params <- new("KEGGHyperGParams", 
                geneIds=entrezId, universeGeneIds=universe, annotation="org.Mm.eg.db", 
                categoryName="KEGG", pvalueCutoff=0.01, testDirection="over")
  over <- hyperGTest(params)
  kegg <- summary(over)
  glist <- geneIdsByCategory(over)
  glist <- sapply(glist, function(x) {y <- mget(x, envir=org.Mm.egSYMBOL); paste(y, collapse=";")})
  kegg$Symbols <- glist[as.character(kegg$KEGGID)]
  mygk[["KEGG"]] <- kegg
  
  return(mygk)
}
