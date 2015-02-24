#function to map ensemblids to gene symbols
ensembl2symbol <- function(x){
  require(org.Hs.eg.db)
  foo <- as.list(org.Hs.egENSEMBL2EG)
  keep <- which(x%in%names(foo))
  bar <- rep(NA,length(x))
  bar[keep] <- sapply(foo[x[keep]],function(x) x[1])
  baz <- org.Hs.egSYMBOL
  # Get the gene symbol that are mapped to an entrez gene identifiers
  mapped_genes <- mappedkeys(baz)
  # Convert to a list
  bax <- as.list(baz[mapped_genes])
  bar[keep] <- unlist(bax[bar[keep]])
  return(bar)
}