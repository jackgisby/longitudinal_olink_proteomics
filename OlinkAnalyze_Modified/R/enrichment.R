#' prot to go terms
#' 
#' Maps gene symbols to GO ids using biomarRt
#' 
#' @export
#' @import topGO biomaRt

prot_to_go <- function(gene_names) {
    # connected to mart
    mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="useast.ensembl.org")
    
    # get mappings
    hgnc_go <- biomaRt::getBM(attributes = c('hgnc_symbol', "go_id"), 
                              filters = 'hgnc_symbol', 
                              values = gene_names, 
                              mart = mart)
    
    hgnc_go <- hgnc_go[hgnc_go$go_id != "",]
    
    # convert mappings to list required by topGO
    go_annot_list <- list()
    for (prot in unique(hgnc_go$hgnc_symbol)) {
        go_annot_list[[prot]] = hgnc_go$go_id[hgnc_go$hgnc_symbol == prot & !is.na(hgnc_go$go_id) & hgnc_go$go_id != ""]
    }
    
    return(go_annot_list)
}

#' colour mapping for topGO
#' 
#' for some reason, topGO is missing a function - here it is
#' 
#' @export
#' @import topGO

colMap <- function(x) {
    .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
    return(.col[match(1:length(x), order(x))])
}

#' run and plot GO analysis using topGO
#' @export
#' @import topGO

run_go <- function(pvals, go_annot_list, p_label=2.7, size_label=1, 
                   plt_title=NULL, return_res=FALSE, test_stat="ks") {
    
    # create a "topGOdata" class to begin using package
    sampleGOdata <- new("topGOdata", ontology = "BP",
                        allGenes = pvals, annotationFun = annFUN.gene2GO,
                        gene2GO=go_annot_list,
                        geneSelectionFun = function(allScore) {return(allScore < 0.05)},  # <0.05 significance
                        nodeSize = 10)
    
    # run an enrichment test using weight01 and the specified test stat
    pvals <- runTest(sampleGOdata, algorithm = "weight01", statistic = test_stat)
    
    # convert the results into a table
    allRes <- GenTable(sampleGOdata, pvals = pvals, orderBy = "pvals", ranksOf = "pvals", 
                       topNodes = length(usedGO(sampleGOdata)))
    
    # calculate other measures of interest
    allRes$gSize <- log2(allRes$Significant / allRes$Expected)
    
    allRes$pvals <- as.numeric(allRes$pvals)  # no fdr!
    allRes$logp <- -log10(allRes$pvals)
    
    print(ggplot(allRes, aes(x=gSize, y=logp, colour=pvals < 0.05)) +
        geom_point(alpha=0.5) +
        ylab("-log10(Pvalue)") +
        xlab("Log2 Fold Change (Significant Annotations - Expected Annotations)") +
        ggtitle(plt_title) +
        geom_text_repel(data=subset(allRes,
                                    logp > p_label | 
                                        (gSize > size_label & pvals < 0.05)),
                        aes(gSize, logp, label = Term), size = 3, color="steelblue"))
    
    return(allRes)
}

#' send topGO GO terms to cytoscape
#' 
#' can be useful to visualise the topGO terms in a network - so here's a function
#' for exporting results to cytoscape
#' 
#' there's actually a built in function for creating a plot based on this network,
#' but it was not suitable as there were so many terms
#' 
#' @export
#' @import RCy3

make_cyt_graph <- function(pvals, go_annot_list, test_stat = "ks", graph_name="myGraph", firstSigNodes=10) {
    sampleGOdata <- new("topGOdata", ontology = "BP",
                        allGenes = pvals, annotationFun = annFUN.gene2GO,
                        gene2GO=go_annot_list,
                        geneSelectionFun = function(allScore) {return(allScore < 0.05)}, 
                        nodeSize = 10)
    
    # re-calculate pvals etc
    pvals <- runTest(sampleGOdata, algorithm = "weight01", statistic = test_stat)
    
    allRes <- GenTable(sampleGOdata, pvals = pvals, orderBy = "pvals", ranksOf = "pvals", 
                       topNodes = length(usedGO(sampleGOdata)))
    
    # generate the unsuitable inbuilt topgo plot
    sig_graph <- showSigOfNodes(sampleGOdata, score(pvals), firstSigNodes = firstSigNodes, useInfo = "np")
    
    # extract the relevant information from the topgo object
    nodeDataDefaults(sig_graph$dag, "p") <- 1
    nodeDataDefaults(sig_graph$dag, "desc") <- ""
    
    for (n in nodes(sig_graph$dag)) {
        nodeData(sig_graph$dag, n, "p") <- as.numeric(allRes$pvals[allRes$GO.ID == n])
        nodeData(sig_graph$dag, n, "desc") <- allRes$Term[allRes$GO.ID == n]
    }
    
    # send the extracted tree to cytoscape
    createNetworkFromGraph(sig_graph$dag, graph_name)
    
    return(sig_graph$dag)
}

#' prot to kegg terms
#' 
#' Maps gene symbols to GO ids using biomarRt
#' 
#' @export
#' @import topGO biomaRt

prot_to_entrez <- function(gene_names, kegg_rno) {
    # connected to mart
    mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="useast.ensembl.org")
    
    # get mappings
    hgnc_go <- biomaRt::getBM(attributes = c('hgnc_symbol', "entrezgene_id"), 
                              filters = 'hgnc_symbol', 
                              values = gene_names, 
                              mart = mart)
    
    # convert mappings to list required by topGO
    go_annot_list <- list()
    for (prot in unique(hgnc_go$hgnc_symbol)) {
        go_annot_list[[prot]] = hgnc_go$entrezgene_id[hgnc_go$hgnc_symbol == prot & !is.na(hgnc_go$entrezgene_id) & hgnc_go$entrezgene_id != ""]
    }
    
    go_annot_list[["IGLC2"]] <- 3538
    go_annot_list[["GIF"]] <- 2694
    go_annot_list[["FIGF"]] <- 2277
    go_annot_list[["IL8"]] <- 3576
    
    return(go_annot_list)
}

#' run kegg DE
#' 
#' @export
#' @import topGO biomaRt

run_kegg <- function(pvals, entrez_annot_list) {
    
    for (i in 1:length(pvals)) {
        names(pvals)[i] <- entrez_annot_list[[names(pvals)[i]]]
    }
    
    enrichKEGG(names(pvals)[pvals < 0.05], keyType="ncbi-geneid", minGSSize=3, universe=names(pvals),pAdjustMethod="none")
}