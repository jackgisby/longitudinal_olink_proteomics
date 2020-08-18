#' prot to go term
#' @export
#' @import topGO biomaRt

prot_to_go <- function(gene_names) {
    mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="useast.ensembl.org")
    
    hgnc_go <- biomaRt::getBM(attributes = c('hgnc_symbol', "go_id"), 
                              filters = 'hgnc_symbol', 
                              values = names(cctrl_pvals), 
                              mart = mart)
    
    hgnc_go <- hgnc_go[hgnc_go$go_id != "",]
    
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

run_go <- function(pvals, go_annot_list, p_label=2.7, size_label=c(0.9, 1.5), 
                   plt_title=NULL, return_res=FALSE, test_stat="ks") {
    
    sampleGOdata <- new("topGOdata", ontology = "BP",
                        allGenes = pvals, annotationFun = annFUN.gene2GO,
                        gene2GO=go_annot_list,
                        geneSelectionFun = function(allScore) {return(allScore < 0.05)},  # <0.05 significance
                        nodeSize = 10)
    
    ks_p <- runTest(sampleGOdata, algorithm = "weight01", statistic = test_stat)
    
    allRes <- GenTable(sampleGOdata, ks_p = ks_p, orderBy = "ks_p", ranksOf = "ks_p", 
                       topNodes = length(usedGO(sampleGOdata)))
    
    allRes$gSize <- allRes$Significant / allRes$Expected
    
    allRes$ks_p <- as.numeric(allRes$ks_p)
    allRes$logp <- -log10(allRes$ks_p)
    
    if (return_res) {
        return(allRes)
    }
    
    ggplot(allRes, aes(x=gSize, y=logp, colour=ks_p < 0.05)) +
        geom_point(alpha=0.5) +
        ylab("-log10(adj_p)") +
        xlab("Significant proteins with term / Expected proteins with term") +
        ggtitle(plt_title) +
        geom_text_repel(data=subset(allRes,
                                    logp > p_label | 
                                        ((gSize < size_label[1] | gSize > size_label[2]) & ks_p < 0.05)),
                        aes(gSize, logp, label = Term), size = 3, color="steelblue") +
        geom_vline(xintercept =  1)
}

#' send topGO GO terms to cytoscape
#' @export
#' @import RCy3

make_cyt_graph <- function(pvals, go_annot_list, test_stat = "ks", graph_name="myGraph", firstSigNodes=10) {
    sampleGOdata <- new("topGOdata", ontology = "BP",
                        allGenes = pvals, annotationFun = annFUN.gene2GO,
                        gene2GO=go_annot_list,
                        geneSelectionFun = function(allScore) {return(allScore < 0.05)}, 
                        nodeSize = 10)
    
    ks_p <- runTest(sampleGOdata, algorithm = "weight01", statistic = test_stat)
    
    allRes <- GenTable(sampleGOdata, ks_p = ks_p, orderBy = "ks_p", ranksOf = "ks_p", 
                       topNodes = length(usedGO(sampleGOdata)))
    
    sig_graph <- showSigOfNodes(sampleGOdata, score(ks_p), firstSigNodes = firstSigNodes, useInfo = "np")
    
    nodeDataDefaults(sig_graph$dag, "p") <- 1
    nodeDataDefaults(sig_graph$dag, "desc") <- ""
    
    for (n in nodes(sig_graph$dag)) {
        nodeData(sig_graph$dag, n, "p") <- as.numeric(allRes$ks_p[allRes$GO.ID == n])
        nodeData(sig_graph$dag, n, "desc") <- allRes$Term[allRes$GO.ID == n]
    }
    
    createNetworkFromGraph(sig_graph$dag, graph_name)
    
    return(sig_graph$dag)
}
