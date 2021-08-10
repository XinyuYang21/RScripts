# Multi-layer network
mypal <- colorRampPalette(pal_nejm("default", alpha = 1)(8))(18)
pathway_col = data.frame(Pathway=c("Cell_Cycle", "Chromatin_modifiers", "Notch", "NRF2", "PI3K", "RTK_RAS", "TGF-Beta", "TP53", "WNT", "Hippo", "Myc",NA), col=c(mypal[1:11],"grey70"))


get_layers = function(net) {
  layers_name = unique(V(net)$value) %>% as.numeric() %>% sort()
  layers_id = data.frame(name=layers_name,layers=length(layers_name):1)
  layer = data.frame(name = as.numeric(V(net)$value)) %>% left_join(layers_id,by="name") %>% .$layers
  return(layer)
}


#layout = layout_with_sugiyama(net, layers=get_layers(net))
check_layout_facet <- function(layout,node,facet=NA,layer1_threshold=0.99,simple=FALSE) {
  x_level = max(layout[,1])
  y_level = length(unique((layout[,facet])))
  df = layout;colnames(df) = c("x_level","y","id","y_level","freq",ccf_type)
  
  if (!simple) {
  area = df  %>% filter(x_level==1) %>%group_by(y_level) %>%
    dplyr::summarise(n=n(),max_freq=max(freq),y_width=sum(freq %>% rescale(c(20,30)))+max(freq)*100)  %>%
    mutate(y_max=cumsum(y_width))
  area$y_min =  c(4,area[-nrow(area),]$y_max+area[-nrow(area),]$max_freq*30)
  
  
  # Avoid overlap
  for (i in area$y_level) {
    #for (j in 1:x_level) {
    df[which(df$y_level==i),]=df[which(df$y_level==i),] %>% arrange(freq) %>% arrange(y) %>% mutate(y=seq(area[which(area$y_level==i),]$y_min,area[which(area$y_level==i),]$y_max,length.out=nrow(.)))
    # df_xlevel = df[which(df$y_level==i),][which(df[which(df$y_level==i),]$x_level==j),]
    #if (nrow(df_xlevel)>0) df[which(df$y_level==i),][which(df[which(df$y_level==i),]$x_level==j),]=df_xlevel  %>% arrange(y) %>% mutate(y=seq(min(y),max(y),length.out=nrow(.)))
    #}
  }
  }
  df %>% mutate(x_level=as.numeric(node$median_ccf )) 
}

layout_to_matrix <- function(layout) {layout %>% arrange(id) %>% .[,1:2] %>% as.matrix()}

net_to_layout <- function(net,node,simple=FALSE){
  cbind(get_layers(net),layout_with_sugiyama(net, layers=V(net)$value)$layout[,1]) %>%
    as.data.frame() %>% cbind(node[,c("id","community_id","freq",ccf_type)]) %>%
    mutate(community_id=as.numeric(as.character(community_id))) %>% 
    check_layout_facet(facet="community_id",node=node,simple=simple)  %>%
    layout_to_matrix()
}


ggraph_plot = function(nodes,edges,facet=FALSE,title="",fill=FALSE,zoom=NA,p_value,no_layout=FALSE,subtitle="") {
  
  mypal <- colorRampPalette(pal_nejm("default", alpha = 0.8)(8))(18)
  #scales::show_col(mypal)
  pathway_col <- data.frame(Pathway=c("Cell_Cycle", "Chromatin_modifiers", "Notch", "NRF2", "PI3K", "RTK_RAS", "TGF-Beta", "TP53", "WNT", "Hippo", "Myc",NA), col=c(mypal[c(1:8,10,15:16)],"grey70"))
  pathway_col_values <- pathway_col$col
  names(pathway_col_values) = pathway_col$Pathway
  
  gene_edges <- edges %>% filter(gene1 %in% nodes$label | gene2 %in% nodes$label) 
  from_id <- sapply( gene_edges$gene1,function(x) nodes %>% filter(label==x) %>% .$id ) %>% as.numeric()
  to_id <- sapply( gene_edges$gene2,function(x) nodes %>% filter(label==x) %>% .$id ) %>% as.numeric()
  
  ### get layout
  edges <- data.frame(from = from_id, to = to_id,
                      color = ifelse( gene_edges$Event == "Co_Occurence", "#4058a3",  "#f18826")) %>% cbind(gene_edges) 
  
  co_nodes <- c(t(edges[which(edges$Event == "Co_Occurence"), c("gene1","gene2")]))  %>% unique(.)
  mu_nodes <- c(t(edges[which(edges$Event == "Mutually_Exclusive"), c("gene1","gene2")]))  %>% unique(.)
  #nodes_mu = nodes %>% filter(label %in% mu_nodes) 
  if (length(co_nodes)>0) {
   nodes_co <- nodes %>% filter(label %in% co_nodes)
   edges_co <- edges %>% filter(Event == "Co_Occurence", gene1 %in% co_nodes, gene2 %in% co_nodes)
   nodes_mu <- nodes %>% filter(label %in% mu_nodes | freq==1) 
   edges_mu <- edges %>% filter(Event == "Mutually_Exclusive",gene1 %in% mu_nodes , gene2 %in% mu_nodes)
   net_co <- graph_from_data_frame(d = edges_co, vertices = nodes_co, directed = F) 
  # Merge group of co-occurrence nodes label to whole nodes
  } else {
    nodes_mu <- nodes %>% filter(label %in% mu_nodes| freq==1) 
    edges_mu <- edges %>% filter(Event=="Mutually_Exclusive", gene1 %in% mu_nodes , gene2 %in% mu_nodes) # nolint
    net_co <- graph_from_data_frame(d=edges_mu, vertices=nodes_mu, directed=F)
  }
  # Calculate groups for nodes with co-occurrences relationship
  tg_co <- tidygraph::as_tbl_graph(net_co) %>% activate(nodes) %>% mutate(label = label, name = as.numeric(name), community_id = group_components())
  gene_nodes = left_join(nodes,data.frame(label=V(tg_co)$label,community_id=V(tg_co)$community_id),by="label") %>%
    mutate(community_id=as.factor(ifelse(is.na(community_id),0,community_id)))
  net = graph_from_data_frame(d = edges, vertices = gene_nodes, directed = F)

  tg <- tidygraph::as_tbl_graph(net) %>% activate(nodes) %>% mutate(label=label)
  if (no_layout) l = net_to_layout(net,gene_nodes,simple=TRUE) else
    l = net_to_layout(net,gene_nodes)

  v.size <- V(tg)$freq_value %>% as.character() %>% as.numeric()
  e.color <- E(tg)$color
  E(tg)$weight <- E(tg)$log10pval %>% rescale(c(0.5,2.5))
  eigenCent <- evcent(tg)$vector
  bins <- unique(quantile(eigenCent, seq(0,1,length.out=30)))

  p <- ggraph(tg, layout = l)
  if (nrow(nodes_co)>0) p = p + geom_edge_link0(aes(filter = (Event=="Co_Occurence")),edge_colour = "#4058a3",width = E(tg)$weight[which(E(tg)$Event == "Co_Occurence")]) 
  if (nrow(nodes_mu)>0) p = p + geom_edge_link0(aes(filter = (Event=="Mutually_Exclusive")),edge_colour = "#f18826",width = E(tg)$weight[which(E(tg)$Event == "Mutually_Exclusive")])

  p <- p +
    geom_node_point(aes(color = group), size = v.size/2, alpha = 0.9) +
    ggtitle(title)+
    scale_color_manual(values = pathway_col_values)+
    #geom_node_point(aes(size=as.numeric(as.character(freq_value))),alpha = 0.9) +
    #guides(shape = guide_legend(override.aes = list(size = 15:23)))+
    geom_node_text(aes(label = label), #repel = TRUE,
                   size = log(v.size),colour = "black")+
    theme_classic() +
    scale_x_reverse() +
    labs(x=ccf_type,y = "", subtitle = subtitle) +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none")
  
  if (fill) {p <- p +
    ggforce::geom_mark_hull(aes(filter=community_id!=0,x=x,y=y,fill = as.factor(community_id),label = paste0("trajectory ", community_id)),
                            colour = NA,con.colour = "grey",show.legend = FALSE,concavity = 4,
                            expand = unit(3, "mm"),alpha = 0.1)
  }
  
  if (facet) {
    p = p+facet_nodes(~community_id,scales="free_y",ncol=1)+
      theme(strip.text.x = element_blank(),strip.background = element_blank())
  }
  p
}

mafToNodeEdges <- function(maf,ccf_filter,manSelect,gistic=NULL,ksByCancer,pathway=FALSE,top_gene,oncoplot.feature=NA,fabcolors=NA,oncoPath,bottom_barplot=TRUE,bottom_feature) {
  
  mafTemp   <- maf@data
  clinicalTemp <- getSampleSummary(x = maf) %>%
    left_join(manSelect, by = c("Tumor_Sample_Barcode"="samplename")) %>%
    distinct(Tumor_Sample_Barcode,.keep_all = TRUE)
  # Calculate Sample frequency
  geneccfFreq <- data.table(ccf_filter)[,ccube_ccf:=ifelse(ccube_ccf>=1,1,ccube_ccf)] %>%
      .[samplename %in% clinicalTemp$Tumor_Sample_Barcode & Variant_Classification %nin%  c('Translation_Start_Site') ] %>%
      inner_join(oncoPath, by=c("Hugo_Symbol"="Gene")) 

  if (!is.null(gistic)) gistic = subset(gistic,variable %in% getSampleSummary(x = maf)$Tumor_Sample_Barcode )
  
  if (pathway) {
    mafTemp <- inner_join(mafTemp, oncoPath, by = c("Hugo_Symbol"="Gene"))
    colnames(mafTemp)[c(1,which(colnames(mafTemp)=="Pathway"))] <- c("Hugo_Symbol_Gene","Hugo_Symbol")
    
    mafClinical = read.maf(maf = mafTemp, clinicalData = clinicalTemp)
    res <- getInteractions(mafClinical, top=10, colPal = "PRGn", returnAll = TRUE, sigSymbolsSize=1.5)
    
    geneccfFreq = geneccfFreq %>%
      group_by(Pathway) %>%
      dplyr::summarise(mean_ccf=mean(ccube_ccf,na.rm=TRUE),median_ccf=median(ccube_ccf,na.rm=TRUE),n=length(unique(samplename))) 
    
    bottom.feature=NULL
  } else {
    mafClinical = read.maf(maf = mafTemp, clinicalData = clinicalTemp,cnTable = gistic)
    res <- getInteractions(maf, top=top_gene, colPal = "PRGn", returnAll = TRUE, sigSymbolsSize=1.5) 
 
    geneccfFreq = geneccfFreq %>%
      group_by(Hugo_Symbol) %>%
      dplyr::summarise(mean_ccf=mean(ccube_ccf,na.rm=TRUE),median_ccf=median(ccube_ccf,na.rm=TRUE),n=length(unique(samplename)))  
  }

  geneccfFreq = data.table(geneccfFreq)[,freq:= n/nrow(clinicalTemp)]
  geneccfFreq[,freq_value:=cut(freq,breaks=c(seq(0,1,length.out=10)),labels=seq(15,50,length.out=9))]
  colnames(geneccfFreq)[1] = "label"

  top_genes = geneccfFreq %>% arrange(desc(freq)) %>% .$label %>% .[1:5] 
    # Draw oncoplot
  if (is.null(bottom_feature)) bottom_barplot = FALSE else bottom_barplot = TRUE
  p_oncoplot = as.ggplot(function()
      if (all(!is.na(oncoplot.feature))) {
        oncoplot(maf = mafClinical, top=top_gene, removeNonMutated = FALSE, #colors =  mycolors, 
                 drawRowBar = TRUE, annotationColor= fabcolors,
                 keepGeneOrder = FALSE, bgCol="#F5F5F5",sortByAnnotation = TRUE,
                 clinicalFeatures=oncoplot.feature,bottom_barplot = bottom_barplot,bottom.feature  = bottom_feature,
                 annotationFontSize = 1.1,legendFontSize = 1.1)
                 
      } else {
        oncoplot(maf = mafClinical, top=top_gene, removeNonMutated = FALSE, #colors = mycolors, 
                 annotationColor=fabcolors,  
                 drawRowBar = TRUE,sortByAnnotation = TRUE, keepGeneOrder = FALSE, bgCol="#F5F5F5",bottom_barplot = bottom_barplot,
                 bottom.feature  = bottom_feature)
      }
    )

    # Draw interaction heatmap
    p_interheatmap <- interaction_heatmap(res)
    res$pairs %>%
      arrange(11)
    # Get nodes
    gene_edges <- res$pairs[pValue < p_value]
    gene_nodes <- cbind(gene_edges$gene1, gene_edges$gene2) %>% t(.) %>% as.character() %>% unique(.)
    gene_nodes_col <- data.frame(Pathway = gene_nodes) %>% left_join(pathway_col, by="Pathway")

    if (!pathway) {
      #gene_nodes= c(gene_nodes,top_genes) %>% unique(.)
      gene_nodes= c(gene_nodes) %>% unique(.)
      gene_nodes_col = data.frame(Gene=gene_nodes) %>%  
        left_join(oncoPath,by="Gene") %>%
        left_join(pathway_col,by="Pathway") %>%
        distinct(Gene,.keep_all = TRUE)
    }

    if (length(gene_nodes) > 0) {
      nodes <- data.frame(id = 1:length(gene_nodes),
                          label = gene_nodes,
                          group = gene_nodes_col$Pathway,
                          color = gene_nodes_col$col,
                          shadow = TRUE) %>%
        inner_join(geneccfFreq, by = "label")
      nodes <- nodes[which(nodes[, ccf_type] > 0),] %>%
        mutate(value = cut(.[,ccf_type],breaks = (1:10) / 10,labels = 1:9))
      nodes[, ccf_type] <- round(nodes[, ccf_type], 3)
      gene_edges = gene_edges %>% 
        filter(gene1 %in% nodes$label,gene2 %in% nodes$label)
      return(list(nodes, gene_edges, p_oncoplot, p_interheatmap))
    } else {
      nodes <- data.frame(id = NA)
      return(list(nodes, nodes, p_oncoplot, p_interheatmap))
    }
 
}

interaction_heatmap <- function(res) {

  melted_cormat <- reshape2::melt(res$interactions, na.rm = TRUE) %>%
    mutate(text=ifelse(abs(value)>-log10(0.01),"*",ifelse(abs(value)>-log10(0.05),"��",NA)))
  melted_cormat %>%
    ggplot()+geom_tile(aes(Var1, Var2, fill = value),color="White")+
    geom_text(aes(x=Var1,y=Var2,label=text),color="white",size=8)+
    scale_fill_gradient2(low = "darkgreen", high = "mediumorchid4", mid = "white", 
                         midpoint = 0, limit = c(-3,3), space = "Lab",
                         name="-log10(P-value)",labels=c(">3 (Co-occurrence)","2","1","0","1","2",">3 (Mutually Exclusive)")) +
    guides(fill = guide_colorbar(title.position = "left"))+
    theme_void()+ 
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 0),
          axis.text.y = element_text(angle = 0,size = 10, hjust = 1,vjust=0.5),
          legend.title = element_text(angle = 90,vjust = 1,hjust=1))+
    coord_fixed()
  
}

popfunction <- function(p,position) {
    eval(parse(text=paste0("pushViewport(viewport(",position,",just=c('left','top')))")))
    grid.draw(as.grob(p))
    popViewport()
  }

Coocrrence_network_plot <- function(types,manSelect,gistic=NULL,ks_path,ccf_path,oncoPath,output=NA,ccf_type,p_value=0.05,top_gene,pathway) {
  
  print(types)

  if (nchar(types)==8) typesMaf = c(substr(types,1,4),substr(types,5,8)) else typesMaf = types
  
  ksByCancer = read.csv(ks_path) %>%
    filter(cancertype %in% paste0(sort(typesMaf),collapse  = "")) %>% 
    filter(significant_ks == 1)
  
  #load ccf files
  for (j in 1:length(ccf_path) ) {
    if (j==1) ccf_all= get(load(ccf_path[1])) else 
      ccf_all = rbind(ccf_all,get(load(ccf_path[j])))
  }

  ccf_filter = ccf_all %>% 
    mutate(samplename = substr(Matched_Norm_Sample_Barcode,1,12)) %>%
    filter(samplename %in% manSelect$samplename)
 
  # Load MAF
  if (types != "All") {
    lstObject <- lapply(typesMaf, function(i) tcga_load(i))
  } else {
    lstObject <- lapply(tcga_available()$Study_Abbreviation[1:33], function(i) tcga_load(i))
  }
  
  ### Load selected gene ccf file
  geneCounts <- ksByCancer %>% count(Gene, sort = TRUE, name = "n_genes_across_cancertype")
 
  mergeMaf<- merge_mafs(lstObject)
  mergeMaf@data$Tumor_Sample_Barcode = substr(as.character(mergeMaf@data$Tumor_Sample_Barcode),1,12)

  ## No filter
  mafEsig4 <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig4")$samplename, genes = geneCounts$Gene,
                        query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  mafEsig3 <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig3")$samplename, genes = geneCounts$Gene,
                        query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  
  # Perform interaction analysis using oncogenic pathway data
  NodeEdges4=mafToNodeEdges(maf=mafEsig4,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoPath=oncoPath)
  NodeEdges3=mafToNodeEdges(maf=mafEsig3,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoPath=oncoPath)
  NodeEdges4_pathway=mafToNodeEdges(mafEsig4,pathway=TRUE,top_gene=10,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoPath=oncoPath)
  NodeEdges3_pathway=mafToNodeEdges(mafEsig3,pathway=TRUE,top_gene=10,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoPath=oncoPath)

  # gene network
  if (nrow(NodeEdges4[[2]])>0 & ncol(NodeEdges4[[2]])>1) 
    p1=ggraph_plot(nodes=NodeEdges4[[1]],edges=NodeEdges4[[2]],fill=TRUE,p_value=0.01,title="ESig4") 
  if (nrow(NodeEdges3[[2]])>0 & ncol(NodeEdges3[[2]])>1) 
    p2=ggraph_plot(nodes=NodeEdges3[[1]],edges=NodeEdges3[[2]],fill=TRUE,p_value=0.01,title="ESig3") 

  # oncoplot
  p31=NodeEdges4[[3]];p41=NodeEdges3[[3]]
  p32=NodeEdges4_pathway[[3]];p42=NodeEdges3_pathway[[3]]
  
  # Pathway network
  if (nrow(NodeEdges4_pathway[[2]])>0 & ncol(NodeEdges4_pathway[[2]])>1) 
    p5=ggraph_plot(nodes=NodeEdges4_pathway[[1]],edges=NodeEdges4_pathway[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="Pathway") 

  if (nrow(NodeEdges3_pathway[[2]])>0 & ncol(NodeEdges3_pathway[[2]])>1) 
    p6=ggraph_plot(nodes=NodeEdges3_pathway[[1]],edges=NodeEdges3_pathway[[2]],fill=TRUE,p_value=0.05,title="",subtitle="Pathway") 


  grid.newpage()
  popfunction(p1,position="x = 0, y=1,width = 0.25,height = 0.5")       # ESig4 - Gene Network
  popfunction(p5,position="x = 0.25, y=1,width = 0.25,height = 0.5")    # ESig4 - Pathway Network
  popfunction(p31,position="x = 0.6, y=1,width = 0.4,height = 0.25")    # ESig4 - Gene Oncoplot
  popfunction(p31,position="x = 0.6, y=0.75,width = 0.4,height = 0.25") # ESig4 - Pathway Oncoplot
  popfunction(p2,position="x = 0, y=0.5,width = 0.25,height = 0.5")     # ESig3 - Gene Network
  popfunction(p6,position="x = 0.25, y=0.5,width = 0.25,height = 0.5")  # ESig3 - Pathway Network
  popfunction(p41,position="x = 0.6, y=0.5,width = 0.4,height = 0.25")  # ESig3 - Gene Network
  popfunction(p42,position="x = 0.6, y=0.25,width = 0.4,height = 0.25") # ESig3 - Pathway Network
  
  caption=paste0("Top ",top_gene," genes with significant ccf transition. \n Edge width represents the strength of p value, only edges with p>",p_value," were shown. \n Circle size represent population frequency of each gene.")
  #subtitle = 
  
  g2 <- grid.grab(wrap.grobs = TRUE) %>% as.ggplot() +
    labs(caption=caption)
  g2

  if (!is.na(output)) { ggsave(output,width = 20, height = 15)} 
  
  }
 
Coocrrence_network_plot_MSI <- function(types,manSelect,gistic=NULL,ks_path,ccf_path,oncoPath,output=NA,ccf_type,p_value=0.05,top_gene,pathway,oncoplot.feature,oncoplot.feature.col, bottom_feature) {
  
  print(types)
  
  if (nchar(types)==8) typesMaf = c(substr(types,1,4),substr(types,5,8)) else typesMaf = types
  
  ksByCancer = fread(ks_path)[cancertype %in% paste0(sort(typesMaf),collapse  = "") & (significant_ks == 1|significant_fisher==1)]

  #load ccf files
  for (j in 1:length(ccf_path) ) {
    ccf_load = get(load(ccf_path[j]))
    if (j==1) ccf_all= ccf_load else ccf_all = rbind(ccf_all,ccf_load)
  }
  
  ccf_filter = data.table(ccf_all)[,samplename := substr(Matched_Norm_Sample_Barcode,1,12)][samplename %in% manSelect$samplename]
    
  # Load MAF
  if (types != "All") {
    lstObject <- lapply(typesMaf, function(i) tcga_load(i))
  } else {
    lstObject <- lapply(tcga_available()$Study_Abbreviation[1:33], function(i) tcga_load(i))
  }
  
  ### Load selected gene ccf file
  geneCounts <- ksByCancer %>% count(Gene, sort = TRUE, name = "n_genes_across_cancertype") %>% rbind(list("POLE",1))

  mergeMaf<- merge_mafs(lstObject)
  mergeMaf@data$Tumor_Sample_Barcode = substr(as.character(mergeMaf@data$Tumor_Sample_Barcode),1,12)
  
  ## No filter
  
  mafEsig4 <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig4")$samplename, genes = geneCounts$Gene,
                        query = "Variant_Classification %nin%  c('Translation_Start_Site')")

  if (nrow(subset(manSelect,subtype=="ESig4" & Subtype=="MSS"))>0) {                      
    mafEsig4MSS <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig4" & Subtype=="MSS")$samplename, genes = geneCounts$Gene,
                          query = "Variant_Classification %nin%  c('Translation_Start_Site')")
    MSS4 = TRUE} else {MSS4 = FALSE}

  if (nrow(subset(manSelect,subtype=="ESig4" & Subtype=="MMR"))>0) {
    mafEsig4MMR <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig4" & Subtype=="MMR")$samplename, genes = geneCounts$Gene,
                            query = "Variant_Classification %nin%  c('Translation_Start_Site')")
    MMR4 = TRUE } else {MMR4 = FALSE}

  if (nrow(subset(manSelect,subtype=="ESig4" & Subtype=="POLE"))>0) {
    mafEsig4POLE <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig4" & Subtype=="POLE")$samplename, genes = geneCounts$Gene,
                            query = "Variant_Classification %nin%  c('Translation_Start_Site')")
    POLE4 = TRUE } else {POLE4 = FALSE}
  
  mafEsig3 <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig3")$samplename, genes = geneCounts$Gene,
                        query = "Variant_Classification %nin%  c('Translation_Start_Site')")

  if (nrow(subset(manSelect,subtype=="ESig3" & Subtype=="MSS"))>0) {         
  mafEsig3MSS <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig3" & Subtype=="MSS")$samplename, genes = geneCounts$Gene,
                           query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  MSS3 = TRUE} else {MSS3 = FALSE}

  if (nrow(subset(manSelect,subtype=="ESig3" & Subtype=="MMR"))>0) {
  mafEsig3MMR <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig3" & Subtype=="MMR")$samplename, genes = geneCounts$Gene,
                           query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  MMR3 = TRUE } else {MMR3 = FALSE}   

  if (nrow(subset(manSelect,subtype=="ESig3" & Subtype=="POLE"))>0) {                 
  mafEsig3POLE <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig3" & Subtype=="POLE")$samplename, genes = geneCounts$Gene,
                            query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  POLE3 = TRUE } else {POLE3 = FALSE}                        

  gistic = subset(gistic,Sample %in% geneCounts$Gene)
  # Perform interaction analysis using oncogenic pathway data
  NodeEdges4=mafToNodeEdges(maf=mafEsig4,top_gene=30,manSelect=manSelect,gistic=gistic,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
  NodeEdges3=mafToNodeEdges(maf=mafEsig3,top_gene=30,manSelect=manSelect,gistic=gistic,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
  NodeEdges4_pathway=mafToNodeEdges(mafEsig4,pathway=TRUE,manSelect=manSelect,top_gene=10,gistic=gistic,ksByCancer=ksByCancer,oncoplot.feature=oncoplot.feature,ccf_filter=ccf_filter,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
  NodeEdges3_pathway=mafToNodeEdges(mafEsig3,pathway=TRUE,manSelect=manSelect,top_gene=10,gistic=gistic,ksByCancer=ksByCancer,oncoplot.feature=oncoplot.feature,ccf_filter=ccf_filter,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
  
  if (MSS4) {
    NodeEdges4MSS =mafToNodeEdges(maf=mafEsig4MSS,manSelect=manSelect,top_gene=30,gistic=gistic,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
    if (nrow(NodeEdges4MSS[[2]])>0 & ncol(NodeEdges4MSS[[2]])>1) p11=ggraph_plot(nodes=NodeEdges4MSS[[1]],edges=NodeEdges4MSS[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MSS") 
      else MSS4=FALSE
  }
  if (MMR4) {
    NodeEdges4MMR =mafToNodeEdges(maf=mafEsig4MMR,manSelect=manSelect,top_gene=30,gistic=gistic,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
    if (nrow(NodeEdges4MMR[[2]])>0 & ncol(NodeEdges4MMR[[2]])>1 & nrow(NodeEdges4MMR[[1]])>0) p12=ggraph_plot(nodes=NodeEdges4MMR[[1]],edges=NodeEdges4MMR[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MMR") 
      else MMR4 =FALSE
  }

  if (POLE4) {
    NodeEdges4POLE=mafToNodeEdges(mafEsig4POLE,manSelect=manSelect,top_gene=30,gistic=gistic,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
     if (nrow(NodeEdges4POLE[[2]])>0 & ncol(NodeEdges4POLE[[2]])>1) p13=ggraph_plot(nodes=NodeEdges4POLE[[1]],edges=NodeEdges4POLE[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="POLE") 
      else POLE4 =FALSE
  }
  if (MSS3) {
    NodeEdges3MSS =mafToNodeEdges(mafEsig3MSS,manSelect=manSelect,top_gene=30,gistic=gistic,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
    if (nrow(NodeEdges3MSS[[2]])>0 & ncol(NodeEdges3MSS[[2]])>1) p21=ggraph_plot(nodes=NodeEdges3MSS[[1]],edges=NodeEdges3MSS[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MSS") 
     else MSS3 =FALSE
  }

  if (MMR3) {
    NodeEdges3MMR =mafToNodeEdges(maf=mafEsig3MMR,manSelect=manSelect,top_gene=30,gistic=gistic,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
    if (nrow(NodeEdges3MMR[[2]])>0 & ncol(NodeEdges3MMR[[2]])>1) p22=ggraph_plot(nodes=NodeEdges3MMR[[1]],edges=NodeEdges3MMR[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MMR") 
      else MMR3 =FALSE
  }
  if (POLE3) {
    NodeEdges3POLE=mafToNodeEdges(maf=mafEsig3POLE,manSelect=manSelect,top_gene=30,gistic=gistic,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
    if (nrow(NodeEdges3POLE[[2]])>0 & ncol(NodeEdges3POLE[[2]])>1) p23=ggraph_plot(nodes=NodeEdges3POLE[[1]],edges=NodeEdges3POLE[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="POLE") 
      else POLE3 =FALSE
  }
  # gene network
  if (nrow(NodeEdges4[[2]])>0 & ncol(NodeEdges4[[2]])>1) p1=ggraph_plot(nodes=NodeEdges4[[1]],edges=NodeEdges4[[2]],fill=TRUE,p_value=p_value,title="ESig4") 
  if (nrow(NodeEdges3[[2]])>0 & ncol(NodeEdges3[[2]])>1) p2=ggraph_plot(nodes=NodeEdges3[[1]],edges=NodeEdges3[[2]],fill=TRUE,p_value=p_value,title="ESig3") 

  # oncoplot
  p31=NodeEdges4[[3]];p41=NodeEdges3[[3]]
  p32=NodeEdges4_pathway[[3]];p42=NodeEdges3_pathway[[3]]
  
  # Pathway network
  if (nrow(NodeEdges4_pathway[[2]])>0 & ncol(NodeEdges4_pathway[[2]])>1) {
    p5=ggraph_plot(nodes=NodeEdges4_pathway[[1]],edges=NodeEdges4_pathway[[2]],fill=TRUE,p_value=0.01,title="Pathway") 
    pathway4 = TRUE } else {pathway4 = FALSE}
  if (nrow(NodeEdges3_pathway[[2]])>0 & ncol(NodeEdges3_pathway[[2]])>1) {
    p6=ggraph_plot(nodes=NodeEdges3_pathway[[1]],edges=NodeEdges3_pathway[[2]],fill=TRUE,p_value=0.01,title="Pathway")
    pathway3 = TRUE
  } else {pathway3 = FALSE}
  
 
  # ESig3 - Gene Network
  grid.newpage()
  popfunction(p1,position="x = 0, y=1,width = 0.25,height = 0.25")       # ESig4 - Gene Network
  if (MSS4) popfunction(p11,position="x = 0, y=0.75,width = 0.2,height = 0.2")    # ESig4 - Subtype1 Gene Network
  if (MMR4) popfunction(p12,position="x = .2, y=0.75,width = 0.2,height = 0.2")   # ESig4 - Subtype2 Gene Network
  if (POLE4) popfunction(p13,position="x = .4, y=0.75,width = 0.2,height = 0.2")   # ESig4 - Subtype3 Gene Network
  if (pathway4) popfunction(p5,position="x = 0.25, y=1,width = 0.3,height = 0.25")     # ESig4 - Pathway Network
  popfunction(p31,position="x = 0.5, y=1,width = 0.5,height = 0.25")    # ESig4 - Gene Oncoplot
  popfunction(p32,position="x = 0.5, y=0.75,width = 0.5,height = 0.25") # ESig4 - Pathway Oncoplot
  popfunction(p2,position="x = 0, y=0.5,width = 0.25,height = 0.25")     # ESig3 - Gene Network
  if (MSS3) popfunction(p21,position="x = 0, y=0.25,width = 0.2,height = 0.2")    # ESig3 - Subtype1 Gene Network
  if (MMR3) popfunction(p22,position="x = .2, y=0.25,width = 0.2,height = 0.2")   # ESig3 - Subtype2 Gene Network
  if (POLE3) popfunction(p23,position="x = .4, y=0.25,width = 0.2,height = 0.2")   # ESig3 - Subtype3 Gene Network
  if (pathway3) popfunction(p6,position="x = 0.25, y=0.5,width = 0.3,height = 0.25")   # ESig3 - Pathway Network
  popfunction(p41,position="x = 0.5, y=0.5,width = 0.5,height = 0.25")  # ESig3 - Gene Oncoplot
  popfunction(p42,position="x = 0.5, y=0.25,width = 0.5,height = 0.25") # ESig3 - Pathway Oncoplot

  caption=paste0("Top ",top_gene," genes with significant ccf transition. \n Edge width represents the strength of p value, only edges with p<",p_value," were shown. \n Circle size represent population frequency of each gene.")
  #subtitle = 
  
  g2 <- grid.grab(wrap.grobs = TRUE) %>% as.ggplot() +
    labs(caption=caption)
  g2
  
  if (!is.na(output)) { ggsave(output,width = 30, height = 30,device = cairo_pdf)} 
}

    

