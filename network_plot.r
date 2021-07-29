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
  pathway_col = data.frame(Pathway=c("Cell_Cycle", "Chromatin_modifiers", "Notch", "NRF2", "PI3K", "RTK_RAS", "TGF-Beta", "TP53", "WNT", "Hippo", "Myc",NA), col=c(mypal[c(1:8,10,15:16)],"grey70"))
  pathway_col_values=pathway_col$col
  names(pathway_col_values) = pathway_col$Pathway
  
  gene_edges = edges %>% filter(gene1 %in% nodes$label,gene2 %in% nodes$label) 
  from_id = sapply( edges$gene1,function(x) nodes %>% filter(label==x) %>% .$id ) %>% as.numeric()
  to_id = sapply( edges$gene2,function(x) nodes %>% filter(label==x) %>% .$id ) %>% as.numeric()
  
  ### get layout
  edges <- data.frame(from = from_id, to = to_id,
                      color = ifelse( gene_edges$Event =="Co_Occurence", "#4058a3",  "#f18826")) %>% cbind(gene_edges) 
  
  co_nodes = c(t(edges[which(edges$Event=="Co_Occurence"),c("gene1","gene2")]))  %>% unique(.)
  mu_nodes = c(t(edges[which(edges$Event=="Mutually_Exclusive"),c("gene1","gene2")]))  %>% unique(.)
  
  #nodes_mu = nodes %>% filter(label %in% mu_nodes) 
  if (length(co_nodes)>0) {
   nodes_co = nodes %>% filter(label %in% co_nodes) 
   edges_co = edges %>% filter(Event=="Co_Occurence",gene1 %in% co_nodes, gene2 %in% co_nodes)
   nodes_mu = nodes %>% filter(label %in% mu_nodes) 
   edges_mu = edges %>% filter(Event=="Mutually_Exclusive",gene1 %in% mu_nodes , gene2 %in% mu_nodes)
   net_co = graph_from_data_frame(d=edges_co, vertices=nodes_co, directed=F) 
  # Merge group of co-occurrence nodes label to whole nodes 
  } else {
    nodes_mu = nodes %>% filter(label %in% mu_nodes) 
    edges_mu = edges %>% filter(Event=="Mutually_Exclusive",gene1 %in% mu_nodes , gene2 %in% mu_nodes)
    net_co = graph_from_data_frame(d=edges_mu, vertices=nodes_mu, directed=F) 
  }
  # Calculate groups for nodes with co-occurrences relationship
  tg_co <- tidygraph::as_tbl_graph(net_co) %>% activate(nodes) %>% mutate(label=label,name=as.numeric(name)) %>% mutate(community_id = group_components()) 
  gene_nodes = left_join(nodes,data.frame(label=V(tg_co)$label,community_id=V(tg_co)$community_id),by="label") %>%
    mutate(community_id=as.factor(ifelse(is.na(community_id),0,community_id)))
  net= graph_from_data_frame(d=edges,vertices=gene_nodes, directed=F)  
  
  tg <- tidygraph::as_tbl_graph(net) %>% activate(nodes) %>% mutate(label=label)
  if (no_layout) l = net_to_layout(net,gene_nodes,simple=TRUE) else
    l = net_to_layout(net,gene_nodes)
  
  v.size <- V(tg)$freq_value %>% as.character() %>% as.numeric() %>% rescale(c(15,30))
  e.color <- E(tg)$color
  E(tg)$weight <- E(tg)$log10pval %>% rescale(c(0.5,2.5))
  eigenCent <- evcent(tg)$vector
  bins <- unique(quantile(eigenCent, seq(0,1,length.out=30)))

  
  p = tg %>%
    ggraph(layout=l)
    
  if (nrow(nodes_co)>0) p = p+geom_edge_link0(aes(filter=(Event=="Co_Occurence")),edge_colour="#4058a3",width=E(tg)$weight[which(E(tg)$Event=="Co_Occurence")]) 
  if (nrow(nodes_mu)>0) p = p+ geom_edge_link0(aes(filter=(Event=="Mutually_Exclusive")),edge_colour="#f18826",width=E(tg)$weight[which(E(tg)$Event=="Mutually_Exclusive")])
  
  p = p +
    geom_node_point(aes(color=group),size=v.size/2,alpha = 0.9) +
    ggtitle(title)+
    scale_color_manual(values=pathway_col_values)+
    #geom_node_point(aes(size=as.numeric(as.character(freq_value))),alpha = 0.9) +
    #guides(shape = guide_legend(override.aes = list(size = 15:23)))+
    geom_node_text(aes(label = label), #repel = TRUE,
                   size=log(v.size),colour='black')+
    theme_classic()+
    scale_x_reverse()+
    labs(x=ccf_type,y="",subtitle=subtitle)+
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y =element_blank(),
          legend.position = "none")
  
  if (fill) {p = p+   
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
mafToNodeEdges <- function(maf,ccf_filter,ksByCancer,pathway=FALSE,top_gene,p_value=0.05,oncoplot.feature=NA,fabcolors) {
  
  if (pathway) {
    mafTemp   <- maf@data
    mafTemp <- inner_join(mafTemp, oncoPath, by = c("Hugo_Symbol"="Gene"))
    colnames(mafTemp)[1] <- "Hugo_Symbol_Gene"
    colnames(mafTemp)[which(colnames(mafTemp)=="Pathway")] <- "Hugo_Symbol"
    
    clinicalTemp <- getSampleSummary(x = maf) 
    clinicalTemp <- left_join(clinicalTemp, manSelect, by = c("Tumor_Sample_Barcode"="samplename")) %>%
      distinct(Tumor_Sample_Barcode,.keep_all = TRUE)
    mafClinical = read.maf(maf = mafTemp, clinicalData = clinicalTemp)
    
    annotationDat = clinicalTemp %>%
      mutate(Tumor_Sample_Barcode = levels(mafClinical@data[["Tumor_Sample_Barcode"]]))

    p_oncoplot = as.ggplot(function()
      if (all(!is.na(oncoplot.feature))) {
        oncoplot(maf = mafClinical, top=top_gene, removeNonMutated = FALSE, #colors = mycolors, 
                 drawRowBar = TRUE,
                 annotationColor=fabcolors,
                 sortByAnnotation = TRUE, keepGeneOrder = FALSE, bgCol="#F5F5F5",
                 clinicalFeatures=oncoplot.feature)
      } else {
        oncoplot(maf = mafClinical, top=top_gene, removeNonMutated = FALSE, #colors = mycolors, 
                 drawRowBar = TRUE,
                 sortByAnnotation = TRUE, keepGeneOrder = FALSE, bgCol="#F5F5F5")
      }
    )
    
    res <- getInteractions(mafClinical, top=10, colPal = "PRGn", returnAll = TRUE, sigSymbolsSize=1.5)
    p_interheatmap = interaction_heatmap(res)
    gene_edges = res$pairs %>% filter(pValue<p_value) 
    gene_nodes = cbind(res$pairs$gene1,res$pairs$gene2) %>% t(.) %>% as.character() %>%unique(.)
    gene_nodes_col = data.frame(Pathway=gene_nodes) %>%  
      left_join( pathway_col,by="Pathway")
    
    geneccfFreq <- ccf_filter %>%
      mutate(ccube_ccf= ifelse(ccube_ccf>=1,1,ccube_ccf))%>%
      filter(Variant_Classification %nin%  c('Translation_Start_Site'),
             Hugo_Symbol %in% ksByCancer$Gene) %>%
      inner_join(oncoPath, by=c("Hugo_Symbol"="Gene")) %>%
      group_by(Pathway) %>%
      dplyr::summarise(mean_ccf=mean(ccube_ccf,na.rm=TRUE),median_ccf=median(ccube_ccf,na.rm=TRUE),n=length(unique(samplename))) %>% 
      mutate(freq=n/nrow(manSelect)) %>%
      mutate(freq_value=cut(.$freq,breaks=c(seq(0,1,length.out=10)),labels=seq(15,50,length.out=9)))
    colnames(geneccfFreq)[1] = "label"
    
  } else {
    mafTemp   <- maf@data
    clinicalTemp <- getSampleSummary(x = maf) 
    clinicalTemp <- left_join(clinicalTemp, manSelect, by = c("Tumor_Sample_Barcode"="samplename")) %>%
      distinct(Tumor_Sample_Barcode,.keep_all = TRUE)
    mafClinical = read.maf(maf = mafTemp, clinicalData = clinicalTemp)
    
    # neoantigen = clinicalTemp %>% select(Tumor_Sample_Barcode,Neoantigen_mutation.sc) %>%
    #   distinct(Tumor_Sample_Barcode,.keep_all = TRUE) %>%
    #   as.data.table()
    # neoantigen[,Neoantigen_mutation.sc:=ifelse(is.na(Neoantigen_mutation.sc),0,Neoantigen_mutation.sc)]  
  
    p_oncoplot = as.ggplot(function()
      if (all(!is.na(oncoplot.feature))) {
      oncoplot(maf = mafClinical, top=top_gene, removeNonMutated = FALSE, #colors = mycolors, 
               drawRowBar = TRUE,
               annotationColor=fabcolors,
               sortByAnnotation = TRUE, keepGeneOrder = FALSE, bgCol="#F5F5F5",
               clinicalFeatures=oncoplot.feature)
      } else {
      oncoplot(maf = mafClinical, top=top_gene, removeNonMutated = FALSE, #colors = #mycolors, 
               drawRowBar = TRUE,
               sortByAnnotation = TRUE, keepGeneOrder = FALSE, bgCol="#F5F5F5")
      }
    )
    
    res <- getInteractions(maf, top=top_gene, colPal = "PRGn", returnAll = TRUE, sigSymbolsSize=1.5) 
    p_interheatmap = interaction_heatmap(res)
    gene_edges = res$pairs %>% filter(pValue<p_value) 
    gene_nodes = cbind(gene_edges$gene1,gene_edges$gene2) %>% t(.) %>% as.character() %>%unique(.)
    gene_nodes_col = data.frame(Gene=gene_nodes) %>%  
      left_join(oncoPath, by="Gene") %>%
      left_join( pathway_col,by="Pathway") %>%
      distinct(Gene,.keep_all = TRUE)
   
    geneccfFreq <- ccf_filter %>%
      mutate(ccube_ccf= ifelse(ccube_ccf>=1,1,ccube_ccf))%>%
      filter(Variant_Classification %nin%  c('Translation_Start_Site'),
             Hugo_Symbol %in% ksByCancer$Gene,
             samplename %in% clinicalTemp$Tumor_Sample_Barcode) %>%
      group_by(Hugo_Symbol) %>%
      dplyr::summarise(mean_ccf=mean(ccube_ccf,na.rm=TRUE),median_ccf=median(ccube_ccf,na.rm=TRUE),n=length(unique(samplename))) %>% 
      mutate(freq=n/nrow(manSelect)) %>%
      mutate(freq_value=cut(.$freq,breaks=c(seq(0,1,length.out=10)),labels=seq(15,31,by=2)))
    colnames(geneccfFreq)[1] = "label"
  }
  
  if (length(gene_nodes)>0) {
    
    nodes <- data.frame(id = 1:length(gene_nodes),
                        label = gene_nodes,
                        group = gene_nodes_col$Pathway,
                        color = gene_nodes_col$col,
                        shadow = TRUE) %>%
      left_join(geneccfFreq,by="label") 
    
    
    nodes = nodes[which(nodes[,ccf_type]>0),] %>%
      mutate(value=cut(.[,ccf_type],breaks=(1:10)/10,labels=1:9))
    
    nodes[,ccf_type] =round(nodes[,ccf_type],3)
    return(list(nodes,gene_edges,p_oncoplot,p_interheatmap))
  } else {
    nodes <- data.frame(id = NA)
    return(list(nodes,nodes,p_oncoplot,p_interheatmap))
  }
 
}
interaction_heatmap <- function(res) {
  library(reshape2)
  melted_cormat <- reshape2::melt(res$interactions, na.rm = TRUE) %>%
    mutate(text=ifelse(abs(value)>-log10(0.01),"*",ifelse(abs(value)>-log10(0.05),"¡¤",NA)))
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

Coocrrence_network_plot <- function(types,manSelect,ks_path,ccf_path,oncoPath,output=NA,ccf_type,oncoplot.feature.col,p_value=0.05,top_gene,pathway,oncoplot.feature) {
  
 
  print(types)

  if (nchar(types)==8) typesMaf = c(substr(types,1,4),substr(types,5,8)) else typesMaf = types
  
  ksByCancer = read.csv(ks_path) %>%
    filter(cancertype %in% paste0(sort(typesMaf),collapse  = "")) %>% 
    filter(significant_ks == 1)
  
  #load ccf files
  for (j in 1:length(ccf_path) ) {
    if (j==1) ccf_all= get(load(paste0(folder_path,"Figures and Data/CCF/",ccf_path[1]))) else 
      ccf_all = rbind(ccf_all,get(load(paste0(folder_path,"Figures and Data/CCF/",ccf_path[j]))))
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
  NodeEdges4=mafToNodeEdges(maf=mafEsig4,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  NodeEdges3=mafToNodeEdges(maf=mafEsig3,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  NodeEdges4_pathway=mafToNodeEdges(mafEsig4,pathway=TRUE,top_gene=10,p_value=0.05,ksByCancer=ksByCancer,oncoplot.feature=oncoplot.feature,ccf_filter=ccf_filter,fabcolors = oncoplot.feature.col)
  NodeEdges3_pathway=mafToNodeEdges(mafEsig3,pathway=TRUE,top_gene=10,p_value=0.05,ksByCancer=ksByCancer,oncoplot.feature=oncoplot.feature,ccf_filter=ccf_filter,fabcolors = oncoplot.feature.col)

  # gene network
  if (nrow(NodeEdges4[[2]])>0 & ncol(NodeEdges4[[2]])>1) 
  p1=ggraph_plot(nodes=NodeEdges4[[1]],edges=NodeEdges4[[2]],fill=TRUE,p_value=0.01,title="ESig4") 
  if (nrow(NodeEdges3[[2]])>0 & ncol(NodeEdges3[[2]])>1) 
  p2=ggraph_plot(nodes=NodeEdges3[[1]],edges=NodeEdges3[[2]],fill=TRUE,p_value=0.01,title="ESig3") 
  
  # oncoplot
  p3=NodeEdges4[[3]]
  p4=NodeEdges3[[3]]
  # Pathway network
  if (nrow(NodeEdges4_pathway[[2]])>0 & ncol(NodeEdges4_pathway[[2]])>1) p5=ggraph_plot(nodes=NodeEdges4_pathway[[1]],edges=NodeEdges4_pathway[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="Pathway") 
  if (nrow(NodeEdges3_pathway[[2]])>0 & ncol(NodeEdges3_pathway[[2]])>1) p6=ggraph_plot(nodes=NodeEdges3_pathway[[1]],edges=NodeEdges3_pathway[[2]],fill=TRUE,p_value=0.05,title="",subtitle="Pathway") 
  
  grid.newpage()
  pushViewport(viewport(x = 0, y=1,width = 0.25,height = 0.5,just = c("left", "top")))
  grid.draw(as.grob(p1))
  popViewport()
  

  pushViewport(viewport(x = 0.25, y=1,width = 0.25,height = 0.5,just = c("left", "top")))
  grid.draw(as.grob(p5))
  popViewport()
  
  pushViewport(viewport(x = 0.5, y=1,width = 0.5,height = 0.5,just = c("left", "top")))
  grid.draw(as.grob(p3))
  popViewport()
  
 
  pushViewport(viewport(x = 0, y=0.5,width = 0.25,height = 0.5,just = c("left", "top")))
  grid.draw(as.grob(p2))
  popViewport()
  # vp = viewport(x=.8, y=.75, width=.3, height=.4)
  # pushViewport(vp)
  # grid.draw(as.grob(p5))
  # upViewport()
  pushViewport(viewport(x = 0.25, y=0.5,width = 0.25,height = 0.5,just = c("left", "top")))
  grid.draw(as.grob(p6))
  popViewport()
  
  pushViewport(viewport(x = 0.5, y=0.5,width = 0.5,height = 0.5,just = c("left", "top")))
  grid.draw(as.grob(p4))
  popViewport()
  
 
  caption=paste0("Top ",top_gene," genes with significant ccf transition. \n Edge width represents the strength of p value, only edges with p>",p_value," were shown. \n Circle size represent population frequency of each gene.")
  #subtitle = 
  
  g2 <- grid.grab(wrap.grobs = TRUE) %>% as.ggplot() +
    labs(caption=caption)
  g2

  if (!is.na(output)) { ggsave(output,width = 20, height = 15)} 
  
  }
 
Coocrrence_network_plot_MSI <- function(types,manSelect,ks_path,ccf_path,oncoPath,output=NA,ccf_type,p_value=0.05,top_gene,pathway,oncoplot.feature,oncoplot.feature.col) {
  
  print(types)
  
  if (nchar(types)==8) typesMaf = c(substr(types,1,4),substr(types,5,8)) else typesMaf = types
  
  ksByCancer = read.csv(ks_path) %>%
    filter(cancertype %in% paste0(sort(typesMaf),collapse  = "")) %>% 
    filter(significant_ks == 1)
  
  #load ccf files
  for (j in 1:length(ccf_path) ) {
    if (j==1) ccf_all= get(load(paste0(folder_path,"Figures and Data/CCF/",ccf_path[1]))) else 
      ccf_all = rbind(ccf_all,get(load(paste0(folder_path,"Figures and Data/CCF/",ccf_path[j]))))
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
  mafEsig4MSS <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig4" & Subtype=="MSS")$samplename, genes = geneCounts$Gene,
                        query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  mafEsig4MMR <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig4" & Subtype=="MMR")$samplename, genes = geneCounts$Gene,
                           query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  mafEsig4POLE <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig4" & Subtype=="POLE")$samplename, genes = geneCounts$Gene,
                           query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  
  mafEsig3 <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig3")$samplename, genes = geneCounts$Gene,
                        query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  mafEsig3MSS <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig3" & Subtype=="MSS")$samplename, genes = geneCounts$Gene,
                           query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  mafEsig3MMR <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig3" & Subtype=="MMR")$samplename, genes = geneCounts$Gene,
                           query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  mafEsig3POLE <- subsetMaf(mergeMaf, tsb=subset(manSelect,subtype=="ESig3" & Subtype=="POLE")$samplename, genes = geneCounts$Gene,
                            query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  
  # Perform interaction analysis using oncogenic pathway data
  NodeEdges4=mafToNodeEdges(maf=mafEsig4,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  NodeEdges3=mafToNodeEdges(maf=mafEsig3,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  NodeEdges4_pathway=mafToNodeEdges(mafEsig4,pathway=TRUE,top_gene=10,p_value=0.05,ksByCancer=ksByCancer,oncoplot.feature=oncoplot.feature,ccf_filter=ccf_filter,fabcolors = oncoplot.feature.col)
  NodeEdges3_pathway=mafToNodeEdges(mafEsig3,pathway=TRUE,top_gene=10,p_value=0.05,ksByCancer=ksByCancer,oncoplot.feature=oncoplot.feature,ccf_filter=ccf_filter,fabcolors = oncoplot.feature.col)
  
  NodeEdges4MSS=mafToNodeEdges(mafEsig4MSS,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  NodeEdges4MMR=mafToNodeEdges(mafEsig4MMR,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  NodeEdges4POLE=mafToNodeEdges(mafEsig4POLE,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  NodeEdges3MSS=mafToNodeEdges(mafEsig3MSS,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  NodeEdges3MMR=mafToNodeEdges(mafEsig3MMR,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  NodeEdges3POLE=mafToNodeEdges(maf=mafEsig3POLE,top_gene=30,p_value=0.05,ksByCancer=ksByCancer,ccf_filter=ccf_filter,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col)
  
  # gene network
  if (nrow(NodeEdges4[[2]])>0 & ncol(NodeEdges4[[2]])>1) p1=ggraph_plot(nodes=NodeEdges4[[1]],edges=NodeEdges4[[2]],fill=TRUE,p_value=0.01,title="ESig4") 
  if (nrow(NodeEdges3[[2]])>0 & ncol(NodeEdges3[[2]])>1) p2=ggraph_plot(nodes=NodeEdges3[[1]],edges=NodeEdges3[[2]],fill=TRUE,p_value=0.01,title="ESig3") 

  # oncoplot
  p31=NodeEdges4[[3]];p41=NodeEdges3[[3]]
  p32=NodeEdges4_pathway[[3]];p42=NodeEdges3_pathway[[3]]
  
  # Pathway network
  if (nrow(NodeEdges4_pathway[[2]])>0 & ncol(NodeEdges4_pathway[[2]])>1) p5=ggraph_plot(nodes=NodeEdges4_pathway[[1]],edges=NodeEdges4_pathway[[2]],fill=TRUE,p_value=0.01,title="Pathway") 
  if (nrow(NodeEdges3_pathway[[2]])>0 & ncol(NodeEdges3_pathway[[2]])>1) p6=ggraph_plot(nodes=NodeEdges3_pathway[[1]],edges=NodeEdges3_pathway[[2]],fill=TRUE,p_value=0.01,title="Pathway")
  
  # Subtype network
  if (nrow(NodeEdges4MSS[[2]])>0 & ncol(NodeEdges4MSS[[2]])>1) p11=ggraph_plot(nodes=NodeEdges4MSS[[1]],edges=NodeEdges4MSS[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MSS") 
  if (nrow(NodeEdges4MMR[[2]])>0 & ncol(NodeEdges4MMR[[2]])>1) p12=ggraph_plot(nodes=NodeEdges4MMR[[1]],edges=NodeEdges4MMR[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MMR") 
  if (nrow(NodeEdges4POLE[[2]])>0 & ncol(NodeEdges4POLE[[2]])>1) p13=ggraph_plot(nodes=NodeEdges4POLE[[1]],edges=NodeEdges4POLE[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="POLE") 
  
  if (nrow(NodeEdges3MSS[[2]])>0 & ncol(NodeEdges3MSS[[2]])>1) p21=ggraph_plot(nodes=NodeEdges3MSS[[1]],edges=NodeEdges3MSS[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MSS") 
  if (nrow(NodeEdges3MMR[[2]])>0 & ncol(NodeEdges3MMR[[2]])>1) p22=ggraph_plot(nodes=NodeEdges3MMR[[1]],edges=NodeEdges3MMR[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MMR") 
  if (nrow(NodeEdges3POLE[[2]])>0 & ncol(NodeEdges3POLE[[2]])>1) p23=ggraph_plot(nodes=NodeEdges3POLE[[1]],edges=NodeEdges3POLE[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="POLE") 
  
  # ESig3 - Gene Network
  grid.newpage()
  pushViewport(viewport(x = 0, y=1,width = 0.3,height = 0.25,just = c("left", "top")))
  grid.draw(as.grob(p1))
  popViewport()
  
  pushViewport(viewport(x = 0, y=0.75,width = 0.2,height = 0.2,just = c("left", "top")))
  grid.draw(as.grob(p11))
  popViewport()
  
  pushViewport(viewport(x = .2, y=0.75,width = 0.2,height = 0.2,just = c("left", "top")))
  grid.draw(as.grob(p12))
  popViewport()
  
  pushViewport(viewport(x = .4, y=0.75,width = 0.2,height = 0.2,just = c("left", "top")))
  grid.draw(as.grob(p13))
  popViewport()
  
  # ESig3 - Pathway Network
  pushViewport(viewport(x = 0.3, y=1,width = 0.3,height = 0.25,just = c("left", "top")))
  grid.draw(as.grob(p5))
  popViewport()
  
  # ESig3 - Oncoplot
  pushViewport(viewport(x = 0.6, y=1,width = 0.4,height = 0.25,just = c("left", "top")))
  grid.draw(as.grob(p31))
  popViewport()
  
  pushViewport(viewport(x = 0.6, y=0.75,width = 0.4,height = 0.25,just = c("left", "top")))
  grid.draw(as.grob(p32))
  popViewport()
  
  # ESig4 - Gene Network
  pushViewport(viewport(x = 0, y=0.5,width = 0.3,height = 0.25,just = c("left", "top")))
  grid.draw(as.grob(p2))
  popViewport()
  
  pushViewport(viewport(x = 0, y=0.25,width = 0.2,height = 0.2,just = c("left", "top")))
  grid.draw(as.grob(p21))
  popViewport()
  
  pushViewport(viewport(x = .2, y=0.25,width = 0.2,height = 0.2,just = c("left", "top")))
  grid.draw(as.grob(p22))
  popViewport()
  
  pushViewport(viewport(x = .4, y=0.25,width = 0.2,height = 0.2,just = c("left", "top")))
  grid.draw(as.grob(p23))
  popViewport()
  
  # ESig4 - Pathway Network
  pushViewport(viewport(x = 0.3, y=0.5,width = 0.3,height = 0.25,just = c("left", "top")))
  grid.draw(as.grob(p6))
  popViewport()
  
  # ESig4 - Ocoplot
  pushViewport(viewport(x = 0.6, y=0.5,width = 0.4,height = 0.25,just = c("left", "top")))
  grid.draw(as.grob(p41))
  popViewport()
  
  pushViewport(viewport(x = 0.6, y=0.25,width = 0.4,height = 0.25,just = c("left", "top")))
  grid.draw(as.grob(p42))
  popViewport()
  
  
  caption=paste0("Top ",top_gene," genes with significant ccf transition. \n Edge width represents the strength of p value, only edges with p>",p_value," were shown. \n Circle size represent population frequency of each gene.")
  #subtitle = 
  
  g2 <- grid.grab(wrap.grobs = TRUE) %>% as.ggplot() +
    labs(caption=caption)
  g2
  
  if (!is.na(output)) { ggsave(output,width = 20, height = 15)} 
}

  
  #legend
    # data.frame(group=pathway_col)
    # pathway_col_values=pathway_col$col
    # names(pathway_col_values) = pathway_col$Pathway
    # pathway = pathway_col %>% ggplot() + geom_point(aes(x=Pathway,y=1,fill=Pathway,color=Pathway),size=9) + 
    #   scale_color_manual(values= values)+
    #   scale_fill_manual(values= values)+
    #   theme_void()
    # legend_pathway = extractLegend(pathway)
    # 
    # data.frame(value = (1:9)/10,fre_labels=seq(15,31,by=2)) %>%
    #   ggplot() + geom_point(aes(x=value,y=1,cex=fre_labels)) 
    # 
    # scales::show_col(pathway_col$col)
    # p2+theme(legend.position = "none")
    # 
      
    

