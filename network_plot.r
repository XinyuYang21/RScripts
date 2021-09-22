# Multi-layer network
mypal <- colorRampPalette(pal_nejm("default", alpha = 1)(8))(18)
pathway_col = data.table(Pathway=c("Cell_Cycle", "Chromatin_modifiers", "Notch", "NRF2", "PI3K", "RTK_RAS", "TGF-Beta", "TP53", "WNT", "Hippo", "Myc",NA), col=c(mypal[1:11],"grey70"))


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
  df %>% mutate(x_level=as.numeric(node$ccf_type)) 
}

layout_to_matrix <- function(layout) {layout %>% arrange(id) %>% .[,1:2] %>% as.matrix()}

net_to_layout <- function(net,node,simple=FALSE){
  if (is.data.table(node)) node = as.data.frame(node)
  
  layout = cbind(get_layers(net),layout_with_sugiyama(net, layers=V(net)$value)$layout[,1]) %>%
    as.data.frame() %>% cbind(node[,c("id","community_id","freq","ccf_type")]) %>%
    mutate(community_id=as.numeric(as.character(community_id))) 
  
  check_layout_facet(layout=layout,facet="community_id",node=node,simple=simple)  %>%
    layout_to_matrix()
}


ggraph_plot = function(nodes,edges,facet=FALSE,title="",fill=FALSE,zoom=NA,p_value=0.01,subtype_var_col=NA,no_layout=FALSE,subtitle="") {
  
  mypal <- colorRampPalette(ggsci::pal_nejm("default", alpha = 0.8)(8))(18)
  #scales::show_col(mypal)
  pathway_col <- data.frame(Pathway=c("Cell_Cycle", "Chromatin_modifiers", "Notch", "NRF2", "PI3K", "RTK_RAS", "TGF-Beta", "TP53", "WNT", "Hippo", "Myc",NA), col=c(mypal[c(1:8,10,15:16)],"grey70"))
  pathway_col_values <- pathway_col$col
  names(pathway_col_values) = pathway_col$Pathway
  
  gene_edges = edges[pValue<=p_value]

  from_id <- gene_edges[,sapply(gene1,function(x) as.numeric(nodes[label == x][,.(id)]) )]
  to_id <- gene_edges[,sapply(gene2,function(x) as.numeric(nodes[label == x][,.(id)]) )]
  
  ### get layout
  edges <- data.table(from = from_id, 
                      to = to_id,
                      color = ifelse(gene_edges[,.(Event)] == "Co_Occurence", "#4058a3","#f18826")) 
  edges <- cbind(edges,gene_edges)
  
  ### get nodes in Co_Occurence/Mutually_Exclusive connection
  co_nodes <- unique(c(t(edges[which(edges[,.(Event)] == "Co_Occurence"), .(gene1,gene2)]))) 
  mu_nodes <- unique(c(t(edges[which(edges[,.(Event)] == "Mutually_Exclusive"), .(gene1,gene2)]))) 

  if (length(co_nodes)>0) {
    nodes_co <- nodes[label %in% co_nodes]
    edges_co <- edges[Event == "Co_Occurence" & gene1 %in% co_nodes & gene2 %in% co_nodes]
    nodes_mu <- nodes[label %in% mu_nodes | freq==1] 
    edges_mu <- edges[Event == "Mutually_Exclusive" & gene1 %in% mu_nodes & gene2 %in% mu_nodes]
    net_co <- graph_from_data_frame(d = edges_co, vertices = nodes_co, directed = F) 
  # Merge group of co-occurrence nodes label to whole nodes
  } else {
    nodes_mu <- nodes[label %in% mu_nodes | freq==1] 
    edges_mu <- edges[Event == "Mutually_Exclusive" & gene1 %in% mu_nodes & gene2 %in% mu_nodes]
    net_co <- graph_from_data_frame(d=edges_mu, vertices=nodes_mu, directed=F)
  }
  # Calculate connected groups for nodes with co-occurrences relationship
  tg_co <- activate(tidygraph::as_tbl_graph(net_co),nodes) %>%
    mutate(label=label,name = as.numeric(name), community_id = group_components())
  
  ## merge connected group information in all nodes, factor 0 means no co-occurrences connection
  gene_nodes = left_join(nodes,data.frame(label=V(tg_co)$label,community_id=V(tg_co)$community_id),by="label") %>%
    mutate(community_id=as.factor(ifelse(is.na(community_id),0,community_id)))
  colnames(gene_nodes)[3:4] = c("group","color")
  
  edges = edges[gene1 %in% nodes$label & gene2 %in% nodes$label]
  net = graph_from_data_frame(d = edges, vertices = gene_nodes, directed = F)
  tg <- tidygraph::as_tbl_graph(net) %>% activate(nodes) %>% mutate(label=label)
  
  if (no_layout) l = net_to_layout(net,gene_nodes,simple=TRUE) else
    l = net_to_layout(net,gene_nodes)

  v.size <- V(tg)$freq_value %>% as.character() %>% as.numeric()
  e.color <- E(tg)$color
  E(tg)$weight <- E(tg)$log10pval %>% rescale(c(0.5,2.5))
 
  #p <- ggraph(tg, layout = l)
  p <- ggraph(tg,'manual', x= V(tg)$ccf_type,y=round(V(tg)$freq,2))
  
  if (nrow(nodes_co)>0) { #"#4058a3"
    p = p +  
      geom_edge_fan0(aes(filter = (Event=="Co_Occurence" & pValue<0.01)),edge_colour = "#006D77",edge_width=1) +
      geom_edge_fan0(aes(filter = (Event=="Co_Occurence" & pValue>=0.01 & pValue<0.05)),edge_colour = "#006D77",edge_width=0.5)
  }
  if (nrow(nodes_mu)>0) { #f18826   
    p = p +
      #geom_edge_fan0(aes(filter = (Event=="Mutually_Exclusive" & pValue<0.01)),edge_colour = "#E29578",edge_width=1) +
      #geom_edge_fan0(aes(filter = (Event=="Mutually_Exclusive" & pValue>=0.01 & pValue<0.05)),edge_colour = "#E29578",edge_width=0.5)
      geom_edge_diagonal0(aes(filter = (Event=="Mutually_Exclusive" & pValue<0.01)),edge_colour = "#E29578",edge_width=1) +
      geom_edge_diagonal0(aes(filter = (Event=="Mutually_Exclusive" & pValue>=0.01 & pValue<0.05)),edge_colour = "#E29578",edge_width=0.5)
  }
  #subtype_var_col=c("#aa93af","#9ecac8","#fef29f")
  # if (is.na(subtype_var_col[1])) {
  #   subtype_var_col=c("mediumorchid4","darkgreen","#fef29f")
  #   names(subtype_var_col) = subtype_var_content
  # }
  subtype_var_content = names(subtype_var_col)

  max_ylim = ifelse(max(V(tg)$freq)>0.9,max(V(tg)$freq),0.9)
  p <- p +
    ggtitle(title)+
    scale_color_manual(name = "Pathway",values = pathway_col_values)+
    scale_y_sqrt(breaks=c(0.05,0.1,0.2,0.4,0.6,0.9),labels = percent,limits=c(0.04, max_ylim))+
    #geom_node_point(aes(size=as.numeric(as.character(freq_value))),alpha = 0.9) +
    #guides(shape = guide_legend(override.aes = list(size = 15:23)))+
    scatterpie::geom_scatterpie(
      aes(x=ccf_type,y=freq,r=freq_value),
      cols = c(subtype_var_content),
      data = as_data_frame(tg, "vertices") %>% mutate(freq_value=rescale(as.numeric(freq_value),c(0.01,0.03))),
      colour = NA,pie_scale = 0.5
    ) +
    geom_node_text(aes(label = label,colour=group), repel = TRUE,
                   size = log(v.size)+0.5,vjust=-2,hjust=-2,fontface="bold")+
    theme_pubclean()+
    scale_x_reverse(limits=c(NA,0.5)) +
    labs(x=ccf_type,y = "Sample Frequency", subtitle = subtitle) +
    scale_fill_manual(name = subtype_var, values = subtype_var_col)+
    theme(legend.position = "bottom")
 
  # plot code here with geom_scatterpie
  
  ## Plot group information
  if (fill) {p <- p +
    ggforce::geom_mark_hull(aes(filter=community_id!=0,x=x,y=y,fill = as.factor(community_id),label = "Co-occurrence"),
                            colour = NA,con.colour = "grey",show.legend = FALSE,concavity = 4,
                            expand = unit(3, "mm"),alpha = 0.1)
  }
  
  if (facet) {
    p = p+facet_nodes(~community_id,scales="free_y",ncol=1)+
      theme(strip.text.x = element_blank(),strip.background = element_blank())
  }
  p
}

mafToNodeEdgesGraph <- function(maf,gistic=NULL,pathway=FALSE,top_gene=30,oncoplot.feature=NA,
                           fabcolors=NA,oncoPath,bottom_barplot=TRUE,bottom_feature,p_value=0.01,
                           subtype_var=NA,ggraph=TRUE) {
    
    # 0.Preprocessing
    ## 0.1 Colors
    pathway_col <- data.frame(Pathway=c("Cell_Cycle", "Chromatin_modifiers", "Notch", "NRF2", "PI3K", "RTK_RAS", "TGF-Beta", "TP53", "WNT", "Hippo", "Myc",NA), col=c(mypal[c(1:8,10,15:16)],"grey70"))
    pathway_col_values <- pathway_col$col
    names(pathway_col_values) = pathway_col$Pathway
    
    ## 0.2 Data format
    if (is.data.frame(oncoPath)) oncoPath = as.data.table(oncoPath)
    mafTemp = maf@data 
    n_sample = nrow(maf@clinical.data)
    colnames(mafTemp)[which(colnames(mafTemp)==subtype_var)] = "pie_subtype"
    ## 0.3 If gistic exist load gistic file 
    if (!is.null(gistic)) gistic = subset(gistic,variable %in% getSampleSummary(x = maf)$Tumor_Sample_Barcode )
    
    # 1.Co-occurrence Analysis
    ##1.1 Run Co-occurrence Analysis
    if (pathway) {
      colnames(mafTemp)[c(1,which(colnames( mafTemp)=="Pathway"))] <- c("Hugo_Symbol_Gene","Hugo_Symbol")
      maf = read.maf(mafTemp,clinicalData = maf@clinical.data)
      res <- getInteractions(maf, top=10, colPal = "PRGn", returnAll = TRUE, sigSymbolsSize=1.5)
      bottom.feature=NULL
    } else {
      maf = read.maf(mafTemp,clinicalData = maf@clinical.data)
      res <- getInteractions(maf, top=top_gene, colPal = "PRGn", returnAll = TRUE, sigSymbolsSize=1.5) 
    }
    
    ##1.2 Draw interaction heatmap
    p_interheatmap <- interaction_heatmap(res)
  
    ##1.3 Calculate gene sample frequency
    mafTemp1 = mafTemp[!is.na(Hugo_Symbol)][order(Tumor_Sample_Barcode,Hugo_Symbol,ccube_ccf, decreasing = TRUE)]
    ### Set keys - this sorts the data based on these values
    setkeyv(mafTemp1, c('Tumor_Sample_Barcode','Hugo_Symbol'))
    ### keep unique observations
    geneccfFreq <- subset(unique(mafTemp1))[,.(.N,mean_ccf=mean(ccube_ccf,na.rm=TRUE),median_ccf=median(ccube_ccf,na.rm=TRUE)), by = .(Hugo_Symbol)]
    
    if (!is.na(subtype_var)) {
      subtypeFreq = mafTemp1[, .(.N),by=.(Hugo_Symbol,pie_subtype)]
      subtypeFreq = data.table::dcast(data = subtypeFreq, formula = Hugo_Symbol~pie_subtype, fill = 0, value.var = 'N')
      subtypeFreq[,(subtype_var_content[which(subtype_var_content %nin% colnames(subtypeFreq))]):= 0 ]
      geneccfFreq <- geneccfFreq[subtypeFreq,on="Hugo_Symbol"]
    }
    
    geneccfFreq[,freq:= N/n_sample]
    geneccfFreq = geneccfFreq[freq>=0.05]
    geneccfFreq[,freq_value:=cut(freq,breaks=c(seq(0,1,length.out=10)),labels=seq(15,50,length.out=9))]
    colnames(geneccfFreq)[1] = "label"
    #top_genes = dplyr::pull(geneccfFreq[order(-freq)][1:5, .(label)])
    top_genes = dplyr::pull(geneccfFreq[order(-freq)][,.(label)])
    
 
    ##1.4 Draw oncoplot
    if (is.null(bottom_feature)) bottom_barplot = FALSE else bottom_barplot = TRUE
    p_oncoplot = as.ggplot(function()
        if (all(!is.na(oncoplot.feature))) {
          oncoplot(maf = maf, top=top_gene, removeNonMutated = FALSE, #colors =  mycolors, 
                   drawRowBar = TRUE, annotationColor= fabcolors,
                   keepGeneOrder = FALSE, bgCol="#F5F5F5",sortByAnnotation = TRUE,
                   clinicalFeatures=oncoplot.feature,bottom_barplot = bottom_barplot,bottom.feature  = bottom_feature,
                   annotationFontSize = 1.1,legendFontSize = 1.1)
                   
        } else {
          oncoplot(maf = maf, top=top_gene, removeNonMutated = FALSE, #colors = mycolors, 
                   annotationColor=fabcolors,  
                   drawRowBar = TRUE,sortByAnnotation = TRUE, keepGeneOrder = FALSE, bgCol="#F5F5F5",bottom_barplot = bottom_barplot,
                   bottom.feature  = bottom_feature)
        }
    )

    ##1.5 Get nodes and edges data
    gene_edges <- res$pairs[pValue < p_value][gene1 %in% geneccfFreq$label & gene2 %in% geneccfFreq$label ]
    #gene_nodes <- unique(as.character(c(t(cbind(gene_edges$gene1, gene_edges$gene2)),top_genes))) 
    gene_nodes <- unique(as.character(t(cbind(gene_edges$gene1, gene_edges$gene2)))) 
    gene_nodes_col <- data.table(Pathway = gene_nodes) 


    if (!pathway) {
      # Gene format
      gene_nodes_col = as.data.table(pathway_col)[oncoPath[data.table(Gene = gene_nodes),on="Gene"],on="Pathway"][!duplicated(Gene)]  
    } else {
      # Pathway format
      gene_nodes_col = data.table(Pathway = gene_nodes)[pathway_col, on="Pathway"]
    }

    if (length(gene_nodes) > 0) {
      nodes <- data.table(id = 1:length(gene_nodes),
                          label = gene_nodes,
                          group = gene_nodes_col[, .(Pathway)],
                          color = gene_nodes_col[, .(col)],
                          shadow = TRUE) 
      nodes <- nodes[geneccfFreq, on = "label",nomatch=0]
      nodes[,label_show := paste0(label,"\n ",round(freq,2)*100,"%")]
      colnames(nodes)[which(colnames(nodes)==ccf_type)] = "ccf_type"
      nodes <- nodes[ccf_type>0]       
      nodes[,value := cut(ccf_type,breaks = (1:10) / 10,labels = 1:9)]
      nodes[,ccf_type := round(ccf_type, 3)] 
      edges = as.data.table(gene_edges)[gene1 %in% dplyr::pull(nodes[,.(label)]) & gene2 %in% dplyr::pull(nodes[,.(label)])][] 
      
      if (ggraph) {
        title = paste0(unique(maf@clinical.data$cancertype),"|",unique(maf@clinical.data$subtype)," - ",nrow(maf@clinical.data)," samples")
        subtype_var_col=  oncoplot.feature.col[[subtype_var]]
       
        p_network = ggraph_plot(nodes=nodes,edges=edges,fill=TRUE,p_value=p_value,
                                subtype_var_col=subtype_var_col,title=title) 
        grid.newpage()
        grid.draw(ggplotGrob(p_network))
        popfunction(p_interheatmap,position="x = 0.4, y=1,width = 0.55,height = 0.5")
        p_network = grid.grab(wrap.grobs = TRUE) %>% as.ggplot() 
      }
      
      return(list(nodes, gene_edges, p_oncoplot, p_interheatmap,p_network,geneccfFreq))
    } else {
      nodes <- data.frame(id = NA)
      return(list(nodes, nodes, p_oncoplot, p_interheatmap,NA,geneccfFreq))
    }
    
    
 
}

interaction_heatmap <- function(res) {

  melted_cormat <- reshape2::melt(res$interactions, na.rm = TRUE) %>%
    mutate(text_0.01=ifelse(abs(value)>-log10(0.01),"*",NA),
           text_0.05=ifelse(abs(value)>-log10(0.05),"·",NA))

  melted_cormat %>%
    ggplot()+geom_tile(aes(Var1, Var2, fill = value),color="White")+
    geom_text(aes(x=Var1,y=Var2,label=text_0.01),color="white",size=6,vjust=0.75)+
    geom_text(aes(x=Var1,y=Var2,label=text_0.05),color="white",size=6)+
    scale_fill_gradient2(low = "#006D77", high = "#E29578", mid = "white", 
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

 
Coocrrence_network_plot_MSI <- function(types,manSelect,gistic=NULL,ks_path,ccf_path,oncoPath,output=NA,ccf_type,p_value=0.05,top_gene,pathway,oncoplot.feature,oncoplot.feature.col, bottom_feature) {
  
  print(types)
  if (nchar(types)==8) typesMaf = c(substr(types,1,4),substr(types,5,8)) else typesMaf = types
  
  # get significant transited gene list
  ksByCancer = data.table::fread(ks_path)[cancertype %in% paste0(sort(typesMaf),collapse  = "") & (significant_ks == 1|significant_fisher==1)]
  geneCounts <- ksByCancer %>% count(Gene, sort = TRUE, name = "n_genes_across_cancertype") %>% rbind(list("POLE",1))
  
  #load ccf files
  for (j in 1:length(ccf_path) ) {
    ccf_load = get(load(ccf_path[j]))
    if (j==1) ccf_all= ccf_load else ccf_all = rbind(ccf_all,ccf_load)
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
  mergeMaf <- merge_mafs(lstObject)
  mergeMaf@data = mergeMaf@data[,Tumor_Sample_Barcode:=substr(Tumor_Sample_Barcode,1,12)]
  mergeMaf = subsetMaf(mergeMaf,tsb=manSelect$samplename)
 
  mafTemp <- as.data.frame(mergeMaf@data) %>%
    mutate(Tumor_Sample_Barcode = substr(as.character(Tumor_Sample_Barcode),1,12),
           Variant_Classification=as.character(Variant_Classification)) %>%
    left_join(ccf_filter[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode","ccube_ccf")],by=c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")) %>%
    mutate(ccube_ccf=ifelse(ccube_ccf>=1,1,ccube_ccf),
           Variant_Classification=as.character(Variant_Classification)) %>%
    filter(Variant_Classification %nin% c('Silent','Translation_Start_Site',"3'Flank",
                                          "3'UTR","5'Flank","5'UTR","Intron","Nonstop_Mutation","RNA","IGR","Splice_Region"),
           !is.na(Tumor_Sample_Barcode),
           Hugo_Symbol %in% geneCounts$Gene) %>%
    arrange(desc(Tumor_Sample_Barcode,Hugo_Symbol,ccube_ccf)) %>%
    distinct(Tumor_Sample_Barcode,Hugo_Symbol,.keep_all = TRUE) %>%  # 526 samples
    left_join(oncoPath, by=c("Hugo_Symbol"="Gene")) %>%
    left_join(manSelect[,c("Tumor_Sample_Barcode",subtype_var)], by = c("Tumor_Sample_Barcode")) 
  
  clinicalTemp <- getSampleSummary(x = mergeMaf)[,Tumor_Sample_Barcode := substr(as.character(Tumor_Sample_Barcode),1,12)] 
  clinicalTemp <- clinicalTemp[manSelect,on="Tumor_Sample_Barcode"]
  
  mergeMaf= read.maf(maf = mafTemp,clinicalData = clinicalTemp)

  ## No filter
  mafEsig4 <- subsetMaf(mergeMaf,tsb=subset(manSelect,subtype=="ESig4")$samplename) 
  mafEsig3 <- subsetMaf(mergeMaf,tsb=subset(manSelect,subtype=="ESig3")$samplename) 
  
  # Perform interaction analysis using oncogenic pathway data

  NodeEdges4 = mafToNodeEdgesGraph(maf=mafEsig4,top_gene=30,gistic=gistic,
                                   oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,
                                   bottom_feature=bottom_feature,subtype_var="Subtype")
  
  NodeEdges3 = mafToNodeEdgesGraph(maf=mafEsig3,top_gene=30,gistic=gistic,
                                   oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,
                                   oncoPath=oncoPath,bottom_feature=bottom_feature,subtype_var="Subtype")
  
  NodeEdges4_pathway = mafToNodeEdgesGraph(maf=mafEsig4,pathway=TRUE,top_gene=10,gistic=gistic,
                                           oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,
                                           oncoPath=oncoPath,bottom_feature=bottom_feature,subtype_var="Subtype")
  NodeEdges3_pathway = mafToNodeEdgesGraph(maf=mafEsig3,pathway=TRUE,top_gene=10,gistic=gistic,
                                           oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,
                                           oncoPath=oncoPath,bottom_feature=bottom_feature,subtype_var="Subtype")
  pathway3 = pathway4 = FALSE 
  
  p_network = ggraph_plot(nodes=nodes,edges=edges,fill=TRUE,p_value=p_value,
                          subtype_var_content=subtype_var_content,subtype_var_col=subtype_var_col) 
  
  ## Plot ccf transition with pie chart
  edges = rbind(NodeEdges4[[2]] %>% mutate(ESig="ESig4"),NodeEdges3[[2]] %>% mutate(ESig="ESig3"))
  nodes = rbind(NodeEdges4[[1]] %>% mutate(ESig="ESig4"),NodeEdges3[[1]] %>% mutate(ESig="ESig3"))
  
  geneccfFreq = rbind(NodeEdges4[[6]] %>% mutate(ESig="ESig4"),NodeEdges3[[6]] %>% mutate(ESig="ESig3")) %>%
    filter(label %in% NodeEdges4[[6]]$label & label %in% NodeEdges3[[6]]$label) 
  genelevels =  geneccfFreq %>%
    group_by(label) %>% dplyr::summarise(avg_freq = mean(freq)) %>% 
    arrange(desc(avg_freq)) 
  genelevels$id =  (1:nrow(genelevels))/40
  geneccfFreq = left_join(geneccfFreq,genelevels,by="label") 
  arrowdf = geneccfFreq %>%
    select(label,mean_ccf,avg_freq,ESig,id) %>%
    reshape2::dcast(label+id+avg_freq~ESig,value.var="mean_ccf") %>%
    mutate(direction = ESig4-ESig3) %>%
    mutate(direction=ifelse(direction>0,direction-0.01,direction+0.01),
           start=ifelse(direction>0,ESig3+0.01,ESig3-0.01))
  ## setting arrow
  subtype_var_col=  oncoplot.feature.col[[subtype_var]]
  
  geneccfFreq %>%
    ggplot() +
    geom_segment(data=arrowdf,aes(x=start,xend=ESig3+direction,y=id,yend=id),arrow=arrow(length=unit(0.30,"cm"),type="closed"),size=1)+
    #geom_text(aes(label=labelg))+
    scatterpie::geom_scatterpie(
      aes(x= mean_ccf,y=id,group=id,r=0.01),
      cols = c(subtype_var_content),
      data = as.data.frame(geneccfFreq),
      colour = NA,
      pie_scale = 0.1
    ) +
    theme_pubr()+
    scale_fill_manual(name = subtype_var, values = subtype_var_col)+
    scale_y_continuous(breaks=genelevels$id,labels=genelevels$label)+
    labs(y="",x="mean of CCF")
  
  NodeEdges4[[5]]/NodeEdges3[[5]]
  # gene network
  if (nrow(NodeEdges4[[2]])>0 & ncol(NodeEdges4[[2]])>1) {
    NodeEdges4[[4]]
    p1 = ggraph_plot(nodes=NodeEdges4[[1]],edges=NodeEdges4[[2]],fill=TRUE,p_value=p_value,subtype_var="Subtype",title="ESig4") 
    grid.draw(ggplotGrob(p1))
    popfunction(NodeEdges4[[4]],position="x = 0.4, y=1,width = 0.6,height = 0.55")
    p1 = grid.grab(wrap.grobs = TRUE) %>% as.ggplot() 
  }
  
  
  if (nrow(NodeEdges3[[2]])>0 & ncol(NodeEdges3[[2]])>1) {
    p2 = ggraph_plot(nodes=NodeEdges3[[1]],edges=NodeEdges3[[2]],fill=TRUE,p_value=p_value,title="ESig3") 
    grid.draw(ggplotGrob(p2))
    popfunction(NodeEdges3[[4]],position="x = 0.4, y=1,width = 0.6,height = 0.55")
    p2 = grid.grab(wrap.grobs = TRUE) %>% as.ggplot() 
  
  }
  
  # Pathway network
  if (nrow(NodeEdges4_pathway[[2]])>0 & ncol(NodeEdges4_pathway[[2]])>1) {
    p5=ggraph_plot(nodes=NodeEdges4_pathway[[1]],edges=NodeEdges4_pathway[[2]],fill=TRUE,p_value=0.01,title="Pathway") 
    pathway4 = TRUE 
  } 
  if (nrow(NodeEdges3_pathway[[2]])>0 & ncol(NodeEdges3_pathway[[2]])>1) {
    p6=ggraph_plot(nodes=NodeEdges3_pathway[[1]],edges=NodeEdges3_pathway[[2]],fill=TRUE,p_value=0.01,title="Pathway")
    pathway3 = TRUE
  } 
  
  if (!is.null(gistic)) gistic = subset(gistic,Sample %in% geneCounts$Gene)
  #gistic = NULL
  
  # oncoplot
  p31=NodeEdges4[[3]];p41=NodeEdges3[[3]]
  p32=NodeEdges4_pathway[[3]];p42=NodeEdges3_pathway[[3]]
  
  
  ################################ Network Analysis for subtype
  if (types %in% c("COADREAD","STAD","UCEC"))  {
  
      ESig4MSS_samples = subset(manSelect,subtype=="ESig4" & Subtype=="MSS")$samplename
      ESig4MMR_samples = subset(manSelect,subtype=="ESig4" & Subtype=="MMR")$samplename
      ESig4POLE_samples = subset(manSelect,subtype=="ESig4" & Subtype=="POLE")$samplename
      ESig3MSS_samples = subset(manSelect,subtype=="ESig3" & Subtype=="MSS")$samplename
      ESig3MMR_samples = subset(manSelect,subtype=="ESig3" & Subtype=="MMR")$samplename
      ESig3POLE_samples = subset(manSelect,subtype=="ESig3" & Subtype=="POLE")$samplename
      MSS4 = MMR4 = POLE4 = MSS3 = MMR3 = POLE3 = FALSE
      
      if (length(ESig4MSS_samples)>0) {                      
        mafEsig4MSS <- subsetMaf(mafEsig4, tsb=ESig4MSS_samples)
        # get nodes and edges
        NodeEdges4MSS =  mafToNodeEdges(maf=mafEsig4MSS,top_gene=30,gistic=gistic,
                                      oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
        if (nrow(NodeEdges4MSS[[2]])>0 & ncol(NodeEdges4MSS[[2]])>1) {
          p11=ggraph_plot(nodes=NodeEdges4MSS[[1]],edges=NodeEdges4MSS[[2]],fill=TRUE,title="",no_layout = TRUE,subtitle="MSS") 
          MSS4 = TRUE
        }
      } 
    
      if (length(ESig4MMR_samples)>0) {
        mafEsig4MMR <- subsetMaf(mafEsig4, tsb=ESig4MMR_samples) 
        NodeEdges4MMR =mafToNodeEdges(maf=mafEsig4MMR,top_gene=30,gistic=gistic,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
        if (nrow(NodeEdges4MMR[[2]])>0 & ncol(NodeEdges4MMR[[2]])>1 & nrow(NodeEdges4MMR[[1]])>0) {
          p12=ggraph_plot(nodes=NodeEdges4MMR[[1]],edges=NodeEdges4MMR[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MMR")
          MMR4 = TRUE
        }
      }
    
      if (length(ESig4POLE_samples)>0) {
        mafEsig4POLE <- subsetMaf(mafEsig4, tsb=ESig4POLE_samples)
        NodeEdges4POLE=mafToNodeEdges(mafEsig4POLE,top_gene=30,gistic=gistic,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
        if (nrow(NodeEdges4POLE[[2]])>0 & ncol(NodeEdges4POLE[[2]])>1) {
          p13=ggraph_plot(nodes=NodeEdges4POLE[[1]],edges=NodeEdges4POLE[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="POLE") 
          POLE4 = TRUE
        }
      } 
      
      if (length(ESig3MSS_samples)>0) {                      
        mafEsig3MSS <- subsetMaf(mafEsig3, tsb=ESig3MSS_samples)
        NodeEdges3MSS =mafToNodeEdges(mafEsig3MSS,top_gene=30,gistic=gistic,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
        if (nrow(NodeEdges3MSS[[2]])>0 & ncol(NodeEdges3MSS[[2]])>1) {
          p21=ggraph_plot(nodes=NodeEdges3MSS[[1]],edges=NodeEdges3MSS[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MSS") 
          MSS3 = TRUE
        }
      } 
      
      if (length(ESig3MMR_samples)>0) {
        mafEsig3MMR <- subsetMaf(mafEsig3, tsb=ESig3MMR_samples) 
        NodeEdges3MMR =mafToNodeEdges(maf=mafEsig3MMR,top_gene=30,gistic=gistic,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
        if (nrow(NodeEdges3MMR[[2]])>0 & ncol(NodeEdges3MMR[[2]])>1) {
          p22=ggraph_plot(nodes=NodeEdges3MMR[[1]],edges=NodeEdges3MMR[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="MMR")
          MMR3 = TRUE
        }
      }
      
      if (length(ESig3POLE_samples)>0) {
        mafEsig3POLE <- subsetMaf(mafEsig3, tsb=ESig3POLE_samples)
        NodeEdges3POLE = mafToNodeEdges(maf=mafEsig3POLE,top_gene=30,gistic=gistic,oncoplot.feature=oncoplot.feature,fabcolors = oncoplot.feature.col,oncoPath=oncoPath,bottom_feature=bottom_feature)
        if (nrow(NodeEdges3POLE[[2]])>0 & ncol(NodeEdges3POLE[[2]])>1) {
          p23 = ggraph_plot(nodes = NodeEdges3POLE[[1]],edges=NodeEdges3POLE[[2]],fill=TRUE,p_value=0.05,title="",no_layout = TRUE,subtitle="POLE") 
          POLE3 = TRUE     
        }
        
      } else {
      
      ## Final Plot
      grid.newpage()
    
      popfunction(p1,position="x = 0, y=1,width = 0.25,height = 0.25")       # ESig4 - Gene Network
      
      if (MSS4) popfunction(p11,position="x = 0, y=0.75,width = 0.16,height = 0.2")    # ESig4 - Subtype1 Gene Network
      if (MMR4) popfunction(p12,position="x = .17, y=0.75,width = 0.16,height = 0.2")   # ESig4 - Subtype2 Gene Network
      if (POLE4) popfunction(p13,position="x = .34, y=0.75,width = 0.16,height = 0.2")   # ESig4 - Subtype3 Gene Network
      if (pathway4) popfunction(p5,position="x = 0.25, y=1,width = 0.3,height = 0.25")     # ESig4 - Pathway Network
      popfunction(p31,position="x = 0.5, y=1,width = 0.5,height = 0.25")    # ESig4 - Gene Oncoplot
      popfunction(p32,position="x = 0.5, y=0.75,width = 0.5,height = 0.25") # ESig4 - Pathway Oncoplot
      popfunction(p2,position="x = 0, y=0.5,width = 0.25,height = 0.25")     # ESig3 - Gene Network
      if (MSS3) popfunction(p21,position="x = 0, y=0.25,width = 0.16,height = 0.2")    # ESig3 - Subtype1 Gene Network
      if (MMR3) popfunction(p22,position="x = .17, y=0.25,width = 0.16,height = 0.2")   # ESig3 - Subtype2 Gene Network
      if (POLE3) popfunction(p23,position="x = .34, y=0.25,width = 0.16,height = 0.2")   # ESig3 - Subtype3 Gene Network
      if (pathway3) popfunction(p6,position="x = 0.25, y=0.5,width = 0.3,height = 0.25")   # ESig3 - Pathway Network
      popfunction(p41,position="x = 0.5, y=0.5,width = 0.5,height = 0.25")  # ESig3 - Gene Oncoplot
      popfunction(p42,position="x = 0.5, y=0.25,width = 0.5,height = 0.25") # ESig3 - Pathway Oncoplot
      
      caption=paste0("Top ",top_gene," genes with significant ccf transition. \n Edge width represents the strength of p value, only edges with p<",p_value," were shown. \n Circle size represent population frequency of each gene.")
      #subtitle = 
      
      g2 <- grid.grab(wrap.grobs = TRUE) %>% as.ggplot() +
        labs(caption=caption)
      
      }  
  }
 
      ############################### Final Plot Template
      grid.newpage()
      popfunction(p1,position="x = 0, y=1,width = 0.25,height = 0.25")       # ESig4 - Gene Network
      if (pathway4) popfunction(p5,position="x = 0.25, y=1,width = 0.3,height = 0.25")     # ESig4 - Pathway Network
      popfunction(p31,position="x = 0.5, y=1,width = 0.5,height = 0.25")    # ESig4 - Gene Oncoplot
      popfunction(p32,position="x = 0.5, y=0.75,width = 0.5,height = 0.25") # ESig4 - Pathway Oncoplot
      popfunction(p2,position="x = 0, y=0.5,width = 0.25,height = 0.25")     # ESig3 - Gene Network
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

    

