

Coocrrence_network_plot <- function(types,output=NA,ksByCancer,manSelect,ccf_type,p_value,top_gene,ccf_all,oncoPath) {
  
  ggraph_plot = function(tg,layer,facet=FALSE,title=types,fill=FALSE,zoom=NA) {
    v.size <- V(tg)$freq_value %>% as.character() %>% as.numeric()
    e.color <- E(tg)$color
    E(tg)$weight <- E(tg)$log10pval %>% rescale(c(1,4))
    #eigenCent <- evcent(tg)$vector
    #bins <- unique(quantile(eigenCent, seq(0,1,length.out=30)))
    #vals <- cut(eigenCent, bins, labels=FALSE, include.lowest=TRUE)
    #colorVals <- rev(heat.colors(length(bins)))[vals]
    
    
    p = tg %>%
      ggraph(layout=layer) +
      geom_edge_link0(aes(filter=(Event=="Co_Occurence")),edge_colour="#4058a3",width=E(tg)$weight[which(E(tg)$Event=="Co_Occurence")]) +
      geom_edge_link0(aes(filter=(Event=="Mutually_Exclusive")),edge_colour="#f18826",width=E(tg)$weight[which(E(tg)$Event=="Mutually_Exclusive")])+
      ggtitle(title)+
      geom_node_point(colour=V(tg)$color,size=v.size/2,alpha = 0.9) +
      geom_node_text(aes(label = label), #repel = TRUE,
                     size=log(v.size),colour='black')+
      theme_classic()+
      scale_x_reverse()+
      labs(x=ccf_type,y="",caption=paste0("Edge width represent the strength of p value, only edges with p>",p_value," were shown. \n Circle size represent population frequency of each gene."),
           subtitle = paste0("Top ",top_gene," genes with significant ccf transition"))+
      theme(axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y =element_blank())
    
    if (fill) {p = p+   
      ggforce::geom_mark_hull(aes(filter=community_id!=0,x=x,y=y,fill = as.factor(community_id),label = paste0("trajectory ", community_id)),
                              colour = NA,con.colour = "grey",show.legend = FALSE,concavity = 4,
                              expand = unit(3, "mm"),alpha = 0.25)
    }
    
    if (facet) {
      p = p+facet_nodes(~community_id,scales="free_y",ncol=1)+
        theme(strip.text.x = element_blank(),strip.background = element_blank())
    }
    p
  }
  
  mypal <- colorRampPalette(pal_nejm("default", alpha = 1)(8))(18)
  pathway_col = data.frame(Pathway=c("Cell_Cycle", "Chromatin_modifiers", "Notch", "NRF2", "PI3K", "RTK_RAS", "TGF-Beta", "TP53", "WNT", "Hippo", "Myc",NA), col=c(mypal[1:11],"grey70"))
  
  print(types)
# 
  if (nchar(types)==8) typesMaf = c(substr(types,1,4),substr(types,5,8)) else typesMaf = types
  
  # Load MAF
  if (types != "All") {
    lstObject <- lapply(typesMaf, function(i) tcga_load(i))
  } else {
    lstObject <- lapply(tcga_available()$Study_Abbreviation[1:33], function(i) tcga_load(i))
  }
  
  ### Load selected gene ccf file
  geneCounts <- ksByCancer %>% count(Gene, sort = TRUE, name = "n_genes_across_cancertype")
 
  geneccfFreq <- ccf_all %>%
    mutate(ccube_ccf= ifelse(ccube_ccf>=1,1,ccube_ccf))%>%
    filter(Variant_Classification %nin%  c('Translation_Start_Site'),
           Hugo_Symbol %in% ksByCancer$Gene) %>%
    group_by(Hugo_Symbol) %>%
    dplyr::summarise(mean_ccf=mean(ccube_ccf,na.rm=TRUE),median_ccf=median(ccube_ccf,na.rm=TRUE),n=length(unique(samplename))) %>% 
    mutate(freq=n/nrow(manSelect)) %>%
    mutate(freq_value=cut(.$freq,breaks=c(seq(0,1,length.out=10)),labels=15:23))
  colnames(geneccfFreq)[1] = "label"
  
  mergeMaf<- merge_mafs(lstObject)
  mergeMaf@data$Tumor_Sample_Barcode = substr(as.character(mergeMaf@data$Tumor_Sample_Barcode),1,12)
  
  ## No filter
  mafEsig4 <- subsetMaf(mergeMaf, tsb=manSelect$samplename, genes = geneCounts$Gene,
                        query = "Variant_Classification %nin%  c('Translation_Start_Site')")
  
  # Perform somatic interaction analysis - limit to top 30 mutated genes (this can be expanded)
  res <- getInteractions(mafEsig4, top=top_gene, colPal = "PRGn", returnAll = TRUE, sigSymbolsSize=1.5) 
  
  colnames(oncoPath)[1] = "gene1"
  gene_edges = res$pairs %>% filter(pValue<p_value) %>%
    left_join(oncoPath, by="gene1") %>%
    left_join( pathway_col,by="Pathway") 
  
  gene_nodes = cbind(gene_edges$gene1,gene_edges$gene2) %>% t(.) %>% as.character() %>%unique(.)
  # Concurrence Network visulization
  gene_nodes_col = data.frame(gene1=gene_nodes) %>%  
    left_join(oncoPath, by="gene1") %>%
    left_join( pathway_col,by="Pathway") %>%
    distinct(gene1,.keep_all = TRUE)
  
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
           
    gene_edges = gene_edges %>% filter(gene1 %in% nodes$label,gene2 %in% nodes$label) 
    from_id = sapply( gene_edges$gene1,function(x) nodes %>% filter(label==x) %>% .$id ) %>% as.numeric()
    to_id = sapply( gene_edges$gene2,function(x) nodes %>% filter(label==x) %>% .$id ) %>% as.numeric()

   
    
    # Multi-layer network
    get_layers = function(net) {
      layers_name = unique(V(net)$value) %>% as.numeric() %>% sort()
      layers_id = data.frame(name=layers_name,layers=length(layers_name):1)
      layer = data.frame(name = as.numeric(V(net)$value)) %>% left_join(layers_id,by="name") %>% .$layers
      return(layer)
    }
    
    
    #layout = layout_with_sugiyama(net, layers=get_layers(net))
    check_layout_facet <- function(layout,node,facet=NA,layer1_threshold=0.99) {
      x_level = max(layout[,1])
      y_level = length(unique((layout[,facet])))
      df = layout;colnames(df) = c("x_level","y","id","y_level","freq",ccf_type)
      
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
     
      df %>% mutate(x_level=as.numeric(node$median_ccf )) 
    }
    layout_to_matrix <- function(layout) {layout %>% arrange(id) %>% .[,1:2] %>% as.matrix()}

    net_to_layout <- function(net,node){
      cbind(get_layers(net),layout_with_sugiyama(net, layers=V(net)$value)$layout[,1]) %>%
        as.data.frame() %>% cbind(node[,c("id","community_id","freq",ccf_type)]) %>%
        mutate(community_id=as.numeric(as.character(community_id))) %>% 
        check_layout_facet(facet="community_id",node=node)  %>%
        layout_to_matrix()
    }
    ### get layout
    # net = graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
    # 
    # V(net)$color <- nodes$color
    # #V(net)$size <- nodes$median_c %>% as.numeric() %>% rescale(.,c(10,20))
    # V(net)$size <- nodes$freq_value %>% as.character() %>% as.numeric()
    # E(net)$width <- edges$oddsRatio %>% rescale(c(2,5))
    # E(net)$corlor <- edges$color
    
 
    # net_mu = graph_from_data_frame(d=edges %>% filter(Event=="Mutually_Exclusive"),
    #                                vertices=nodes_mu, directed=F)  
    # tg_mu <- tidygraph::as_tbl_graph(net_mu) %>% activate(nodes) %>% 
    #   mutate(label=label,
    #          name=as.numeric(name),
    #          community_id=as.numeric(community_id))
    # Create an edges dataframe from the significant pairwise co-occurrences.
    # Get nodes with co-occurrences and Mutually_Exclusive relationship relatively.
    # Create an edges dataframe from the significant pairwise co-occurrences.
    edges <- data.frame(from = from_id, to = to_id,
                        color = ifelse( gene_edges$Event =="Co_Occurence", "#4058a3",  "#f18826")) %>% cbind(gene_edges ) 
    
    co_nodes = c(t(edges[which(edges$Event=="Co_Occurence"),c("gene1","gene2")]))  %>% unique(.)
    mu_nodes = c(t(edges[which(edges$Event=="Mutually_Exclusive"),c("gene1","gene2")]))  %>% unique(.)
    
    #nodes_mu = nodes %>% filter(label %in% mu_nodes) 
    if (length(co_nodes)>0) {
      nodes_co = nodes %>% filter(label %in% co_nodes) 
      edges_co = edges %>% filter(Event=="Co_Occurence",gene1 %in% co_nodes, gene2 %in% co_nodes)
      
      net_co = graph_from_data_frame(d=edges_co, vertices=nodes_co, directed=F) 
      
      # Calculate groups for nodes with co-occurrences relationship
      tg_co <- tidygraph::as_tbl_graph(net_co) %>% activate(nodes) %>% mutate(label=label,name=as.numeric(name)) %>% mutate(community_id = group_components()) 
      
      # Merge group of co-occurrence nodes label to whole nodes 
      gene_nodes = left_join(nodes,data.frame(label=V(tg_co)$label,community_id=V(tg_co)$community_id),by="label") %>%
        mutate(community_id=as.factor(ifelse(is.na(community_id),0,community_id)))
      
      net= graph_from_data_frame(d=edges,vertices=gene_nodes, directed=F)  
      
      tg <- tidygraph::as_tbl_graph(net) %>% activate(nodes) %>% mutate(label=label)
      l = net_to_layout(net,gene_nodes)
    
      if (!is.na(output)) {
        ggraph_plot(tg,l,fill=TRUE) 
        ggsave(output,width = 12, height = 10)} else
          ggraph_plot(tg,l,fill=TRUE)
      
      }
    } else{
      
      gene_nodes=nodes
      net= graph_from_data_frame(d=edges,vertices=gene_nodes, directed=F)  
      tg <- tidygraph::as_tbl_graph(net) %>% activate(nodes) %>% mutate(label=label)
      l=cbind(get_layers(net),layout_with_sugiyama(net, layers=V(net)$value)$layout[,1]) 
      
      if (!is.na(output)) {
        ggraph_plot(tg,l) 
        ggsave(output,width = 12, height = 10)} else
          ggraph_plot(tg,l)
      
      }
      
    }
    
  
    
    
    #l_co = net_to_layout(net_co,nodes_co) 
    #l_mu = net_to_layout(net_mu,nodes_mu) 
  
    ## Cluster nodes
  
    #par(mfrow=c(1,1), mar=c(2,1,2,1))
 
    ######## igraph
    # pdf(file=paste0(types,"_igrpah_network_legend.pdf"),width=15,height = 12)
    # plot(net, edge.arrow.size=.4,
    #      vertex.label=V(net)$label,
    #      vertex.label.color="white",
    #      vertex.label.cex=0.6,
    #      vertex.size = V(net)$size,
    #      #edge.color="orange",
    #      edge.curved=.1,
    #      #layout=layout_with_lgl,
    #      layout=l
    # )
    # legend("left", legend = pathway_col$Pathway, pch=21,
    #        col=pathway_col$col, pt.bg=pathway_col$col, pt.cex=3, cex=1, bty="n", ncol=1)
    # legend("bottom", legend = c("Co-Occurence","Mutually Exclusive"), lty=rep(1,2),
    #        col=c("#4058a3","#f18826"), pt.bg=c("#4058a3","#f18826"), pt.cex=2, cex=1, bty="n", ncol=1)
    # legend("right", legend = paste0(round(seq(0.1,1,length.out=10),2)),pch=21,pt.bg="grey90",pt.cex=seq(4,10,length.out=10), bty="n", ncol=1)
    # 
    # dev.off()

    
    

