#' Correlation calculation function 
#' @name cor.table
#' @param df exposure file merged with measures of interest
#' @param va Column indexs of variables A in df
#' @param vb Column indexs of variables B in df
#' @return A table including correlation r,p,n between variables A and variables B
#' @import dplyr
#' @import reshape2
#' @export
cor.table = function(df,var_idx,sig_idx,delete_0=TRUE,min_n=30){
  
  if (delete_0) df[,c(sig_idx,var_idx)] <- apply(df[,c(sig_idx,var_idx)],2,function(x) replace(x,x==0,NA)) # avoid 0 and zero-deviation

  for (i in var_idx)
    for (j in sig_idx) {
        summary = rbind(summary,c(i,j,length(which(!apply(is.na(df[,c(i,j)]),1,any)))))
    }
  colnames(summary) = c("var_idx","sig_idx","n")

  summary = as.data.frame(summary[-1,]) %>% dcast(var_idx~sig_idx,value.var = "n") %>%
    mutate(min=apply(.[,-1],1,min)) %>%
    filter(min>min_n)
  var_idx= summary$var_idx

  res <- psych::corr.test(df[,c(var_idx)],df[,c(sig_idx)],method = "spearman",use = "pairwise",adjust="BH",ci=FALSE)
  if (!is.null(res)) {
    r = melt(res['r']) %>% dplyr::rename(r = value) %>% dplyr::select(-L1)
    p = melt(res['p']) %>% dplyr::rename(adj.p = value) %>% dplyr::select(-L1)
    n = melt(res['n']) %>% dplyr::rename(n = value) %>% dplyr::select(-L1)

    cor_table = left_join(r,p,by=c("Var1","Var2"))
    if (ncol(n)>1) {
      cor_table %>%
        left_join(n,by=c("Var1","Var2")) %>%
        mutate(n=as.numeric(n),significant=ifelse(adj.p<0.05,1,0))
    } else {
      cor_table %>%
        mutate(n=as.numeric(n)) %>%
        mutate(significant=ifelse(adj.p<0.05,1,0))
    }
  }
}

#' Plot scatterplot for specific measure with evo sig and subtype - 20200915
#' @name scatter_facet
#' @param df exposure file merged with measures of interest
#' @param sig Column indexs of evo signature
#' @param vb Column indexs of all measures in association heatmap in df
#' @param facet name of facet variable consistent with association heatmap
#' @return Scatterplots
#' @import dplyr
#' @import ggpubr
#' @import scales
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
#' @example scatter_facet(df=file,vb=var_index,va=sig_idx,var=x,facet_name="cancertype",x_lab="Signature Exposure (count / MB)",legend.position = "",highlight_color="lightsalmon")),ncol=1))
scatter_facet <- function(var,df,vb,var_idx,facet_name="cancertype",cor_table,label.x=NA,label.y=NA,bin.width.x=NA,bin.width.y=NA,
                          x_lab="",legend.position="right",highlight_color="#1A9993FF",zero_delete=TRUE,subtype_panel=FALSE,varwidth=FALSE,log2_scale=FALSE,
                          shape_subtype=FALSE,compare_subtype=FALSE,side_boxplot=TRUE) {
  sig_var = colnames(df)[vb]
  df=  df[,c(var_idx,vb,which(colnames(df) %in% c("subtype",facet_name)))]
  var_idx = which(colnames(df)==var)
  colnames(df)[which(colnames(df)==facet_name)]="facet"

  # if (is.na!(cor_table)) {
  #   cor_table = plyr::ddply(df, .(facet) , .fun =cor.table,va=var_idx,vb=vb) 
  # }
  
  cor_table = cor_table %>%
    mutate(p_label=paste0("r =",round(r,2),"\n p = ",round(adj.p,3),"\n n=",n),significant=as.factor(significant),Var2=as.character(Var2),
           group=as.character(group)) %>%
    filter(Var1==var) 
  colnames(cor_table)[1] ="facet"  

  df$Var=df[,var_idx]
 
  df_new <- df %>%
    reshape2::melt(id=colnames(df)[-which(colnames(df) %in% sig_var)]) %>%
    dplyr::rename(Var2=variable) %>%
    drop_na("Var") %>%
    mutate(Var=as.numeric(Var),facet=as.character(facet),Var2=as.character(Var2)) %>%
    left_join(cor_table,by=c("Var2","facet")) 
    #left_join(subset(cor_table,facet=="All") %>% 
    #            dplyr::rename(all_p_label= p_label,all_significant=significant) %>% .[,c(2,7:9)],by="Var1")
 
  if (zero_delete) df_new <- subset(df_new,Var!=0 & value!=0)
  
  if (log2_scale) 
    df_new = df_new %>%
      mutate(value=ifelse(value>0,log2(value),value),
             Var=ifelse(Var>0,log2(Var),Var))
  
  fun_median_y <- function(x){
    return(data.frame(y=median(x),label=round(median(x,na.rm=T),2)))}
  
  xmax = max(df_new$value); xmin=min(df_new$value);xmean=mean(df_new$value)
  ymin = min(df_new$Var); ymax=max(df_new$Var);ymean <- mean(df_new$Var)
  
  if (is.na(label.x)) label.x=xmax*0.75
  if (is.na(label.y)) label.y=ymax
  if (is.na(bin.width.x)) bin.width.x= (ymax-ymin)*0.05
  if (is.na(bin.width.y)) bin.width.y= (xmax-xmin)*0.05

  scatter_plot_function = function(df) {
    wilcox =  compare_means(Var ~ subtype,subset(df,subtype %in% c("ESig3","ESig4")),group.by='facet',p.adjust.method = "BH")
    df= left_join(df,wilcox,"facet")
    
    if (shape_subtype) 
      p = ggplot(df)+ geom_point(aes(x=value,y=Var,colour = significant,shape=shape),size=0.8)
    else
      p = ggplot(df)+ geom_point(aes(x=value,y=Var,colour = significant),size=0.8)   # y_axis boxplot
    
    if (side_boxplot) {
      if (compare_subtype) 
        p=p+
        geom_boxplot(data=subset(df,subtype %in% c("ESig3","ESig4")),aes(x=xmin,y=Var,fill=subtype),width=bin.width.y,varwidth=varwidth,outlier.size=0.5)+
        geom_boxplot(data=subset(df),aes(x=value,y=ymin),orientation="y",width=bin.width.x,outlier.size=0.5) # x_axis boxplot
      else
        p=p+
          geom_boxplot(data=df,aes(x=xmin,y=Var,fill=subtype),width=bin.width.y,varwidth=varwidth,outlier.size=0.5)+
          geom_boxplot(data=subset(df),aes(x=value,y=ymin),orientation="y",width=bin.width.x,outlier.size=0.5) # x_axis boxplot
    }
    p +
      scale_colour_manual(values=c("grey",highlight_color))+
      geom_smooth(method = "lm",aes(value,Var,color=significant))+
      #scale_x_continuous(trans='sqrt',breaks=trans_breaks("sqrt",function(x) x^2),labels=trans_format("sqrt",function(x) x^2))+
      geom_text(aes(x=xmin,y=ymax,label=p.signif),position = "dodge",size=3,vjust = "inward",size=3,check_overlap = TRUE)+ # Add check_overlap, otherwise produce low print quality
      facet_grid(cols=vars(facet),rows=vars(Var2),scale="free")+
      labs(y=var,x=x_lab,title=var)+
      theme_pubclean()+
      geom_text(aes(y=label.y,x=label.x,label=p_label),position = "dodge",size=3,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
      theme(legend.position = legend.position,plot.margin = unit( c(0.2,0,0,0) , units = "lines" ) )+
      scale_fill_manual(values=c("#eccbae","#abddde"))
  }
 
  if (xmin>=0) df_new = df_new %>% filter(value>0)
  if (ymin>=0) df_new = df_new %>% filter(Var>0)
 
  p1 <-   df_new %>% scatter_plot_function()
  
  if (log2_scale) p1 = p1+labs(y=paste0("log2(",var,")"),x=paste0("log2(",x_lab,")"))
  # if (all(df_new$Var>=0) & log2_scale) p1 <- p1+scale_y_continuous(trans='sqrt',breaks=trans_breaks("sqrt",function(x) x^2),labels=trans_format("sqrt",function(x) x^2))
  
 
  return(p1)

}
#' Correlation calculation function by group and plot association heatmap - 20200915
#' @name cor_facet
#' @param df exposure file merged with measures of interest
#' @param va Column indexs of variables A in df
#' @param vb Column indexs of variables B in df
#' @param facet name of group variable
#' @param heatmap specify whether plot association heatmap
#' @param title Title for association heatmap
#' @param empty_row_delete whether to delete rows without any significant correlation value
#' @param flip whether to flip facet and x
#' @param keep_all whether to show all facet
#' @param col_low color for correlation r=-1
#' @param col_high color for correlation r=1
#' @return A table including correlation r,p,n between variables A and variables B during group variable and association heatmap (optional)
#' @import dplyr
#' @importFrom plyr ddply
#' @import reshape2
#' @import ggpubr
#' @import ggplot2
#' @export
#' @Example va:var vb:sig

cor_facet = function(df,var_idx,sig_idx,facet,heatmap=FALSE,title="",empty_row_delete=FALSE,flip=FALSE,keep_all=TRUE,col_low="#4a7b94",col_high="#bb5a39",scatter=FALSE,scatter_var=NA,
                     label.x=NA,label.y=NA,bin.width.x=NA,bin.width.y=NA){

  var_va = colnames(df)[var_idx];var_vb = colnames(df)[sig_idx]
  df = df[,c(var_idx,sig_idx,which(colnames(df) %in% c(facet,"subtype","samplename")))]
  colnames(df)[which(colnames(df)==facet)] = "group"
  var_idx = which(colnames(df) %in% var_va)
  sig_idx = which(colnames(df) %in% var_vb)
  
  sig_idx = sig_idx[which(apply(df[,sig_idx],2,function(x) length(unique(x))==1)==FALSE)]
  var_idx = var_idx[which(apply(df[,var_idx],2,function(x) length(unique(x))==1)==FALSE)]
  df[,c(var_idx,sig_idx)] = apply(df[,c(var_idx,sig_idx)],2,as.numeric)
  
  if (length(sig_idx)!=0) {
    
    cor_table_all <- cor.table(df,var_idx=var_idx,sig_idx=sig_idx) %>% mutate(group='All')
 
    cor_table <- plyr::ddply(df, .(group) , .fun =cor.table,var_idx=var_idx,sig_idx=sig_idx) %>% 
      rbind(cor_table_all)  %>%
      mutate(r=ifelse(n<30,0,r)) %>%
      mutate(r=replace(r,is.na(r),0),adj.p=replace(adj.p,is.na(adj.p)|is.nan(adj.p),1)) %>% # set correlation with na or p>0.05 invisible
      mutate(r=round(r,2),Var2=as.character(Var2),label=paste0(group," \n (n=",n,")"))
    
    if (scatter) {
      p_scatter = scatter_facet(df=df,va=var_idx,vb=sig_idx,var=scatter_var,facet_name="group",highlight_color="lightsalmon",
                  label.x=label.x,label.y=label.y,bin.width.x=bin.width.x,bin.width.y=bin.width.y,subtype_panel=TRUE)
      
      # p_esig3 = scatter_facet(df= df %>% filter(subtype=="ESig3"),va=va,vb=vb,var=scatter_var,facet_name="group",highlight_color="lightsalmon",
      #                         label.x=label.x,label.y=label.y,bin.width.x=bin.width.x,bin.width.y=bin.width.y,subtype_panel=TRUE)
      # p_esig4 = scatter_facet(df=df %>% filter(subtype=="ESig4"),va=va,vb=vb,var=scatter_var,facet_name="group",highlight_color="lightsalmon",
      #                         label.x=label.x,label.y=label.y,bin.width.x=bin.width.x,bin.width.y=bin.width.y,subtype_panel=TRUE)
      # p_subtype = p_esig3/p_esig4
    }
    
    
    if (heatmap) {
        cor_table_plot <- cor_table %>% mutate(r=replace(r,adj.p>0.05,0))  %>% 
          left_join( cor_table %>% group_by(Var1,group) %>% dplyr::summarise(n_show=min(n)),by=c("Var1","group")) 
    
    # delete empty rows
        if (keep_all==FALSE)  cor_table_plot <- subset( cor_table_plot,group!="All")
        
        if (empty_row_delete==TRUE) {
          non_empty_gene <- cor_table_plot %>% group_by(Var2) %>% 
            dplyr::summarise(empty_r=all(r==0)) %>% filter(empty_r==FALSE) %>% as.data.frame() %>% .[,1] %>% as.character()
          
          cor_table_plot <- subset(cor_table_plot,Var2 %in% non_empty_gene) 
        }
    
        if (flip==TRUE) {
          
          n_facet <- length(non_empty_gene)
          
          p_heatmap <- cor_table_plot %>% 
            ggplot(aes(y=Var1,x=group,fill=r)) +
            geom_tile()+ geom_text(aes(y=Var1,x=group,label=r),size=4,col='#ffffff')+
            facet_grid(cols=vars(Var2),switch="y",scale="free",space="free")
          
        } else {
          
            p_heatmap <- cor_table_plot %>%
              ggplot(aes(y=Var1,x=Var2,fill=r)) +
              geom_tile()+ 
              geom_text(aes(y=Var1,x=Var2,label=r),size=4,col='#ffffff')
          
              if (length(unique(cor_table_plot$n))<=5) {
                
                p_heatmap <- p_heatmap + facet_grid(cols=vars(label),switch="y",scale="free",space="free")
                
              } else {
                
                n_var1 <- length(unique(cor_table_plot$Var2))
                p_heatmap <- p_heatmap + 
                  geom_text(aes(y=Var1,x=n_var1+1,label=n_show),size=3,col="grey50")+
                  facet_grid(cols=vars(group),switch="y",scale="free",space="free")+ 
                  coord_cartesian( xlim=c(1,n_var1+0.5),clip = "off")
              }
        }
  
       
        
    p_heatmap <- p_heatmap +
      scale_fill_gradient2(low=col_low,mid="#ffffff",high=col_high,midpoint = 0,name="Spearman Correlation")+
      labs(x="",y="",title=title)+ theme_pubclean()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),plot.margin = unit(c(1, 7, 1, 1), "lines"),legend.position = "right")
    
    } 
    
    # if (heatmap & scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_heatmap,p_scatter,p_subtype,non_empty_gene)),return(list(cor_table,p_heatmap,p_scatter,p_subtype)))
    # if (!heatmap & scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_scatter,p_subtype,non_empty_gene)),return(list(cor_table,p_scatter,p_subtype)))
    # if (heatmap & !scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_heatmap,non_empty_gene)),return(list(cor_table,p_heatmap)))
    # if (!heatmap & !scatter) return(cor_table)
    
    if (heatmap & scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_heatmap,p_scatter,non_empty_gene)),return(list(cor_table,p_heatmap,p_scatter)))
    if (!heatmap & scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_scatter,non_empty_gene)),return(list(cor_table,p_scatter)))
    if (heatmap & !scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_heatmap,non_empty_gene)),return(list(cor_table,p_heatmap)))
    if (!heatmap & !scatter) return(cor_table)
    
    
    
    
  } else {print("the standard deviation is zero")}
}




#' Load multiple format ccf file
#' @name load_ccf
#' @param samplename cancer type
#' @param input 
#' @export
#' @return ssm
load_ccf <- function(samplename,input){
  Check <- ArgumentCheck::newArgCheck()
  suppressWarnings(rm(ssm,res,ccubeRes))
  
  format <- paste0(input,samplename,"/ccubeRes.RData")
  
  if (file.exists(format )) load_ssm = get(load(format4)) 
  return(load_ssm$ssm) 
}

#' create multiple dir
#' @name multi.dir.create
#' @param list list of directory
#' @return create multiple directory
#' @export
multi_dir_create <- function(dirlist){
  for (i in dirlist) {if (!dir.exists(i)) dir.create(i,recursive = T)}
}

#' Unify format of data frame 
#' @name file_format
#' @param filename data frame
#' @param samplenamecol column index of samplename
#' @export
file_format <- function(filename=filename,samplenamecol){
    names <- colnames(filename)
    names[samplenamecol] <- "samplename"
    names -> colnames(filename)
    filename$samplename <- substr(filename$samplename,1,12)
    filename$samplename <- gsub("[.]","-",filename$samplename)
    return(filename)
  }

#' Merge ssm and cna files
#' @name ParseSnvCnaPcawgFormat
#' @param ssm ssm
#' @param cna cna
#' @export
#' @import dplyr
ParseSnvCnaPcawgFormat <- function (ssm, cna) {
    
    ssm <- ssm %>%
      mutate(chr = substr(chr,4,length(chr)),
             cn_frac = NA,
             major_cn = NA,
             minor_cn = NA,
             mutation_id = NA)
    
    for (jj in seq_len(nrow(cna)) ) {
      cc = cna[jj,]
      
      idx = which(ssm$chr == cc$chromosome &  (ssm$Start_Position >= cc$start & ssm$End_Position <= cc$end) )
      
      if (length(idx) > 0) {
        ssm[idx,] <- ssm[idx,] %>% 
          mutate( major_cn=cc$major_cn,minor_cn =cc$minor_cn, cn_frac = 1)
      }
    }
    
    ssm$mutation_id = paste0(ssm$chr, "_", ssm$Start_Position )
    
    ssm <- ssm %>%
      select(-chr,-Start_Position,-End_Position,-df.n_alt_count,-n_ref_count) %>%
      rename(var_counts='t_alt_count',ref_counts='t_ref_count') %>%
      mutate(total_counts=var_counts+ref_counts,normal_cn=2) %>%
      filter(!is.na(major_cn) ,!is.na(minor_cn),!is.na(cn_frac),major_cn > 0)
    
    return(ssm)
}

#' Customize top strip color
#' @name strip_color
#' @param p plot
#' @param col customized color
#' @param draw whether to display in plots panal
#' @param direction strip location
#' @import ggplot2 
#' @import grid
#' @export
strip_color <- function(p,col1=signature_col,col2=NULL,draw=FALSE,direction='top'){
  p1 <- ggplot_gtable(ggplot_build(p))
  k <- 1
  
  if (direction=='top') strip_col <- which(grepl('strip-t', p1$layout$name))
  if (direction=='bottom') strip_col <- which(grepl('strip-b', p1$layout$name))
  if (direction=='left') strip_col <- which(grepl('strip-l', p1$layout$name))
  if (direction=='right') strip_col <- which(grepl('strip-r', p1$layout$name))
  if (direction=='both'){
    strip_col_t <- which(grepl('strip-t', p1$layout$name))
    strip_col_l <- which(grepl('strip-l', p1$layout$name))
    
    k <- 1
    for (j in strip_col_t) {
      j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
      p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- col1[k]
      k <- k+1
    }
    
    k <- 1
    for (j in strip_col_l) {
      j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
      p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- col2[k]
      k <- k+1
    }
  } else {
    k <- 1
    for (j in strip_col) {
      j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
      p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- col1[k]
      k <- k+1
    }
  }
  
  if (draw) grid.draw(p1)
  return(p1)
}


