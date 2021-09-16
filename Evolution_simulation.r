#### Get mutation label
mut_t = function(t_init,t,lambda=log(2),mu=10,max=10000){
  t_next = function(t_now,t_sub_init,lambda,mu){t_now+1/lambda*log(lambda/mu/2/n_cell_s(t_sub_init,t_now,lambda)+1)}
  t_now = t_init;t_idx = t_init
  while (t_now<t & length(t_idx)<= max) {
    t_now = t_next(t_now,t_init,lambda=lambda,mu=mu)
    t_idx = c(t_idx,t_now)
  }
  t_idx[t_idx<=t & t_idx>t_init]
}
#### Cell growth
n_cell_s = function(t_init,t,lambda,type="exp"){
  exp(lambda*(t-t_init))
}
#### Beta-binomial - 对于reads来说不够大
betabinom <- function(p,n,ro=0){
  shape1 = p*((1/ro)-1)
  shape2 = (1-p)*((1/ro)-1)
  rbinom(1,as.integer(n),rbeta(1, shape1, shape2)) # p_binomial = rbeta(n, shape1, shape2)
}
####    
nested_2sub_mut <- function(t0=1,t1,t2,t,l0=log(2),s1,s2,mu=10,n_clonal=10){
  # subclonal mutations during tumor neutral grow between t0=0 to t1
  mut_sub_n=data.table(mut_t=mut_t(0,t1,l0)) 
  mut_sub_n[,":="(mut_idx=(1:.N+n_clonal),
                  subclone="n",type="n-passenger",s=0,growth_rate=l0,n_cells=floor(n_cell_s(mut_sub_n$mut_t,t1,l0)))][]
  
  ##sample subclonal mutations turn into s1 clonal mutations       
  sample_idx_n = sample(1:nrow(mut_sub_n),floor(mu*2*log(2)*t1))
  s1_clonal_mut_id = mut_sub_n[sample_idx_n,mut_idx]
  mut_clonal_s1 = mut_sub_n[sample_idx_n,]
  mut_clonal_s1[,':='(subclone="s1",type="s1-clonal",n_cells=1,s=s1,growth_rate=(1+s1)*l0)]
  mut_sub_n[sample_idx_n,":="(n_cells = n_cells-1, type="s1-clonal")]
  mut_0_t1 = rbind(mut_sub_n,mut_clonal_s1) 
  
  #tumor grow during time t1 to t2
  mut_sub_s1=data.table(mut_t=mut_t(t1,t2,(1+s1)*l0)) 
  mut_sub_s1[,":="(mut_idx=1:.N+max(mut_0_t1$mut_idx),
                   subclone="s1",type="s1-passenger",s=s1,growth_rate=(1+s1)*l0,n_cells=floor(n_cell_s(mut_t,t2,(1+s1)*l0)))]
  
  sample_idx_s1 = sample(1:nrow(mut_sub_s1),floor((2*log(2)*(s1+2)*(t2-t1))*mu)) # attention of this 
  s2_clonal_mut_id = mut_sub_s1[sample_idx_s1,]
  mut_0_t1[,n_cells:=n_cells*n_cell_s(t1,t2,growth_rate)]
  mut_0_t1[mut_idx %in% s1_clonal_mut_id & subclone=="s1",':='(n_cells=n_cells-1,type="s2-clonal")]
  mut_clonal_s2 = filter(rbind(mut_sub_s1[sample_idx_s1,], mut_0_t1[mut_idx %in% s1_clonal_mut_id & subclone=="s1",]), n_cells>=1)
  mut_clonal_s2[,':='(subclone="s2",type="s2-clonal",n_cells=1,s=s2,growth_rate=(1+s2)*l0)]
  mut_sub_s1[sample_idx_s1,":="(n_cells = n_cells-1, type="s2-clonal")]
  mut_t1_t2 = filter(rbind(mut_sub_s1,mut_clonal_s2) ,n_cells>=1)
  mut_0_t2= filter(rbind(mut_0_t1, mut_t1_t2),n_cells>=1)
  
  # during t2 to t
  mut_sub_s2=data.table(mut_t=mut_t(t2,t,(1+s2)*l0)) 
  mut_sub_s2[,":="(mut_idx=1:.N+max(mut_0_t2$mut_idx),
                   subclone="s2",type="s2-passenger",s=s2,
                   growth_rate=(1+s2)*l0,n_cells=floor(n_cell_s(mut_t,t,(1+s2)*l0)))][]
  mut_0_t = copy(mut_0_t2)
  mut_0_t = rbind(mut_0_t[,n_cells:=n_cells*n_cell_s(t2,t,growth_rate)],mut_sub_s2) %>%
    mutate(type=factor(type,levels=c("s2-clonal","s2-passenger","s1-clonal","s1-passenger","n-passenger","clonal"))) %>% arrange(mut_idx,type)
  mut_annotation = mut_0_t  %>% distinct(mut_idx,.keep_all = TRUE)
  
  
  n_cells_0_t_all = n_cell_s(0,t,l0)+n_cell_s(t1,t,(1+s1)*l0)+n_cell_s(t2,t,(1+s2)*l0)
  
  mut_0_t_all = mut_0_t %>%
    group_by(mut_idx) %>%
    dplyr::summarise(n_cells_mut_all=floor(sum(n_cells))) %>%
    mutate(ccf=n_cells_mut_all/n_cells_0_t_all) %>%
    left_join(mut_annotation,by="mut_idx") %>% 
    rbind(data.frame(mut_idx=1:n_clonal,n_cells_mut_all=n_cells_0_t_all,ccf=1,mut_t=0,
                     subclone="s2",type="clonal",s=s2,
                     growth_rate=(1+s2)*l0,n_cells=n_cells_0_t_all))
  # sythetic data algorithms
  mut_0_t_all %>%
    as.data.frame() %>%
    filter(ccf>=0.01) %>%
    mutate(n_cells_all_sample = rpois(nrow(.),n_cells_0_t_all)) %>%
    mutate(mut_sample_cells=apply(.,1,function(x) betabinom(p=as.numeric(x[3]),n=as.integer(x[2]),ro=0.2))) %>% #p=ccf,n=n_cells_all_sample
    mutate(ccf_sample=mut_sample_cells/n_cells_all_sample) 
}

### v2.0 - inclde t<t1 or t2

nested_2sub_mut <- function(t0=1,t1,t2,t,l0=log(2),s1,s2,mu=10,n_clonal=10){
  # subclonal mutations during tumor neutral grow between t0=0 to t1
  ##sample subclonal mutations turn into s1 clonal mutations 
  if (t<t1) {
    mut_sub_n=data.table(mut_t=mut_t(0,t,l0)) 
    mut_sub_n[,":="(mut_idx=(1:.N+n_clonal),
                    subclone="n",type="n-passenger",s=0,growth_rate=l0,n_cells=floor(n_cell_s(mut_sub_n$mut_t,t,l0)))][]
    mut_0_t = mut_sub_n
  } else{
    mut_sub_n=data.table(mut_t=mut_t(0,t1,l0)) 
    mut_sub_n[,":="(mut_idx=(1:.N+n_clonal),
                    subclone="n",type="n-passenger",s=0,growth_rate=l0,n_cells=floor(n_cell_s(mut_sub_n$mut_t,t1,l0)))][]
    sample_idx_n = sample(1:nrow(mut_sub_n),floor(mu*2*log(2)*t1))
    s1_clonal_mut_id = mut_sub_n[sample_idx_n,mut_idx]
    mut_clonal_s1 = mut_sub_n[sample_idx_n,]
    mut_clonal_s1[,':='(subclone="s1",type="s1-clonal",n_cells=1,s=s1,growth_rate=(1+s1)*l0)]
    mut_sub_n[sample_idx_n,":="(n_cells = n_cells-1, type="s1-clonal")]
    mut_0_t1 = rbind(mut_sub_n,mut_clonal_s1) 
    n_cells_0_t_all = floor(n_cell_s(t0,t1,l0))
    
    if(t<t2){
      mut_sub_s1=data.table(mut_t=mut_t(t1,t,(1+s1)*l0)) 
      mut_sub_s1[,":="(mut_idx=1:.N+max(mut_0_t1$mut_idx),
                       subclone="s1",type="s1-passenger",s=s1,growth_rate=(1+s1)*l0,n_cells=floor(n_cell_s(mut_t,t,(1+s1)*l0)))]
      mut_0_t1[,n_cells:=floor(n_cells*n_cell_s(t1,t,growth_rate))]
      mut_0_t = filter(rbind(mut_sub_s1, mut_0_t1[mut_idx %in% s1_clonal_mut_id & subclone=="s1",]), n_cells>=1)
      n_t1=floor(n_cell_s(t0,t1,l0))
      n_cells_0_t_all = (n_t1-1)*n_cell_s(t1,t,l0) + n_cell_s(t1,t,(1+s1)*l0)
    } else {
      #tumor grow during time t1 to t2
      mut_sub_s1=data.table(mut_t=mut_t(t1,t2,(1+s1)*l0)) 
      mut_sub_s1[,":="(mut_idx=1:.N+max(mut_0_t1$mut_idx),
                       subclone="s1",type="s1-passenger",s=s1,growth_rate=(1+s1)*l0,n_cells=floor(n_cell_s(mut_t,t2,(1+s1)*l0)))]
      sample_idx_s1 = sample(1:nrow(mut_sub_s1),floor((2*log(2)*(s1+2)*(t2-t1))*mu)) # attention of this 
      s2_clonal_mut_id = mut_sub_s1[sample_idx_s1,]
      mut_0_t1[,n_cells:=n_cells*n_cell_s(t1,t2,growth_rate)]
      mut_0_t1[mut_idx %in% s1_clonal_mut_id & subclone=="s1",':='(n_cells=n_cells-1,type="s2-clonal")]
      mut_clonal_s2 = filter(rbind(mut_sub_s1[sample_idx_s1,], mut_0_t1[mut_idx %in% s1_clonal_mut_id & subclone=="s1",]), n_cells>=1)
      mut_clonal_s2[,':='(subclone="s2",type="s2-clonal",n_cells=1,s=s2,growth_rate=(1+s2)*l0)]
      mut_sub_s1[sample_idx_s1,":="(n_cells = n_cells-1, type="s2-clonal")]
      mut_t1_t2 = filter(rbind(mut_sub_s1,mut_clonal_s2) ,n_cells>=1)
      mut_0_t2= filter(rbind(mut_0_t1, mut_t1_t2),n_cells>=1)
      
      # during t2 to t
      mut_sub_s2=data.table(mut_t=mut_t(t2,t,(1+s2)*l0)) 
      mut_sub_s2[,":="(mut_idx=1:.N+max(mut_0_t2$mut_idx),
                       subclone="s2",type="s2-passenger",s=s2,
                       growth_rate=(1+s2)*l0,n_cells=floor(n_cell_s(mut_t,t,(1+s2)*l0)))][]
      mut_0_t = copy(mut_0_t2)
      mut_0_t = rbind(mut_0_t[,n_cells:=n_cells*n_cell_s(t2,t,growth_rate)],mut_sub_s2) %>%
        mutate(type=factor(type,levels=c("s2-clonal","s2-passenger","s1-clonal","s1-passenger","n-passenger","clonal"))) %>% arrange(mut_idx,type)
      
      n_t1=floor(n_cell_s(t0,t1,l0))
      n_t2=floor((n_cell_s(t0,t1,l0)-1)*n_cell_s(t1,t2,l0)+n_cell_s(t1,t2,(1+s1)*l0))
      n_neutral = (n_t1-1)*n_cell_s(t1,t,l0)
      n_s1 =  (n_cell_s(t1,t2,(1+s1)*l0)-1)*n_cell_s(t2,t,(1+s1)*l0)
      n_s2 =   n_cell_s(t2,t,(1+s2)*l0)
      n_cells_0_t_all = n_neutral + n_s1+n_s2
    }
    
  }
  mut_annotation = mut_0_t  %>% distinct(mut_idx,.keep_all = TRUE)
  
  mut_0_t_all = mut_0_t %>%
    group_by(mut_idx) %>%
    dplyr::summarise(n_cells_mut_all=floor(sum(n_cells))) %>%
    mutate(ccf=n_cells_mut_all/n_cells_0_t_all) %>%
    left_join(mut_annotation,by="mut_idx") 
  # sythetic data algorithms
  mut_0_t_all %>%
    as.data.frame() %>%
    filter(ccf>=0.01) %>%
    mutate(n_cells_all_sample = rpois(nrow(.),n_cells_0_t_all)) %>%
    mutate(mut_sample_cells=apply(.,1,function(x) betabinom(p=as.numeric(x[3]),n=as.integer(x[2]),ro=0.2))) %>% #p=ccf,n=n_cells_all_sample
    mutate(ccf_sample=mut_sample_cells/n_cells_all_sample) %>% 
    rbind(data.frame(mut_idx=1:n_clonal,n_cells_mut_all=n_cells_0_t_all,ccf=1,mut_t=0,
                     subclone="s2",type="clonal",s=s2,
                     growth_rate=(1+s2)*l0,n_cells=n_cells_0_t_all,n_cells_all_sample=n_cells_0_t_all,mut_sample_cells=n_cells_0_t_all,ccf_sample=1))
}

#### Main ###################################
library(plyr)
library(dplyr)
library(ggpubr)
library(data.table)
library(gganimate)
library(ggplot2)
install.packages("gganimate")

type_col=c("#008000","#00B4E0","#330080","#f7931e","#F6CBCB","#d40000")
names(type_col) = c("n-passenger","s1-passenger","s2-passenger","s2-clonal","s1-clonal","clonal")
### plot - gif
t_idx = data.frame(tend=c(10:20)) 
gif = ddply(t_idx,.(tend), .fun= function(x) nested_2sub_mut(t0=0,t1=5,t2=7,t=x$t,s1=0.2,s2=0.8)   )

res = gif  %>% 
  filter(ccf_sample>=0.1) %>%
  ggplot()+geom_histogram(aes(x=ccf_sample,fill=type))+
  theme_pubr()+
  theme(legend.position = "right")+
  scale_fill_manual(values=type_col)+
  labs(x="Cancer Cell Fraction",caption ="t1=5,t2=7 \n s1=0.2,s2=0.8 \n n_clonal=10,mu=10")   +
  transition_time(tend)+
  labs(title = "tend: {frame_time}")
animate(res)
### single plot
nested_2sub_mut(t0=0,t1=5,t2=7,s1=0.2,s2=0.8,t=15) %>%
  filter(ccf_sample>=0.1) %>%
  ggplot()+geom_histogram(aes(x=ccf_sample,fill=type))+
  theme_pubr()+
  theme(legend.position = "right")+
  scale_fill_manual(values=type_col)+
  labs(x="Cancer Cell Fraction",caption ="t1=5;t2=6;s1=0.2;s2=0.8;l0=log(2);number=10;mu=10;t=15")


### COADREAD
### We get the probability of observing certain mutation in the sample after selection (Fitness)
install.packages("TCGAmutations")
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("TCGAmutations")
library(TCGAmutation)
# Load MAF
if (types != "All") {
  lstObject <- lapply(typesMaf, function(i) tcga_load(i))
} else {
  lstObject <- lapply(tcga_available()$Study_Abbreviation[1:33], function(i) tcga_load(i))
}

