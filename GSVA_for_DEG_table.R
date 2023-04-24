
# GSVA to t test 
get_gene_df=function(gmt){
    genes=list()
    gene_df = read.delim(gmt,sep='\t',header=F)
    for (i in 1:nrow(gene_df)){
        g = as.character(gene_df[i,3:ncol(gene_df)])
        g = g[g!=""]
        genes[[as.character(gene_df$V1[i])]] = g
    }
    return(genes)
}


main=function(path_gmt,gct_file,cls_file,genelist_file,outpath,outname){
    library('GSVA')
    genelist=list()
    for (c in list.files(path_gmt))
    {
    gmt_list =paste0(path_gmt,c)
    genelist = c(genelist,get_gene_df(gmt_list))
    }

    deg = read.table(gct_file, sep="\t",skip=2,header=T,row.names = 1)
    deg$DESCRIPTION=NULL

    genes_top3 = read.table(genelist_file,header=F)$V1

    class_vector = t(read.table(cls_file,sep=' ',skip=2))
    class_vector = as.character(class_vector)
    deg = deg[rownames(deg)%in%genes_top3,]

    row = rownames(deg)
    deg = as.data.frame(sapply(deg, as.numeric))
    rownames(deg) = row

    deg = as.matrix(deg)


    gsva_es <- gsva(deg, genelist,class_vector,method = 'gsva')


    gsva_es = as.data.frame(gsva_es)
    gsva_es[c('group'),]=  as.numeric(class_vector)
    gsva_es = t(gsva_es)
    gsva_es = as.data.frame(gsva_es)
    gsva_es = unique(gsva_es)


    ttest_ls=c()
    cols = colnames(gsva_es)
    cols = cols[cols!='group']
    for (i in cols){
        group1 = as.numeric(gsva_es[gsva_es$group==1,i])
        group0 = as.numeric(gsva_es[gsva_es$group==0,i])
        tt = t.test(group1,group0)
        tt = tt$p.value
        ttest_ls = c(ttest_ls,tt)
    }
    gsva_es$group=NULL
    gsva_es[c('PVAL'),] = ttest_ls
    gvsa_es = t(gsva_es)
    gvsa_es = as.data.frame(gsva_es)
    write.table(gsva_es[,c(gsva_es[c('PVAL'),]<0.05)],paste0(outpath,outname,'.0.05.tsv',collapse=''),sep='\t',quote=FALSE)
    write.table(gsva_es,paste0(outpath,outname,'.tsv',collapse=''),sep='\t',quote=FALSE)

    return(0)
}



PATH_GMT =
GCT_FILE =
CLS_FILE =
GENELIST_FILE =
OUTPATH =
OUTNAME =


main(PATH_GMT,GCT_FILE,CLS_FILE,GENELIST_FILE,OUTPATH,OUTNAME)
