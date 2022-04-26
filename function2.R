#function2
download_TCGA<-function(cancer_type,dir){
  # library(TCGAbiolinks)
  # library(dplyr)
  # library(stringr)
  # library(data.table)
  # library(DT)
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # 
  # BiocManager::install("miRBaseConverter")
  library(miRBaseConverter)
  
  # cancer_type<-'STAD'
  # dir<-'D:\\cancer_type_classification\\tmp'
  
  request_cancer=cancer_type
  setwd(dir)
  cancer_type=str_c('TCGA-',request_cancer)
  query1 <- GDCquery(project = cancer_type,
                     data.category = "Transcriptome Profiling", 
                     data.type = "Gene Expression Quantification", 
                     workflow.type = "STAR - Counts")
  GDCdownload(query1, method = "client", files.per.chunk = 25)
  expdat1 <- GDCprepare(query = query1)
  count_matrix=assay(expdat1)
  write.csv(count_matrix,file = paste(cancer_type,"Counts.csv",sep = "-"))
  
  query <- GDCquery(project = cancer_type,
                    data.category = "Transcriptome Profiling", 
                    data.type = "Isoform Expression Quantification") 
  GDCdownload(query, method = "client", files.per.chunk = 25)
  expdat <- GDCprepare(query)
  
  expdat = expdat[-grep("precursor|stemloop", expdat$miRNA_region),]
  expdat$miRNA_region = gsub(pattern = c(as.character("mature,")), replacement = "", x = expdat$miRNA_region)
  x = aggregate(read_count ~ miRNA_region + barcode, data=expdat, sum)
  x<-dcast(x, miRNA_region ~ barcode, value.var="read_count")
  miRNA_region<-miRNA_AccessionToName(x$miRNA_region)
  x$miRNA_region<-miRNA_region$TargetName
  x<-x[!(is.na(x$miRNA_region)),]
  x[is.na(x)]<-0
  
  write.csv(x,file = paste(cancer_type,"miRNA.csv",sep = "-"),row.names = F)
  
  clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")
  write.csv(clinical,file = paste(cancer_type,"clinical.csv",sep = "-"))
  
  rm(cancer_type,expdat,query,expdat1,count_matrix,query1,clinical)
  cat("Download is finished")
  
} #download tcga transcriptome data (rna expression)
##################################################################################################
download_TCGA_planB<-function(){
  setwd("D:\\tmp")
  library('TCGAbiolinks')
  project_name <- "TCGA-STAD"
  
  # Defines the query to the GDC
  query <- GDCquery(project = project_name,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    experimental.strategy = "RNA-Seq",
                    workflow.type = "STAR - Counts")
  
  # Get metadata matrix
  metadata <- query[[1]][[1]]
  
  # Download data using api
  GDCdownload(query, method = "api")
  
  # Get main directory where data is stored
  main_dir <- file.path("GDCdata", project_name)
  # Get file list of downloaded files
  file_list <- file.path("GDCdata", project_name,list.files(main_dir,recursive = TRUE)) 
  
  # Read first downloaded to get gene names
  test_tab <- read.table(file = file_list[1], sep = '\t', header = TRUE)
  # Delete header lines that don't contain usefull information
  test_tab <- test_tab[-c(1:4),]
  # STAR counts and tpm datasets
  tpm_data_frame <- data.frame(test_tab[,1])
  count_data_frame <- data.frame(test_tab[,1])
  
  # Append cycle to get the complete matrix
  for (i in c(1:length(file_list))) {
    # Read table
    test_tab <- read.table(file = file_list[i], sep = '\t', header = TRUE)
    # Delete not useful lines
    test_tab <- test_tab[-c(1:4),]
    # Column bind of tpm and counts data
    tpm_data_frame <- cbind(tpm_data_frame, test_tab[,7])
    count_data_frame <- cbind(count_data_frame, test_tab[,4])
    # Print progres from 0 to 1
    print(i/length(file_list))
  }
}
##################################################################################################
tnm_statistic<-function(dir){
  # dir<-"D:\\cancer_type_classification\\tcga_data\\BRCA\\TCGA-BRCA-clinical.csv"
  
  
  clinic<-read.csv(dir,row.names = 1,check.names = F)
  clinic<-clinic[,c('submitter_id','ajcc_pathologic_m','ajcc_pathologic_t','ajcc_pathologic_n')]
  # write.csv(clinic,dir,quote = F)
  #remove X stage
  clinic=clinic[which(str_sub(clinic$ajcc_pathologic_m,1,2)!='MX' &
                      str_sub(clinic$ajcc_pathologic_t,1,2)!='TX' &
                      str_sub(clinic$ajcc_pathologic_n,1,2)!='NX'),]
  #n-metastasis
  meta=clinic[which(str_sub(clinic$ajcc_pathologic_n,1,2)!='N0'),] 
  meta=meta[which(str_sub(meta$ajcc_pathologic_t,1,2)!='T0'),] 
  meta=meta[which(str_sub(meta$ajcc_pathologic_m,1,2)!='M1'),]
  #non-n-metastasis
  nonmeta=clinic[which(str_sub(clinic$ajcc_pathologic_n,1,2)=='N0'),] 
  nonmeta=nonmeta[which(str_sub(nonmeta$ajcc_pathologic_t,1,2)!='T0'),] 
  nonmeta=nonmeta[which(str_sub(nonmeta$ajcc_pathologic_m,1,2)!='M1'),] 
  cat(str_c("\nnum of meta: ",nrow(meta)," ,num of non-meta: ",nrow(nonmeta),"\n"))
  res<-c(nrow(meta),nrow(nonmeta))
  return(res)
}
##################################################################################################
pre_clinic<-function(dir){
  # dir<-'D:\\cancer_type_classification\\tcga_data\\ACC\\TCGA-ACC-clinical.csv'
  
  dt<-read.csv(dir,check.names = F,row.names = 1)
  dt<-dt[,c('submitter_id','ajcc_pathologic_m','ajcc_pathologic_t','ajcc_pathologic_n')]
  dt<-dt[is.na(dt$ajcc_pathologic_n)==F,]
  write.csv(dt,dir,quote = F)
}
##################################################################################################
check_csv<-function(dir){
  #dir<-"D:\\cancer_type_classification\\tcga_data\\THYM\\TCGA-THYM-Counts.csv"
  ref<-read.csv(dir)
  row<-nrow(ref)
  col<-ncol(ref)-1
  cat(str_c(dir,':\n'))
  ind<-c(row,col)
  cat(str_c('[',row,',',col,']'))
  cat('\n')
  rm(dir,ref,row,col)
  return(ind)
} # check dimension of csv file
##################################################################################################
manipul_mat<-function(dir,sep,cont,row,col){
  # dir<-'D:\\cancer_type_classification\\tcga_data\\BLCA\\mRNAMatrix2.csv'
  # cont<-'id'
  # row<-1
  # col<-1
  # sep<-','
  
  que<-read.table(dir,header = F, sep = sep)
  que[row,col]<-cont
  write.table(que,dir,row.names = F, quote = F, col.names = F, sep = sep)
}
##################################################################################################
pickcom<-function(dir, op){
  # library(data.table)
  # library(limma)
  # library(stringr)
  # cancer_type='KIRC'
  setwd(dir)
  miRNA_exp=fread('miRNAMatrix2.csv',sep=",",header=T,check.names=F)
  lncRNA_exp=fread('mRNAMatrix2.csv',sep=",",header=T,check.names=F)
  
  miRNA=data.frame(miRNA_exp,check.names=F)
  lncRNA=data.frame(lncRNA_exp,check.names=F)
  
  #mRNA=avereps(mRNA[,-1],ID=mRNA$id)#avereps? and delete column 1, make column2 as id
  miRNA=avereps(miRNA[,-1],ID=miRNA$id)
  lncRNA=avereps(lncRNA[,-1],ID=lncRNA$id)
  #cat("searching common samples.\n")
  ss=c()
  for(i in colnames(miRNA)){
    col_names = str_sub(i,start = 1, end = 15) #id must be fixed like this:TCGA-J4-A83J-11A-11R-A36G-07->TCGA-EJ-7794-11
    ss=c(ss,col_names)
  }
  colnames(miRNA)<-ss
  ss=c()
  for(i in colnames(lncRNA)){
    col_names = str_sub(i,start = 1, end = 15)
    ss=c(ss,col_names)
  }
  colnames(lncRNA)<-ss
  
  if (op==T){
    heji=intersect(colnames(lncRNA),colnames(miRNA))#find common parts (mrna-lncrna<->mirna)
    #mRNA=mRNA[,heji]
    miRNA=miRNA[,heji]
    miRNA[is.na(miRNA)]<-0
    lncRNA=lncRNA[,heji]
  }
  
  #write.csv(miRNA,'miRNAmatrix.csv')
  write.csv(lncRNA,'mRNAMatrix2.csv')
  write.csv(miRNA,'miRNAMatrix2.csv')
}
##################################################################################################
make_intersect<-function(dirs,sep,row,col){
  # dirs<-c('D:\\cancer_type_classification\\tcga_data\\BLCA\\miRNAMatrix2.csv','D:\\cancer_type_classification\\tcga_data\\KIRC\\miRNAMatrix2.csv')
  # sep<-','
  # row<-0
  # col<-1
  
  li<-list()
  if (row==0){
    for (dir in dirs){
      tmp<-read.table(dir,sep=sep,check.names = F,header = T)[,col]
      li<-append(li,list(tmp))
    }
  }else{
    for (dir in dirs){
      tmp<-read.table(dir,sep=sep,check.names = F,header = T)[row,]
      li<-append(li,list(tmp))
    }
  }
  return(Reduce(intersect, li))
}
##################################################################################################
normalization_deseq2<-function(dir,sep=','){
  # dir<-'D:\\cancer_type_classification\\tcga_data\\KIRC\\mRNAMatrix2.csv'
  # sep<-','
  library("DESeq2")
  data<-read.table(dir,sep = sep,check.names = F,header = T)
  row.names(data)<-data$id
  data<-data[,-1]
  data<-as.matrix(data)
  meta<-cbind(colnames(data),c(rep(1,71),rep(0,509)))
  dds <- DESeqDataSetFromMatrix(countData = data,colData = meta,design = ~V2)
  dds <- estimateSizeFactors(dds)
  table_counts_normalized <- counts(dds, normalized=TRUE)
}
##################################################################################################
cpm<-function(dir,sep=','){
  # dir<-"D:\\cancer_type_classification\\tcga_data\\PAAD\\miRNAMatrix2.csv"
  # sep<-','
  
  data<-read.table(dir,sep = sep,header = T,check.names = F)
  row.names(data)<-data[,1]
  data<-data[,-1]
  res<-log(t(t(data)/colSums(data))*(10^6)+1,base=2)
  return(res)
}
##################################################################################################
log2x<-function(dir,sep=','){
  data<-read.table(dir,sep = sep,header = T,check.names = F)
  row.names(data)<-data[,1]
  data<-data[,-1]
  res<-log(data+1,base=2)
  return(res)
}
##################################################################################################
tmm<-function(dir,sep=',', op="pass"){
  # dir<-"D:\\lymph_single\\tcga_data\\STAD\\miRNAMatrix2.csv"
  # sep<-','
  
  library(psych)
  library(corpcor)
  library(parallel)
  library(stringr)
  library(data.table)
  library("edgeR")
  
  data<-read.table(dir,sep = sep,header = T,check.names = F)
  row.names(data)<-data[,1]
  data<-data[,-1]
  
  coln<-colnames(data)
  n1<-length(coln[which(str_sub(coln,start = 14, end = 15)=='11')])
  t1<-length(coln[which(str_sub(coln,start = 14, end = 15)=='01')])
  
  
  exp=as.matrix(data)
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  data=avereps(data)
  data=data[rowMeans(data)>1,]
  
  if (op=="pass"){
    return(data)
  }else{
    group=c(rep("normal",n1),rep("tumor",t1))                         
    design <- model.matrix(~group)
    y <- DGEList(counts=data,group=group)
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y,pair = c("normal","tumor"))
    topTags(et)
    ordered_tags <- topTags(et, n=100000)
    
    allDiff=ordered_tags$table
    allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
    diff=allDiff
    # newData=y$pseudo.counts
    newData <- edgeR::cpm(y, log=F)
    
    normalizeExp=newData
    normalizeExp=as.data.frame(normalizeExp)
    
    return(normalizeExp)
  }
}
##################################################################################################
integrate<-function(dirs,sep=',',bind="c"){
  library(utils)
  # dirs<-c('D:\\cancer_type_classification\\tcga_data\\KIRC\\normalized_g.csv', 'D:\\cancer_type_classification\\tcga_data\\KIRP\\normalized_g.csv')
  # sep<-','
  # bind<-"c"
  # library(progress)
  
  li<-list()
  
  cat("\nloading files\n")
  pb = txtProgressBar(min=0, max=length(dirs), style = 1, char="=")
  for (i in 1:length(dirs)){
    data<-read.table(dirs[i],sep = sep,header = T,check.names = F, row.names = 1)
    li<-append(li,list(data))
    setTxtProgressBar(pb, value = i)
  }
  close(pb)
  
  # res<-li[[1]]
  # res<-res[,FALSE]
  cat("\nintegrating matrices\n")
  # pb = txtProgressBar(min=0, max=length(dirs), style = 1, char="=")
  # for (i in 1:length(li)){
  #   setTxtProgressBar(pb, value = i)
  #   if (bind=="c"){
  #     res<-cbind(res,li[[i]])
  #   }else{
  #     res<-rbind(res,li[[i]])
  #   }
  # }
  # close(pb)

  if (bind=="c"){
    res<-do.call(cbind,li)
  }else{
    res<-do.call(rbind,li)
  }
  rm(li)
  return(res)
}
##################################################################################################
mircode<-function(dbdir,gene,mirna){
  # library(doParallel)
  # library(data.table)
  # library(limma)
  # library(stringr)
  # library(stringi)
  # library(readxl)
  # gene<-row.names(read.csv('D:\\cancer_type_classification\\tcga_data\\meta_ceRNA_g.csv',row.names=1,check.names=F))
  # mirna<-row.names(read.csv('D:\\cancer_type_classification\\tcga_data\\meta_mi.csv',row.names=1,check.names=F))
  # dbdir='D:\\cancer_type_classification\\target\\mircode_highconsfamilies.txt'
  
  db=fread(dbdir,header = T) 
  db=data.frame(db,check.names = F)
  db_new<-data.table() # blank data.table
  n <- nrow(db) # total number of rows
  
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(stringr))
  
  clusterExport(cl,c("db_new","n","db"),envir=environment())
  mydata1 <- parSapply(#delete delimiter
    
    cl,      
    1:n,    #whole number of combinations 
    function(i) {
      
      target_index <- as.character(db[i, 2]) 
      id_index <- as.character(db[i, 4]) 
      
      id_index_split_temp <- data.frame(strsplit(id_index, split = '/'))
      mart_temp <- data.frame(cbind(id_index_split_temp, target_index))
      
      names(mart_temp) <- c("id", "Target gene")
      
      db_new <- rbind(db_new, mart_temp)
      
      return(db_new)
    }  
  )
  stopCluster(cl)
  mamydata1=as.matrix(mydata1)
  id=unlist(mamydata1[1,])
  Targetgene=unlist(mamydata1[2,])
  db_new<-data.frame(id,Targetgene)
  rm(id,Targetgene,mydata1) # delete temp dataset
  cat("delete delimiter.\n")
  
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(stringr))
  
  clusterExport(cl,c("db_new"),envir=environment())
  mydata1 <- parSapply(#add hsa prefix
    
    cl,      
    1:nrow(db_new),    #whole number of combinations 
    function(i) {
      if (str_detect(db_new[i,1],'^[m][i][R]')==FALSE) {
        if (str_detect(db_new[i,1],'^[l][e][t]')==TRUE) {
          str_c("hsa-",db_new[i,1])
        }else{
          str_c("hsa-miR-",db_new[i,1])
        }
      }else{
        str_c("hsa-",db_new[i,1])
      }
      
    }  
  )
  stopCluster(cl)
  mamydata1=as.matrix(mydata1)
  mamydata1=as.data.frame(mamydata1)
  db_new[,1]<-mamydata1
  db_new<-db_new[!(db_new$id == "NULL"), ] 
  cat("attaching hsa- prefix.\n")
  # numCores <- parallel::detectCores() - 1
  
  db<-db_new
  db<-db[db[,2]%in%gene,]
  rm(db_new)
  
  
  id<-c()
  Targetgene<-c()
  cat("selecting valid interactions.\n")
  
  st<-Sys.time()
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(stringr))
  
  clusterExport(cl,c("db","id","Targetgene","mirna"),envir=environment())
  mydata1 <- parSapply(#add hsa prefix
    
    cl,      
    1:nrow(db),    #whole number of combinations 
    function(i) {
      k1<-str_split(db[i,1],'-',simplify = T)
      key1<-k1[3]
      digit1<-str_extract(key1,'\\d+')
      alpha1<-str_extract(key1,'[a-z]+')
      for (j in 1:length(mirna)){
        k2<-str_split(mirna[j],'-',simplify = T)
        key2<-k2[3]
        digit2<-str_extract(key2,'\\d+')
        alpha2<-str_extract(key2,'[a-z]')
        if (length(k1)==4 & length(k2)==4){
          if (is.na(alpha1)==FALSE & is.na(alpha2)==FALSE & k1[2]==k2[2] & k1[4]==k2[4]){
            if (digit1==digit2 & str_detect(alpha1,alpha2)==TRUE){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }else if (is.na(alpha1)==TRUE & is.na(alpha2)==TRUE & k1[2]==k2[2] & k1[4]==k2[4]){
            if (digit1==digit2){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }
        }else if (length(k1)==5 & length(k2)==5){
          if (is.na(alpha1)==FALSE & is.na(alpha2)==FALSE & k1[2]==k2[2] & k1[4]==k2[4] & k1[5]==k2[5]){
            if (digit1==digit2 & str_detect(alpha1,alpha2)==TRUE){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }else if (is.na(alpha1)==TRUE & is.na(alpha2)==TRUE & k1[2]==k2[2] & k1[4]==k2[4] & k1[5]==k2[5]){
            if (digit1==digit2){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }
        }else if (length(k1)==3 & length(k2)==3){
          if (is.na(alpha1)==FALSE & is.na(alpha2)==FALSE & k1[2]==k2[2]){
            if (digit1==digit2 & str_detect(alpha1,alpha2)==TRUE){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }else if (is.na(alpha1)==TRUE & is.na(alpha2)==TRUE & k1[2]==k2[2]){
            if (digit1==digit2){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }
        }
      }
      return(data.frame(id,Targetgene))
    }  
  )
  stopCluster(cl)
  et<-Sys.time()
  cat(et-st)
  cat("\n")
  db<-do.call(rbind,mydata1)
  
  db<-unique(db)
  
  db<-db[db[,1]%in%mirna==TRUE,]
  cat("sorted interactions in mircode.\n")
  return(db)
} 
##################################################################################################
targetscan<-function(dbdir,gene,mirna){
  # library(doParallel)
  # library(data.table)
  # library(limma)
  # library(stringr)
  # library(stringi)
  # library(readxl)
  # gene<-row.names(read.csv('D:\\cancer_type_classification\\tcga_data\\meta_ceRNA_g.csv',row.names=1,check.names=F))
  # mirna<-row.names(read.csv('D:\\cancer_type_classification\\tcga_data\\meta_mi.csv',row.names=1,check.names=F))
  # dbdir='D:\\cancer_type_classification\\target\\Conserved_Family_Info.txt'
  
  db=fread(dbdir,header = T) 
  db=data.frame(db,check.names = F)
  db<-db[db$`Species ID`==9606,]
  db_new<-data.table() # blank data.table
  n <- nrow(db) # total number of rows
  
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(stringr))
  
  clusterExport(cl,c("db_new","n","db"),envir=environment())
  mydata1 <- parSapply(#delete delimiter
    
    cl,      
    1:n,    #whole number of combinations 
    function(i) {
      
      target_index <- as.character(db[i, 3]) 
      id_index <- as.character(db[i, 1]) 
      
      id_index_split_temp <- data.frame(strsplit(id_index, split = '/'))
      mart_temp <- data.frame(cbind(id_index_split_temp, target_index))
      
      names(mart_temp) <- c("id", "Target gene")
      
      db_new <- rbind(db_new, mart_temp)
      
      return(db_new)
    }  
  )
  stopCluster(cl)
  mamydata1=as.matrix(mydata1)
  id=unlist(mamydata1[1,])
  Targetgene=unlist(mamydata1[2,])
  db_new<-data.frame(id,Targetgene)
  rm(id,Targetgene,mydata1) # delete temp dataset
  cat("delete delimiter.\n")
  
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(stringr))
  
  clusterExport(cl,c("db_new"),envir=environment())
  mydata1 <- parSapply(#add hsa prefix
    
    cl,      
    1:nrow(db_new),    #whole number of combinations 
    function(i) {
      if (str_detect(db_new[i,1],'^[m][i][R]')==FALSE) {
        if (str_detect(db_new[i,1],'^[l][e][t]')==TRUE) {
          str_c("hsa-",db_new[i,1])
        }else{
          str_c("hsa-miR-",db_new[i,1])
        }
      }else{
        str_c("hsa-",db_new[i,1])
      }
      
    }  
  )
  stopCluster(cl)
  mamydata1=as.matrix(mydata1)
  mamydata1=as.data.frame(mamydata1)
  db_new[,1]<-mamydata1
  db_new<-db_new[!(db_new$id == "NULL"), ] 
  cat("attaching hsa- prefix.\n")
  # numCores <- parallel::detectCores() - 1
  
  db<-db_new
  db<-db[db[,2]%in%gene,]
  rm(db_new)
  id<-c()
  Targetgene<-c()
  cat("selecting valid interactions.\n")

  st<-Sys.time()
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(stringr))
  
  clusterExport(cl,c("db","id","Targetgene","mirna"),envir=environment())
  mydata1 <- parSapply(#add hsa prefix
    
    cl,      
    1:nrow(db),    #whole number of combinations 
    function(i) {
      k1<-str_split(db[i,1],'-',simplify = T)
      key1<-k1[3]
      digit1<-str_extract(key1,'\\d+')
      alpha1<-str_extract(key1,'[a-z]+')
      for (j in 1:length(mirna)){
        k2<-str_split(mirna[j],'-',simplify = T)
        key2<-k2[3]
        digit2<-str_extract(key2,'\\d+')
        alpha2<-str_extract(key2,'[a-z]')
        if (length(k1)==4 & length(k2)==4){
          if (is.na(alpha1)==FALSE & is.na(alpha2)==FALSE & k1[2]==k2[2] & k1[4]==k2[4]){
            if (digit1==digit2 & str_detect(alpha1,alpha2)==TRUE){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }else if (is.na(alpha1)==TRUE & is.na(alpha2)==TRUE & k1[2]==k2[2] & k1[4]==k2[4]){
            if (digit1==digit2){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }
        }else if (length(k1)==5 & length(k2)==5){
          if (is.na(alpha1)==FALSE & is.na(alpha2)==FALSE & k1[2]==k2[2] & k1[4]==k2[4] & k1[5]==k2[5]){
            if (digit1==digit2 & str_detect(alpha1,alpha2)==TRUE){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }else if (is.na(alpha1)==TRUE & is.na(alpha2)==TRUE & k1[2]==k2[2] & k1[4]==k2[4] & k1[5]==k2[5]){
            if (digit1==digit2){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }
        }else if (length(k1)==3 & length(k2)==3){
          if (is.na(alpha1)==FALSE & is.na(alpha2)==FALSE & k1[2]==k2[2]){
            if (digit1==digit2 & str_detect(alpha1,alpha2)==TRUE){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }else if (is.na(alpha1)==TRUE & is.na(alpha2)==TRUE & k1[2]==k2[2]){
            if (digit1==digit2){
              id<-append(id,mirna[j])
              Targetgene<-append(Targetgene,db[i,2])
            }
          }
        }
      }
      return(data.frame(id,Targetgene))
    }  
  )
  stopCluster(cl)
  et<-Sys.time()
  cat(et-st)
  cat("\n")
  db<-do.call(rbind,mydata1)
  
  db<-unique(db)
  
  db<-db[db[,1]%in%mirna==TRUE,]
  cat("sorted interactions in targetscan.\n")
  return(db)
}
##################################################################################################
mirtarbase<-function(dbdir,gene,mirna){
  # library(doParallel)
  # library(data.table)
  # library(limma)
  # library(stringr)
  # library(stringi)
  library(readxl)
  # gene<-row.names(read.csv('D:\\cancer_type_classification\\tcga_data\\meta_ceRNA_g.csv',row.names=1,check.names=F))
  # mirna<-row.names(read.csv('D:\\cancer_type_classification\\tcga_data\\meta_mi.csv',row.names=1,check.names=F))
  # dbdir='D:\\cancer_type_classification\\target\\hsa_MTI.xlsx'
  
  db=read_xlsx(dbdir)
  db=data.frame(db,check.names = F)
  db<-cbind(db[,2],db[,4])
  colnames(db)<-c('id','Targetgene')

  db<-db[db[,2]%in%gene,]
  
  db<-unique(db)
  
  db<-db[db[,1]%in%mirna,]
  
  cat("sorted interactions in mirtarbase.\n")
  return(db)
}
##################################################################################################
target_post<-function(pair_tab){
  
  # pair_tab<-as.data.frame(pair_tab)
  # pair_tab<-unique(pair_tab)
  b<-pair_tab
  
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(stringr))
  clusterEvalQ(cl,library(data.table))
  
  clusterExport(cl,c("pair_tab","b"),envir=environment())
  # 
  mydata1 <- parSapply(#捍纺贸府
    
    cl,      
    1:nrow(pair_tab),    #whole number of combinations 
    function(i) {
      res<-c()
      for (j in 1:nrow(b)){
        if (pair_tab[i,1]==b[j,1]){
          res<-append(res,str_c(pair_tab[i,1],"//",pair_tab[i,2],"//",b[j,2]))
        }
      }
      return(res)
    }  
  )
  stopCluster(cl)
  mamydata1=as.matrix(mydata1)
  mamydata1=as.data.frame(mamydata1)
  mamydata1<-unlist(mamydata1)
  mamydata1<-as.data.frame(mamydata1)
  test<-mamydata1
  DUPLICATED_fix = test[!duplicated(test$mamydata1),]
  name_strsplit<-data.frame(do.call('rbind', strsplit(as.character(test$mamydata1), split = '//', fixed = TRUE)))
  name_strsplit<-name_strsplit[!(name_strsplit[,2]==name_strsplit[,3]),]
  name_strsplit<-unique(name_strsplit)
  
  res<-unique(t(apply(name_strsplit, 1, sort)))
  return(res)
}
##################################################################################################
geneid_converter<-function(vec,query='en to hugo',setdir){
  ref<-read.csv(str_c(setdir,"\\ensembl\\annotation.csv"),header = T,check.names = F)
  # vec<-c("ENSG00000000003","ENSG00000000419")
  
  if (query=="en to hugo"){
    li<-list()
    for (i in 1:length(vec)){
      tmp<-ref[ref[,1]==vec[i],]
      li<-append(li,list(tmp))
    }
    res<-do.call(rbind,li)
    return(res[,2:3])
  }else{
    li<-list()
    for (i in 1:length(vec)){
      tmp<-ref[ref[,2]==vec[i],]
      li<-append(li,list(tmp))
    }
    res<-do.call(rbind,li)
    return(cbind(res[,1],res[,3]))
  }
}
##################################################################################################
spearman<-function(mat,filter){#calculate scc
  # mat<-read.csv('D:\\cancer_type_classification\\tcga_data\\inte_mi.csv',row.names=1)[1:120,]
  # filter<-0.6
  library(parallel)
  library(utils)
  library(Hmisc)
  library(igraph)
  
  s_t<-Sys.time()
  scc<-cor(scale(as.matrix(t(mat))),method="pearson")
  
  scc[is.na(scc)] <- 0
  
  cat("\ntransforming scc matrix to adjacency matrix\n")
  
  cl = makeCluster(parallel::detectCores() - 9)
  clusterEvalQ(cl,library(ggm))
  clusterEvalQ(cl,library(corpcor))
  
  clusterExport(cl,c("scc","filter"),envir=environment())
  
  mydata1 <- parSapply(#捍纺贸府
    
    cl,      
    1:nrow(scc),    #whole number of combinations 
    function(i) {
      rowvec<-ifelse(scc[i,]>filter,1,0)
      return(rowvec)
    }  
  )
  stopCluster(cl)
  adj=data.frame(t(mydata1),check.names = F)
  rm(mydata1,scc)
  ###
  # adj<-adj[-1,]
  row.names(adj)<-colnames(adj)
  
  adj<-adj[rowSums(adj)>1,]
  adj<-adj[,colSums(adj)>1]
  
  g<-graph.adjacency(as.matrix(adj))
  rm(adj)
  res<-get.edgelist(g)
  colnames(res)<-c("source","target")
  res<-as.data.frame(res)
  res<-res[!(res$source==res$target),]
  
  e_t<-Sys.time()
  print(s_t-e_t)
  
  return(res)
}
##################################################################################################
spearman_ceRNA<-function(interaction_tab,ex1,ex2,filter){
  # interaction_tab<-read.csv('D:\\cancer_type_classification\\target\\mirtarbase_meta.csv')
  # ex_dir1<-'D:\\cancer_type_classification\\tcga_data\\meta_mi.csv'
  # ex_dir2<-'D:\\cancer_type_classification\\tcga_data\\meta_ceRNA_g.csv'
  # filter<-0.5
  
  # ex1<-read.csv(ex_dir1,row.names = 1,check.names = F)
  # ex2<-read.csv(ex_dir2,row.names = 1,check.names = F)
  cat('calculating scc\n')
  
  if (ncol(interaction_tab)==2){
    # st<-Sys.time()
    cl = makeCluster(parallel::detectCores() - 1)
    clusterEvalQ(cl,library(ggm))
    clusterEvalQ(cl,library(corpcor))
    
    clusterExport(cl,c("ex1","ex2","interaction_tab"),envir=environment())
    
    scc <- parSapply(#捍纺贸府
      cl,      
      1:nrow(interaction_tab),    #whole number of combinations 
      function(i) {
        # cox_all=matrix(nrow = 1, ncol = 1)
        # i=1
        
        # ce1_1<-as.character(interaction_tab[i,1])
        # ce2_1<-as.character(interaction_tab[i,2])
        # s1<-cbind(t(ex1[ce1_1,]), t(ex2[ce2_1,]))
        xcor=cor(t(ex1[interaction_tab[i,1],]),t(ex2[interaction_tab[i,2],]), method = "pearson")
        return(xcor)
      }  
    )
    stopCluster(cl)
    # et<-Sys.time()
    # print(et-st)
    res<-cbind(interaction_tab,scc)
    res<-res[abs(res[,3])>filter,]
    return(res)
  }else if (ncol(interaction_tab)==3){
    cl = makeCluster(parallel::detectCores() - 1)
    clusterEvalQ(cl,library(ggm))
    clusterEvalQ(cl,library(corpcor))
    
    clusterExport(cl,c("ex1","ex2","interaction_tab"),envir=environment())
    
    mydata1 <- parSapply(#捍纺贸府
      cl,      
      1:nrow(interaction_tab),    #whole number of combinations 
      function(i) {
        cox_all=matrix(nrow = 3, ncol = 1)
        
        ce1_1= as.character(interaction_tab[i,1])
        ce2_1= as.character(interaction_tab[i,2])
        miRNA1= as.character(interaction_tab[i,3])
        
        s1<-cbind(t(ex2[ce1_1,]), t(ex2[ce2_1,]), t(ex1[miRNA1,]))
        xcor=cor(s1,method = "pearson")
        cox_all[1,1]=xcor[2,1]
        cox_all[2,1]=xcor[3,1]
        cox_all[3,1]=xcor[3,2]
        return(cox_all)
      }  
    )
    stopCluster(cl)
    scc<-data.frame(mydata1)
    scc<-t(scc)
    res<-cbind(interaction_tab,scc)
    colnames(res)<-c('x','y','miRNA','x_y','mi_x','mi_y')
    
    #post process of scc
    res<-res[res$x_y>filter,]#select triplets with |pcc|>filter
    res<-res[abs(res$mi_x)>filter & abs(res$mi_y)>filter & (res$mi_y)*(res$mi_x)>0,]
    # res<-res[abs(res$mi_x)>filter,]
    # res<-res[abs(res$mi_y)>filter,]
    return(res)
  }
}
##################################################################################################
delta_pcc_ceRNA<-function(pair_tab,normal_g,normal_mi,target_g,target_mi){
  # i<-"BRCA"
  # pair_tab<-read.csv(str_c(setdir,"\\tcga_data\\",i,"\\admat.csv"),check.names=F)
  # target_g<-read.csv(str_c(setdir,"\\tcga_data\\",i,"\\meta_normalized_g.csv"),row.names = 1, check.names = F)
  # vec<-row.names(target_g)
  # id_tab<-geneid_converter(vec,setdir = setdir)
  # row.names(target_g)<-id_tab[,1]
  # target_mi<-read.csv(str_c(setdir,"\\tcga_data\\",i,"\\meta_normalized_mi.csv"),row.names = 1, check.names = F)
  # normal_g<-read.csv(str_c(setdir,"\\tcga_data\\",i,"\\normal_normalized_g.csv"),row.names = 1, check.names = F)
  # normal_mi<-read.csv(str_c(setdir,"\\tcga_data\\",i,"\\normal_normalized_mi.csv"),row.names = 1, check.names = F)
  
  nodes<-unique(c(pair_tab[,1],pair_tab[,2]))
  
  normal_g<-normal_g[row.names(normal_g)%in%nodes,]
  normal_mi<-normal_mi[row.names(normal_mi)%in%nodes,]
  target_g<-target_g[row.names(target_g)%in%nodes,]
  target_mi<-target_mi[row.names(target_mi)%in%nodes,]
  
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(ggm))
  clusterEvalQ(cl,library(corpcor))
  
  clusterExport(cl,c("normal_g","normal_mi","target_g","target_mi","pair_tab"),envir=environment())
  
  scc <- parSapply(#捍纺贸府
    cl,      
    1:nrow(pair_tab),    #whole number of combinations 
    function(i) {
      delta<-numeric(ncol(target_g))
      
      exp_g<-c(t(normal_g[pair_tab[i,2],]))
      exp_mi<-c(t(normal_mi[pair_tab[i,1],]))
      cor=cor(exp_g, exp_mi, method = "pearson")
      
      for (j in 1:ncol(target_g)){
        # i<-1
        # j<-2
        
        exp_g<-c(t(normal_g[pair_tab[i,2],]),target_g[pair_tab[i,2],j])
        exp_mi<-c(t(normal_mi[pair_tab[i,1],]),target_mi[pair_tab[i,1],j])
        xcor=cor(exp_g, exp_mi, method = "pearson")
        
        delta[j]<-(xcor-cor)
      }
      return(delta)
    }  
  )
  stopCluster(cl)
  res<-cbind(pair_tab,t(scc))
  colnames(res)<-c("mi","rna",colnames(target_g))
  # res<-na.omit(res)
  
  return(res)
}
##################################################################################################
delta_wilcox_test<-function(group1,group2,p_filter){
  # group1<-meta_delta
  # group2<-nonmeta_delta
  # p_filter<-0.00001
  
  pair_tab<-group1[,1:2]
  group1<-group1[,3:ncol(group1)]
  group2<-group2[,3:ncol(group2)]
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(ggm))
  clusterEvalQ(cl,library(corpcor))
  
  clusterExport(cl,c("group1","group2"),envir=environment())
  
  p <- parSapply(#捍纺贸府
    cl,      
    1:nrow(pair_tab),    #whole number of combinations 
    function(i) {
      # i<-1
      # res<-try({res<-wilcox.test(t(group1[i,]),t(group2[i,]))},silent = T)
      
      Pvalue=tryCatch({
        Pvalue=wilcox.test(t(group1[i,]),t(group2[i,]))$p.value
      }, warning = function(w){
        return(2)
      }, 
      error = function(e) {
        return(3)
      }
      )
      return(Pvalue)
    }  
  )
  stopCluster(cl)
  res<-cbind(pair_tab,p)
  res<-res[res[,3]<p_filter,]
  return(res)
}
##################################################################################################
pspearman<-function(mat,filter){
  # mat<-read.csv('D:\\cancer_type_classification\\tcga_data\\inte_g.csv',row.names=1)[1:120,]
  # filter<-0.6
  library(parallel)
  st<-Sys.time()
  x<-row.names(mat)
  com<-t(combn(x,2))
  
  cl = makeCluster(parallel::detectCores() - 1)
  clusterEvalQ(cl,library(ggm))
  clusterEvalQ(cl,library(corpcor))
  
  clusterExport(cl,c("mat","com"),envir=environment())
  
  mydata1 <- parSapply(#捍纺贸府
    
    cl,      
    1:nrow(com),    #whole number of combinations 
    function(i) {
      
      gene1= as.character(com[i,1])
      gene2= as.character(com[i,2])
      
      s1=cbind(t(mat[gene1,]),t(mat[gene2,]))
      xcor=cor(s1,method = "spearman")[1,2]
      return(xcor)
    }  
  )
  stopCluster(cl)
  res=data.frame(mydata1)
  et<-Sys.time()
  print(et-st)
  res<-cbind(com,res)
  colnames(res)<-c('gene1','gene2','scc')
  return(res)
}
##################################################################################################
minmax<-function(mat){
  # mat<-read.csv('D:\\cancer_type_classification\\tcga_data\\BLCA\\normalized_mi.csv',row.names = 1)
  # mat<-gctm[,1]
  
  library(caret)
  
  if(is.null(ncol(mat))==T){
    norm_scale<-(mat-min(mat))/(max(mat)-min(mat))
  }else{
    process <- preProcess(mat, method=c("range"))
    norm_scale <- predict(process, mat)
  }
  return(norm_scale)
}
##################################################################################################
net_info<-function(edge_tab){
  cat('\nnum of edges: \n')
  cat(nrow(edge_tab))
  cat('\nnum of nodes: \n')
  cat(length(unique(c(edge_tab[,1],edge_tab[,2]))))
}
##################################################################################################
round_csv<-function(dir,r){
  # dir<-'D:\\cancer_type_classification\\tcga_data\\inte_mi.csv'
    
  dt<-read.csv(dir,row.names = 1,check.names = F)
  res<-round(dt,r)
  write.csv(res,dir)
  rm(dt,res)
}
##################################################################################################
integerize<-function(vec){
  # vec<-row.names(dt)
  
  integerized<-1:length(vec)
  mapping<-cbind(vec,integerized)
  return(mapping)
}
##################################################################################################
transition<-function(row,mapping){
  #row<-admat[1,]
  mapping<-id_mapping
  row[,1]<-id_mapping[id_mapping[,1]==row[,1],2]
  row[,2]<-id_mapping[id_mapping[,1]==row[,2],2]
  return(row)
}
##################################################################################################
remove_normal<-function(dir){
  # dir<-"D:\\cancer_type_classification\\tcga_data\\KIRC\\miRNAMatrix2.csv"
  
  dt<-read.csv(dir,header = T,row.names = 1,check.names = F)
  dt=dt[,which(str_sub(colnames(dt),14,14)!="1")]
  write.csv(dt,dir,quote = F)
}
##################################################################################################
separate_normal<-function(dir){
  # spl<-str_split(dir,'\\\\')
  # outputname<-tail(spl[[1]],1)
  
  dt<-read.csv(dir,header = T,row.names = 1,check.names = F)
  dt=dt[,which(str_sub(colnames(dt),14,14)=="1")]
  # write.csv(dt,str_c(dir,"\\normal.csv"),quote = F)
  return(dt)
}
##################################################################################################
pick_solid<-function(dir){
  # dir<-"D:\\cancer_type_classification\\tcga_data\\BRCA\\miRNAMatrix2.csv"
  
  dt<-read.csv(dir,header = T,row.names = 1,check.names = F)
  dt=dt[,which(str_sub(colnames(dt),14,15)=="11" | str_sub(colnames(dt),14,15)=="01")]
  write.csv(dt,dir,quote = F)
}
##################################################################################################
colname_slice<-function(dir,st,ed){
  # dir<-""
  # st<-1
  # ed<-12
  
  dt<-read.csv(dir,header = T,row.names = 1,check.names = F)
  cn<-colnames(dt)
  for (i in 1:length(cn)){
    cn[i]<-str_sub(cn[i],start = st,end = ed)
  }
  colnames(dt)<-cn
  write.csv(dt,dir,quote = F)
}
##################################################################################################
colname_add<-function(dir,str){
  dt<-read.csv(dir,header = T,row.names = 1,check.names = F)
  cn<-colnames(dt)
  for (i in 1:length(cn)){
    cn[i]<-str_c(cn[i],str)
  }
  colnames(dt)<-cn
  write.csv(dt,dir,quote = F)
}
##################################################################################################
lymp_meta_separator<-function(dir,clinic_dir,setwd, op){
  # dir<-'D:\\cancer_type_classification\\tcga_data\\KIRC\\normalized_mi.csv'
  # clinic_dir<-'D:\\cancer_type_classification\\tcga_data\\KIRC\\TCGA-KIRC-clinical.csv'
  
  setwd(setwd)
  
  spl<-str_split(dir,'\\\\')
  outputname<-tail(spl[[1]],1)
  
  dt<-read.csv(dir,row.names = 1,check.names = F)
  
  clinic<-read.csv(clinic_dir,row.names = 1,check.names = F)
  clinic<-clinic[which(str_sub(clinic$ajcc_pathologic_m,1,2)!='MX' &
                       str_sub(clinic$ajcc_pathologic_t,1,2)!='TX' &
                       str_sub(clinic$ajcc_pathologic_n,1,2)!='NX'),]
  
  if (op=="bol"){
    #metastasis
    meta=clinic[which(str_sub(clinic$ajcc_pathologic_n,1,2)!='N0'),] 
    meta=meta[which(str_sub(meta$ajcc_pathologic_t,1,2)!='T0'),] 
    meta=meta[which(str_sub(meta$ajcc_pathologic_m,1,2)!='M1'),]
    #non-n-metastasis
    nonmeta=clinic[which(str_sub(clinic$ajcc_pathologic_n,1,2)=='N0'),] 
    nonmeta=nonmeta[which(str_sub(nonmeta$ajcc_pathologic_t,1,2)!='T0'),] 
    nonmeta=nonmeta[which(str_sub(nonmeta$ajcc_pathologic_m,1,2)!='M1'),] 
  }else if (op=="men"){
    #metastasis
    meta=clinic[which(str_sub(clinic$ajcc_pathologic_n,1,2)!='N0'),] 
    meta=meta[which(str_sub(meta$ajcc_pathologic_t,1,2)!='T0'),] 
    #non-n-metastasis
    nonmeta=clinic[which(str_sub(clinic$ajcc_pathologic_n,1,2)=='N0'),] 
    nonmeta=nonmeta[which(str_sub(nonmeta$ajcc_pathologic_t,1,2)!='T0'),] 
  }
  
  dt_meta<-dt[,colnames(dt)%in%(meta$submitter_id)]
  dt_nonmeta<-dt[,colnames(dt)%in%(nonmeta$submitter_id)]
  
  fwrite(dt_meta,str_c('meta_',outputname),sep = ',',row.names = T,quote = F)
  fwrite(dt_nonmeta,str_c('nonmeta_',outputname),sep = ',',row.names = T,quote = F)
}
##################################################################################################
tocyto<-function(triplet){
  triplet<-triplet[,1:3]
  tmp1<-triplet[,2:3]
  tmp2<-cbind(triplet[,1],triplet[,3])
  colnames(tmp1)<-c('g','mi')
  colnames(tmp2)<-c('g','mi')
  res<-unique(rbind(tmp1,tmp2))
  return(res)
}
##################################################################################################
cyto_type<-function(vec,setdir){
  # vec<-c("CDC20","PCLAF","hsa-miR-18a-5p")
  
  ref<-read.csv(str_c(setdir,"\\ensembl\\Homo_sapiens1.csv"),check.names = F, header = F)
  ref<-unique(ref[,2:3])
  type<-c()
  for(i in 1:length(vec)){
    if(vec[i]%in%ref[,1]){
      type<-append(type,ref[ref[,1]==vec[i],2])
    }else{
      type<-append(type,"miRNA")
    }
  }
  res<-cbind(vec,type)
  return(res)
}
##################################################################################################