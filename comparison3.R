accepted <- c("0|1", "1|0", "1|2", "2|1")
sswap <- function(str){
  lst <- strsplit(str,"\\|")
  s1 <- lst[[1]][1]
  s2 <- lst[[1]][2]
  return (paste(s2,"|",s1, sep=""))
}
simcode <- function(str){
  lst <- strsplit(str,"_")
  s1 <- lst[[1]][1]
  s2 <- lst[[1]][2]
  if(s1==s2)
    return(1)
  else if (sswap(s1)==s2)
    return(2)
  #otherwise
  return(3)
}
findblocks <- function(v, el){
  #in vector v finds index  of blocks of el
  idx <-numeric(length=0)
  block_start <- -1
  for(i in 1:length(v)){
    if (v[i]==el){
      if (block_start == -1)
         block_start <- i
    }else{
      if (block_start != -1){
        idx <- c(idx,block_start)
        block_start <- -1
      }
    }
  }
  if(block_start != -1)
    idx <- c(idx,block_start)
  
  return (idx)
  
}
switcherror <- function(v){
  vv <- sapply(v, FUN=simcode)  
   #count number of blocks of 1,2 and 3
  #example if we have 11122111 then assuming phase switch error will be 1, assuming opposite-phase then switch error will be 2, so the final switch error is min(1,2)=1
  #assuming phase with reference then we count blocks of 2 and 3 (none 1):
  blocks1 <- findblocks(vv,1)
  blocks2 <- findblocks(vv,2)
  blocks3 <- findblocks(vv,3)
  phased <- c(blocks2 , blocks3)
  #assuming unphase with reference then we count blocks of 1 and 3 (none 2):
  
  unphased <- c(blocks1 , blocks3)
  
  return (min(length(phased), length(unphased)))
}



n50 <- function(contig_lens){
  contig_lens <- sort(contig_lens, decreasing = T)
  n50_value <- 0.5 * sum(contig_lens);
  cum_lens <- cumsum(contig_lens)
  for(i in 1:length(contig_lens)){
    if (cum_lens[i]>= n50_value)
      return (contig_lens[i])
  }
  
  return (-1)
  
}

dfvcf<-read.table("hg19/NA12878_sorted_noheader.vcf", header=F, colClasses=c("character","numeric","NULL","NULL","NULL","NULL","NULL","NULL","NULL","character"))
#only keep het positiosn
dfvcf <- dfvcf[dfvcf$V10 %in% c('1|0','0|1','1|2','2|1','0|2','2|0'),]
dfvcf$snp_num <- c(1:nrow(dfvcf))
dim(dfvcf)
#2,490,062  4
dfhaplominer <-read.table("NA12878_hifi_haplominer.hap", colClasses=c("character", "numeric","character","character","character","character") ,header=F, fill=T)
dfhaplominer[dfhaplominer=='']<-NA
head(dfhaplominer)
dfhaplominer_rows <-dfhaplominer[!is.na(dfhaplominer$V6), ]
factors <- dfhaplominer[is.na(dfhaplominer$V6), 4]
dfhaplominer_rows$rownum <- c(1:nrow(dfhaplominer_rows))
dfmerged <-merge(dfhaplominer_rows,dfvcf, by=c(1,2), all.x=T)
dfmerged <- dfmerged[order(dfmerged$rownum),]
rowfactor <- rep(c(1:length(factors)), factors)
dfmerged$rowfactor <- rowfactor
dfmerged$both <- paste(dfmerged$V6,'_',dfmerged$V10,sep="")
agg<-aggregate(both~rowfactor, dfmerged, FUN=function(b){switcherror(b)})

tmp <- dfhaplominer[is.na(dfhaplominer$V6), 2:3]
haplominer_lens <- as.numeric(tmp$V3)-as.numeric(tmp$V2) + 1
n50(haplominer_lens)
#31,638
sum(haplominer_lens)
# 1,540,619,001
length(haplominer_lens)
#100,606
mean(haplominer_lens)
#15,313
max(haplominer_lens)
#315,905

###########################################
#now estimate haplominer swtich rate per kbp of long haplotype
length(which(haplominer_lens<1000))
#33,664
length(which(haplominer_lens>=1000))
#149,746
agg$len_bp <- haplominer_lens
agg_long <- agg[agg$len_bp>=2,]
(sum(agg_long$both)*1000)/sum(agg_long$len_bp)
#0.00267
##########################################################################
############################################################################
#hapcut2
dfhapcut <- read.table("NA12878_hifi_hapcut2_pruned.hap",sep = "\t" , fill=T)
dfhapcut[dfhapcut=='']<-NA
dfhapcut_rows <- dfhapcut[!is.na(dfhapcut$V2) ,c(4,5,2,3)]
dfhapcut_rows$V4 <- as.character(dfhapcut_rows$V4)
dfhapcut_rows$V5 <- as.numeric(dfhapcut_rows$V5)
summary_vec <-which(dfhapcut$V1=="******** ")
factors3<-  summary_vec - c(0, summary_vec[-length(summary_vec)])  - 2
factors3 <- c(factors3,2)
sum(factors3)
#1,978,496
nrow(dfhapcut_rows)
#1,978,496
dfhapcut_rows$V6 <- paste(dfhapcut_rows$V2,"|", dfhapcut_rows$V3, sep="")
dfhapcut_rows$rownum <- c(1:nrow(dfhapcut_rows))
dfmerged3 <-merge(dfhapcut_rows,dfvcf, by=c(1,2), all.x=T)
dfmerged3 <- dfmerged3[order(dfmerged3$rownum),]
dup_idx <- which(duplicated(dfmerged3$rownum))
rowfactor <- rep(c(1:length(factors3)), factors3)
rowfactor <- c(rowfactor, rowfactor[length(rowfactor)])
dfmerged3$rowfactor <- rowfactor
dfmerged3$both <- paste(dfmerged3$V6,'_',dfmerged3$V10,sep="")
agg3<-aggregate(both~rowfactor, dfmerged3, FUN=function(b){switcherror(b)})
#work out hapcut_lens
hapcut_lens <- numeric(length=0)
sidx <- 0
for(hl in factors3){
  tmp <- dfmerged3[sidx+hl,"V5"] - dfmerged3[sidx+1,"V5"] + 1
  hapcut_lens <- c(hapcut_lens, tmp)
  sidx <- sidx + hl
}
###########################################
#now estimate hapcut swtich rate per kbp of long haplotype
length(which(hapcut_lens<1000))
#2079
length(which(hapcut_lens>=1000))
#37,299
agg3$len_bp <- hapcut_lens
agg3_long <- agg3[agg3$len_bp>=2,]
(sum(agg3_long$both)*1000)/sum(agg3_long$len_bp)
#0.0010
#############################################
#because hapcut2 reports non-continuous SNP as haplotype: we write a code to detect only conitinuous ones:
v<-dfmerged3$snp_num
series3<-numeric(length=0)
hapcut_lens <- numeric(length=0)
start_idx <- 1
for(fact in factors3){
  
  split_vec <- numeric(length=0)
  for (i in 0:(fact-2)){
    if(v[start_idx+i+1]-v[start_idx+i] >1){
      split_vec <- c(split_vec,i)
    }
  }
  if(length(split_vec)==0){
    series3 <- c(series3, fact)
    hapcut_lens <- c(hapcut_lens, dfhapcut_rows[start_idx + fact-1,"V5"]- dfhapcut_rows[start_idx,"V5"] + 1)
  }else{
    last <- -1
    s_idx <- start_idx
    for(s in split_vec){
      series3 <- c(series3, s - last)
      hapcut_lens <- c(hapcut_lens, dfhapcut_rows[s_idx + s - last-1,"V5"]- dfhapcut_rows[s_idx,"V5"] + 1)
      s_idx <- s_idx + s - last
      last <- s
      
      
    }     
    series3 <- c(series3, fact - last -1)
    hapcut_lens <- c(hapcut_lens, dfhapcut_rows[s_idx + fact - last -2,"V5"]- dfhapcut_rows[s_idx,"V5"] + 1)
  }
  start_idx <- start_idx + fact
}
n50(hapcut_lens)
#8,288
sum(hapcut_lens)
#1,149,388,869
mean(hapcut_lens)
#3750.
max(hapcut_lens)
#92,696
max(series3)
#110
length(hapcut_lens)
#306,442

###########################################################################################################
###########################################################################################################
#WhatsHap

dfwhatshap<-read.table("NA12878_hifi_whatshap.vcf", header=F, colClasses=c("character","numeric","NULL","NULL","NULL","NULL","NULL","NULL","NULL","character"))
hap.list <- list()
haplotype_started <- F
last_pos <- 0 # last position
last_start_pos <- 0

haplotype_num <- 0
cnt = 1;
for(i in 1: nrow(dfwhatshap)){
  ss <- strsplit(dfwhatshap[i,"V10"], ":")
  if(length(ss[[1]])>1){
    #haplotype 
    start_pos <- as.numeric(ss[[1]][2])
    if(!haplotype_started){
      haplotype_started <- T
      #start_pos <- dfwhatshap[i,"V2"]#pos
      haplotype_num <- haplotype_num + 1
    }else if (last_start_pos!=0 && start_pos != last_start_pos ){
      #a new haplotype must have been started
      haplotype_num <- haplotype_num + 1
    }
    # now write haplotype line 
    #dfwhatshap_rows <- rbind(dfwhatshap_rows, list(chr=dfwhatshap[i,"V1"], pos=dfwhatshap[i,"V2"], GN=ss[[1]][1], rowfactor=haplotype_num))
    #dfwhatshap_rows[cnt,] <- list(dfwhatshap[i,"V1"],dfwhatshap[i,"V2"],ss[[1]][1],haplotype_num)
    hap.list[[cnt]] <- list(dfwhatshap[i,"V1"],dfwhatshap[i,"V2"],ss[[1]][1],haplotype_num)
    cnt <- cnt + 1
    last_start_pos <- start_pos
  }else{
    haplotype_started <- F
  }
  last_pos <- dfwhatshap[i,"V2"]
  
}
dfhwhathaps_rows <- do.call("rbind", hap.list) # this will return a matrix
dfhwhathaps_rows <- as.data.frame(dfhwhathaps_rows)
dfhwhathaps_rows$V1 <- as.character(dfhwhathaps_rows$V1)
dfhwhathaps_rows$V2 <- as.integer(dfhwhathaps_rows$V2)
dfhwhathaps_rows$V3 <- as.character(dfhwhathaps_rows$V3)
dfhwhathaps_rows$V4 <- as.integer(dfhwhathaps_rows$V4)
colnames(dfhwhathaps_rows)[4] <- "rowfactor"
dim(dfhwhathaps_rows)
#[1] 2,407,735      4
hap_lens<-aggregate(V1~rowfactor, dfhwhathaps_rows, FUN=length)
length(which(hap_lens$V1==1))
#[1] 145,407  haplotypes are with length 1
#now filter rows that contain hap with length 1
dfhwhathaps_rows <- dfhwhathaps_rows[!dfhwhathaps_rows$rowfactor %in% hap_lens[hap_lens$V1==1,"rowfactor"] , ]
dim(dfhwhathaps_rows)
# 2,262,328      4

dfhwhathaps_rows$rownum <- c(1:nrow(dfhwhathaps_rows))
dfmerged4 <-merge(dfhwhathaps_rows,dfvcf, by=c(1,2), all.x=T)
dfmerged4 <- dfmerged4[order(dfmerged4$rownum),]
dup_idx <- which(duplicated(dfmerged4$rownum))
dfagg4_len<-aggregate(both~rowfactor, dfmerged4, FUN=length)
# As WhatsHap reports non-continuous SNPs as haplotype: we require code 
# to detect only continuous ones

v<-dfmerged4$snp_num
series4<-numeric(length=0)#new haplotype lengths are stored here
whatshap_lens <- numeric(length=0)
start_idx <- 1

for(fact in agg4_len$both){
  
  split_vec <- numeric(length=0)
  for (i in 0:(fact-2)){
    if(v[start_idx+i+1]-v[start_idx+i] >1){
      split_vec <- c(split_vec,i)
    }
  }
  if(length(split_vec)==0){
    series4 <- c(series4, fact)
    whatshap_lens <- c(whatshap_lens, dfhwhathaps_rows[start_idx + fact-1,"V2"]- dfhwhathaps_rows[start_idx,"V2"] + 1)
  }else{
    last <- -1
    s_idx <- start_idx
    for(s in split_vec){
      series4 <- c(series4, s - last)
      whatshap_lens <- c(whatshap_lens, dfhwhathaps_rows[s_idx + s - last-1,"V2"]- dfhwhathaps_rows[s_idx,"V2"] + 1)
      s_idx <- s_idx + s - last
      last <- s
      
      
    }     
    series4 <- c(series4, fact - last -1)
    whatshap_lens <- c(whatshap_lens, dfhwhathaps_rows[s_idx + fact - last -2,"V2"]- dfhwhathaps_rows[s_idx,"V2"] + 1)
  }
  start_idx <- start_idx + fact
}

n50(whatshap_lens)
#12,542
sum(whatshap_lens)
#1,148,739,180
mean(whatshap_lens)
#4,369
min(whatshap_lens)
#2
max(whatshap_lens)
#183,560

sum(series4)
#[1] 2,262,328
max(series4)
#374
table(series4)
length(which(series4==1))
#[1] 0
#after splitting based on coninuous SNP, some haplotype with length 1 will generated, we remove them again from dfmerged4

dfmerged4$rowfactor <- rep(1:length(series4), series4)
dfmerged4$hlen <- rep(series4, series4)
dfmerged4$hlen_bp <- rep(whatshap_lens, series4)
dfmerged4 <- dfmerged4[dfmerged4$hlen>1,]
dim(dfmerged4) 
#[1] 2,262,328      10
#now calculate accuracy 
agg4<-aggregate(both~rowfactor, dfmerged4, FUN=function(b){switcherror(b)})

#find switch error rate per kbp for all haplotypes
agg4$len_bp <- whatshap_lens
agg4_long <- agg4[agg4$len_bp>=2,]
(sum(agg4_long$both)*1000)/sum(agg4_long$len_bp)
#0.00453976
length(which(whatshap_lens<1000))
#101331
length(which(whatshap_lens>=1000))
#161554

