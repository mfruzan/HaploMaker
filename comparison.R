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
dfvcf <- read.table("commonmerge_noheader.vcf", header = F, colClasses = c("character","numeric","NULL","NULL","NULL","NULL","NULL","NULL","NULL","character"))
dfvcf$snp_num <- c(1:nrow(dfvcf))
dim(dfvcf)
#1,890,198  4
############################## HaploMiner
dfhaplominer <-read.table("haplominer_commonmerge.hap", colClasses = c("character", "numeric","character","character","character","character"), header = F, fill = T)
dfhaplominer[dfhaplominer=='']<-NA
dfhaplominer_rows <- dfhaplominer[!is.na(dfhaplominer$V6), ]
factors <- dfhaplominer[is.na(dfhaplominer$V6), 4]
dfhaplominer_rows$rownum <- c(1:nrow(dfhaplominer_rows))
dfmerged <-merge(dfhaplominer_rows,dfvcf, by=c(1,2), all.x=T)
dfmerged <- dfmerged[order(dfmerged$rownum),]
dup_idx <- which(duplicated(dfmerged$rownum))
dfmerged <- dfmerged[-dup_idx,]
rowfactor <- rep(c(1:length(factors)), factors)
dfmerged$rowfactor <- rowfactor
dfmerged$both <- paste(dfmerged$V6,'_',dfmerged$V10,sep="")
agg <- aggregate(both ~ rowfactor, dfmerged, FUN=function(b)switcherror(b))
# Estimate N50
tmp <- dfhaplominer[is.na(dfhaplominer$V6), 2:3]
haplominer_lens <- as.numeric(tmp$V3)-as.numeric(tmp$V2) + 1
n50(haplominer_lens)
#326
mean(haplominer_lens)
# 251
max(haplominer_lens)
# 2770

now estimate haplominer swtich rate per kbp of long haplotype
length(which(haplominer_lens<1000))
#320142
length(which(haplominer_lens>=1000))
#2318
agg$len_bp <- haplominer_lens
agg_long <- agg[agg$len_bp>=2,]
(sum(agg_long$both)*1000)/sum(agg_long$len_bp)
#0.0321 
############################### HapCompass
dfhapcompass <- read.table("hapcompass_commonmerge_q5.hap", colClasses=c("character","numeric","numeric","numeric","numeric","NULL","character"), header=F, fill=T)
#1,173,435
dfhapcompass[dfhapcompass=='']<-NA 
summary_vec <-which(!is.na(dfhapcompass$V7))
summary <-dfhapcompass[which(dfhapcompass$V1=="BLOCK"),c("V2","V3")]
hapcompass_lens <- summary$V3 - summary$V2 +1 
factors2<- c(summary_vec[-1], (length(dfhapcompass$V1)+1)) - summary_vec - 1
chr_vec <- dfhapcompass$V7[!is.na(dfhapcompass$V7)]
dfhapcompass_rows <- dfhapcompass[is.na(dfhapcompass$V7) ,1:5]
dfhapcompass_rows$V1 <- rep(chr_vec, factors2)
dfhapcompass_rows$V8 <- paste(dfhapcompass_rows$V4,"|", dfhapcompass_rows$V5, sep="")
dfhapcompass_rows$rownum <- c(1:nrow(dfhapcompass_rows))
dfmerged2 <-merge(dfhapcompass_rows,dfvcf, by=c(1,2), all.x=T)
dfmerged2 <- dfmerged2[order(dfmerged2$rownum),]
dup_idx <- which(duplicated(dfmerged2$rownum))
dfmerged2 <-dfmerged2[-dup_idx,]
rowfactor <- rep(c(1:length(factors2)), factors2)
dfmerged2$rowfactor <- rowfactor
dfmerged2$both <- paste(dfmerged2$V8,'_',dfmerged2$V10,sep="")
agg2<-aggregate(both~rowfactor, dfmerged2, FUN=function(b){switcherror(b)})
#now estimate hapcompass swtich rate per kbp of long haplotype
length(which(hapcompass_lens<1000))
#330295
length(which(hapcompass_lens>=1000))
#2829
agg2$len_bp <- hapcompass_lens
agg2_long <- agg2[agg2$len_bp>=2,]
(sum(agg2_long$both)*1000)/sum(agg2_long$len_bp)
#0.1439

# As HapCompass reports non-continuous SNPs as haplotype: we require code 
# to detect only continuous ones
v<-dfhapcompass_rows$V3
series<-numeric(length=0)
hapcompass_lens <- numeric(length=0)
start_idx <- 1
for(fact in factors2){
  
  split_vec <- numeric(length=0)
  for (i in 0:(fact-2)){
    if(v[start_idx+i+1]-v[start_idx+i] >1){
      split_vec <- c(split_vec,i)
    }
  }
  if(length(split_vec)==0){
    series <- c(series, fact)
    hapcompass_lens <- c(hapcompass_lens, dfhapcompass_rows[start_idx + fact-1,"V2"]- dfhapcompass_rows[start_idx,"V2"] + 1)
  }else{
    last <- -1
    s_idx <- start_idx
    for(s in split_vec){
      series <- c(series, s - last)
      hapcompass_lens <- c(hapcompass_lens, dfhapcompass_rows[s_idx + s - last-1,"V2"]- dfhapcompass_rows[s_idx,"V2"] + 1)
      s_idx <- s_idx + s - last
      last <- s
      
      
    }     
    series <- c(series, fact - last -1)
    hapcompass_lens <- c(hapcompass_lens, dfhapcompass_rows[s_idx + fact - last -2,"V2"]- dfhapcompass_rows[s_idx,"V2"] + 1)
  }
  start_idx <- start_idx + fact
}
# Estimate N50
n50(hapcompass_lens)
# 329
mean(hapcompass_lens)
# 246
max(hapcompass_lens)
# 2770
############################### HapCUT2
dfhapcut <- read.table("hapcut2_commonmerge_q5.hap",sep = "\t" , fill=T)
dfhapcut[dfhapcut=='']<-NA
dfhapcut_rows <- dfhapcut[!is.na(dfhapcut$V2) ,c(4,5,2,3)]
dfhapcut_rows$V4 <- as.character(dfhapcut_rows$V4)
dfhapcut_rows$V5 <- as.numeric(dfhapcut_rows$V5)
summary_vec <-which(dfhapcut$V1=="******** ")
factors3<-  summary_vec - c(0, summary_vec[-length(summary_vec)])  - 2
factors3 <- c(factors3,2)
sum(factors3)
nrow(dfhapcut_rows)
#667676
dfhapcut_rows$V6 <- paste(dfhapcut_rows$V2,"|", dfhapcut_rows$V3, sep="")
dfhapcut_rows$rownum <- c(1:nrow(dfhapcut_rows))
dfmerged3 <-merge(dfhapcut_rows,dfvcf, by=c(1,2), all.x=T)
dfmerged3 <- dfmerged3[order(dfmerged3$rownum),]
dup_idx <- which(duplicated(dfmerged3$rownum))
rowfactor <- rep(c(1:length(factors3)), factors3)
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

#now estimate hapcut swtich rate per kbp of long haplotype
length(which(hapcut_lens<1000))
#270,488
length(which(hapcut_lens>=1000))
#1721
agg3$len_bp <- hapcut_lens
agg3_long <- agg3[agg3$len_bp>=2,]
(sum(agg3_long$both)*1000)/sum(agg3_long$len_bp)
#0.0376

# As HapCUT2 reports non-continuous SNPs as haplotype: we require code 
# to detect only continuous ones
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
# Estimate N50
n50(hapcut_lens)
# 313
mean(hapcut_lens)
# 224
max(hapcut_lens)
# 2581

############################### WhatsHap

dfwhatshap<-read.table("whatshap_phased_commonmerge.vcf", header=F, colClasses=c("character","numeric","NULL","NULL","NULL","NULL","NULL","NULL","NULL","character"))
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
#[1] 828,878      4
#then we realize that whatshap reports many haplotypes with length 1 SNP!!!
hap_lens<-aggregate(V1~rowfactor, dfhwhathaps_rows, FUN=length)
length(which(hap_lens$V1==1))
#[1] 8957  haplotypes are with length 1
#now filter rows that contain hap with length 1
dfhwhathaps_rows <- dfhwhathaps_rows[!dfhwhathaps_rows$rowfactor %in% hap_lens[hap_lens$V1==1,"rowfactor"] , ]
dim(dfhwhathaps_rows)
# 815,976      4
dfhwhathaps_rows$rownum <- c(1:nrow(dfhwhathaps_rows))
dfmerged4 <-merge(dfhwhathaps_rows,dfvcf, by=c(1,2), all.x=T)
dfmerged4 <- dfmerged4[order(dfmerged4$rownum),]
dup_idx <- which(duplicated(dfmerged4$rownum))
#rowfactor <- rep(c(1:length(factors3)), factors3)
#dfmerged3$rowfactor <- rowfactor
dfmerged4$both <- paste(dfmerged4$V3,'_',dfmerged4$V10,sep="")

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
#329
mean(whatshap_lens)
#252
max(whatshap_lens)
#2770

#after splitting based on continuous SNP, some haplotype with length 1 will be generated, we remove them again from dfmerged4

dfmerged4$rowfactor <- rep(1:length(series4), series4)
dfmerged4$hlen <- rep(series4, series4)
dfmerged4 <- dfmerged4[dfmerged4$hlen>1,]
#now calculate accuracy 
agg4<-aggregate(both~rowfactor, dfmerged4, FUN=function(b){switcherror(b)})
#find switch error rate per kbp for all haplotypes
agg4$len_bp <- whatshap_lens
agg4_long <- agg4[agg4$len_bp>=2,]
(sum(agg4_long$both)*1000)/sum(agg4_long$len_bp)
#0.0842