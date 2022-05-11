library(plotly)
library(gplots)
library(RColorBrewer)
library(sqldf)
library(DESeq2)
library(stringr)



#############################################################################
# Load the data
#############################################################################


cond <-read.csv("nanostring_cond.csv",stringsAsFactors = FALSE)
rownames(cond) <- cond$sample

read_cnt <- function(fname, panel){
  cnt <-read.csv(fname,stringsAsFactors = FALSE)
  cnt_lab <- cnt[,1:2]
  cnt <- cnt[,-(1:3)]
  cnt_lab$panel <- panel
  cnt
  rownames(cnt) <- sprintf("%s_%s",cnt_lab$panel,cnt_lab$symbol)
  colnames(cnt) <- rownames(cond)
  cnt_lab$measurement <- rownames(cnt)
  
  return(list(cnt=cnt, cnt_lab=cnt_lab, cond=cond))
}

dat_cart <- read_cnt("raw_cart.csv","cart")
dat_metabolism <- read_cnt("raw_metabolism.csv","metabolism")

# controller: a mouse that has the tumour under control



#############################################################################
# Normalization, deseq-style
#############################################################################

normalize_nstring <- function(dat_cnt){
  cnt <- dat_cnt$cnt
  cnt_lab <- dat_cnt$cnt_lab
  
  print(table(rowSums(cnt)>500))
  sf <- estimateSizeFactorsForMatrix(cnt[rowSums(cnt)>500,])
  ncnt <- t(t(cnt)/sf)
  #plot(log(ncnt[,1]),log(ncnt[,2]),pch=19,cex=0.5,col=df$the_col)
  log_cnt <- log2(ncnt)
  
  df <- data.frame(
    gene=cnt_lab$symbol,
    mean=(rowMeans(log_cnt)),
    sd=(apply(log_cnt,1,sd)),
    the_col="black",
    stringsAsFactors = FALSE
  )
  df$the_col[cnt_lab$gene_class=="Negative"] <- "red"
  df$the_col[cnt_lab$gene_class=="Positive"] <- "blue"
  
  plot(log_cnt[,1],log_cnt[,2],pch=19,cex=0.5,col=df$the_col)
  
  dat_cnt$ncnt <- ncnt
  dat_cnt$log_cnt <- log_cnt
  return(dat_cnt)
}

dat_cart <- normalize_nstring(dat_cart)
dat_metabolism <- normalize_nstring(dat_metabolism)


### What genes are in common between panels that are not controls?
control_probes <- dat_cart$cnt_lab$symbol[dat_cart$cnt_lab$gene_class %in% c("Positive","Negative")]
common_genes <- setdiff(intersect(dat_cart$cnt_lab$symbol,dat_metabolism$cnt_lab$symbol),control_probes)

#Normalize these together; repeated measurements should agree
corr_factor <- mean(dat_cart$log_cnt[sprintf("cart_%s",common_genes),] - dat_metabolism$log_cnt[sprintf("metabolism_%s",common_genes),])
dat_cart$log_cnt <- dat_cart$log_cnt - corr_factor

plot(dat_cart$log_cnt[sprintf("cart_%s",common_genes),],dat_metabolism$log_cnt[sprintf("metabolism_%s",common_genes),],pch=19,cex=0.1)

#Produce common objects
cnt_lab <- rbind(dat_cart$cnt_lab, dat_metabolism$cnt_lab)
log_cnt <- rbind(dat_cart$log_cnt, dat_metabolism$log_cnt)

#Ignore some genes
cnt_lab$gene_class[
  str_starts(cnt_lab$symbol,"TRAC") | 
    str_starts(cnt_lab$symbol,"TRAV") | 
    str_starts(cnt_lab$symbol,"TRB") | 
    str_starts(cnt_lab$symbol,"TRC") | 
    str_starts(cnt_lab$symbol,"TRD") |
  str_starts(cnt_lab$symbol,"TRG") 
] <- "ignore"


#############################################################################
# Check replicability after normalization
#############################################################################

df <- data.frame(
  gene=cnt_lab$symbol,
  mean=(rowMeans(log_cnt)),
  sd=(apply(log_cnt,1,sd)),
  col="black",
  stringsAsFactors = FALSE
)

#Show replicated measurements
#rep_genes <- common_genes#names(table(cnt_lab$symbol)[table(cnt_lab$symbol)>1])
#col <- rep("gray",nrow(cnt))
df$col[df$gene %in% common_genes] <- "red"
plot(df$mean, df$sd, col=df$col, pch=19)  

#Plot lines to show how repeatable the measuments are
df <- df[df$gene %in% common_genes,]
df <- df[order(df$gene),]
for(i in 1:(nrow(df)/2)){
  j <- (i-1)*2+1
  lines(df$mean[j:(j+1)], df$sd[j:(j+1)], col="blue")
}


#############################################################################
# Average technical reps
#############################################################################

new_lab <- unique(cnt_lab[,c("gene_class","symbol")])
new_lab$panel <- "merged"
new_lab$measurement <- sprintf("merged_%s",new_lab$symbol)

new_cnt <- matrix(NA, nrow=nrow(new_lab), ncol=ncol(log_cnt))
for(i in 1:nrow(new_lab)){
  new_cnt[i,] <- colMeans(log_cnt[cnt_lab$symbol == new_lab$symbol[i],,drop=FALSE])
}
rownames(new_cnt) <- new_lab$measurement

cnt_lab <- rbind(cnt_lab, new_lab)
log_cnt <- rbind(log_cnt, new_cnt)

#merged_AKT1  merged_CD247   merged_CD8A merged_CMKLR1  merged_CXCL9   merged_GNLY 
#2             2             2             2             2             2 

#############################################################################
# Estimate variance - bio reps
#############################################################################

# Estimate variance
#fit a linear model ... estimate sd for each point
# log(x / y) = log10(x) - log10(y)      +-    sqrt(2) * sd(for one point ... same as gene)  use the approximation
#   p-values by using normal assumption... 

# df <- data.frame(
#   mean=rowMeans(log_cnt),
#   sd=apply(log_cnt,1,sd)
# )

#df$mean2 <- 1/df$mean#)#*df$mean
#df$mean3 <- df$mean*df$mean*df$mean

df <- data.frame(
  gene=cnt_lab$symbol,
  mean=(rowMeans(log_cnt)),
  sd=(apply(log_cnt,1,sd)),
  col="black",
  stringsAsFactors = FALSE
)

df$mean2 <- df$mean*df$mean
df$mean3 <- df$mean*df$mean*df$mean
df$mean4 <- df$mean*df$mean*df$mean*df$mean

df2 <- df#[df$sd<2 & df$mean<12,]
thelm <- lm(sd~mean+mean2, df2)
plot(df2$mean, df2$sd, col=df2$col, pch=19)  
points(df2$mean, thelm$fitted.values, col="blue")

cnt_lab$est_sd <- thelm$fitted.values


#############################################################################
# Alternative way, focus on technical replicability
#############################################################################

# reps <- data.frame(
#   a=log_cnt[cnt_lab$symbol %in% common_genes & cnt_lab$panel=="cart"],
#   b=log_cnt[cnt_lab$symbol %in% common_genes & cnt_lab$panel=="metabolism"])
# reps$mean <- rowMeans(reps)
# reps$diff <- abs(reps$a-reps$b)
# 
# thelm <- lm(diff ~ mean,reps)
# plot(reps$mean, abs(reps$diff))  # sqrt(2) sigma
# points(reps$mean, thelm$fitted.values,col="blue")
# 
# cnt_lab$est_sd <- (thelm$coefficients[1] + rowMeans(log_cnt)*thelm$coefficients[2])/sqrt(2)






#############################################################################
# Run comparisons
#############################################################################

## to compare groupA/groupB
fun_comp <- function(set_groupA, set_groupB){
  #Fit and estimate sd's
  # df <- data.frame(
  #   mean=rowMeans(log_cnt),
  #   sd=apply(log_cnt,1,sd)
  # )
  # thelm <- lm(sd~mean, df)
  # df$sd_est <- thelm$fitted.values

  # Calculate means in both groups to compare
  a <- rowMeans(log_cnt[,set_groupA,drop=FALSE])
  b <- rowMeans(log_cnt[,set_groupB,drop=FALSE])
  sd_est <- sqrt(sqrt(length(set_groupA)) + sqrt(length(set_groupB))) * cnt_lab$est_sd   ###df$sd_est
  
  #zscore:  x-u / sd
  #Var(x-y) = 2*var[x]   ... sd[x-y] = sqrt(2)*sd[x]
  #more generally: two groups of size n vs size m.   then variance is sqrt(sqrt(n) + sqrt(m))*sd(x)
  
  # Object to return with fold change and z-score
  comp <- data.frame(
    measurement = cnt_lab$measurement,
    gene = cnt_lab$symbol,
    panel = cnt_lab$panel,
    fc = a-b,
    z = (a-b) / sd_est,
    gene_class = cnt_lab$gene_class,
    stringsAsFactors = FALSE
  )
  
  comp$pvalue2sided <- 2*pnorm(-abs(comp$z))
  
  comp <- comp[order(abs(comp$z), decreasing = TRUE),]
  return(comp)
}
# 
# #### Plot genes as a scatter plot
# plot_topgenes_xy <- function(res, set_groupA, set_groupB, n=10, 
#                              xlab=do.call(paste, c(as.list(set_groupA), sep = ",")), 
#                              ylab=do.call(paste, c(as.list(set_groupB), sep = ","))){
#   x <- rowMeans(log_cnt[,set_groupA,drop=FALSE])
#   y <- rowMeans(log_cnt[,set_groupB,drop=FALSE])
#   
#   sig_measurement <- res$measurement[res$pvalue2sided < 1e-3 & !(res$gene_class %in% c("Negative","Positive"))]
#   print(sig_measurement)
#   
#   col <- rep("gray",nrow(res))
#   col[cnt_lab[rownames(res),]$gene_class %in% c("Negative")] <- "blue"
#   #col[cnt_lab$gene_class %in% c("Positive")] <- "red"
#   keep <- cnt_lab[rownames(res),]$gene_class %in% c("Negative","Endogenous")
#   plot(x[keep],y[keep], pch=19, cex=0.3, xlab=xlab, ylab=ylab, col=col[keep])
#   
#   #col[rownames(log_cnt) %in% sig_measurement] <- "red"
#   if(length(sig_measurement)>0){
#     col[rownames(log_cnt) %in% sig_measurement] <- "sig"
#     if(any(col=="sig")){
#       text(x[col=="sig"], y[col=="sig"], labels = cnt_lab$symbol[col=="sig"],cex=0.5)
#     }
#   }
# }



#### Plot genes as a scatter plot
plot_topgenes_xy <- function(res, set_groupA, set_groupB, n=10, 
                             xlab=do.call(paste, c(as.list(set_groupA), sep = ",")), 
                             ylab=do.call(paste, c(as.list(set_groupB), sep = ","))){
  
  coords <- data.frame(
    row.names = rownames(log_cnt),
    measurement = rownames(log_cnt),
    x = rowMeans(log_cnt[,set_groupA,drop=FALSE]),
    y = rowMeans(log_cnt[,set_groupB,drop=FALSE])
  )
  coords <- merge(coords,cnt_lab)
  coords <- coords[coords$gene_class %in% c("Negative","Endogenous"),]
  coords$col <- "gray"
  coords[coords$gene_class %in% c("Negative"),]$col <- "blue"
  plot(coords$x, coords$y, pch=19, cex=0.3, xlab=xlab, ylab=ylab, col=coords$col)
  
  sig_measurement <- res$measurement[res$pvalue2sided < 1e-3 & !(res$gene_class %in% c("Negative","Positive"))]
  print(sig_measurement)
  
  coords$sig <- FALSE
  coords$sig[coords$measurement %in% sig_measurement] <- TRUE
  if(any(coords$sig)){
    coords <- coords[coords$sig,]
    text(coords$x, coords$y, labels = coords$symbol,cex=0.5)
  }
  
}


#rownames(cnt_lab) <- sprintf("merged_%s",cnt_lab$measurement)

#res <- fun_comp(c("M28z_day0_BC"),c("M28z_sick_M14"))
#plot_topgenes_xy(res,c("M28z_day0_BC"),c("M28z_sick_M14"))

### Compare all groups
tocomp <- read.csv("comparisons5.csv", stringsAsFactors = FALSE)
allres <- list()
for(i in 1:nrow(tocomp)){
  grp_a <- colnames(tocomp)[tocomp[i,]=="a"]
  grp_b <- colnames(tocomp)[tocomp[i,]=="b"]
  print(i)
  print(tocomp[i,1])

  res <- fun_comp(grp_a, grp_b)
  res <- res[res$gene_class=="Endogenous",]
  #write.csv(res,sprintf("out/all_%s.csv",tocomp[i,1]))

  pdf(sprintf("out/metabolismpanel_%s.pdf",tocomp[i,1]))
  plot_topgenes_xy(res[res$panel=="metabolism",],grp_a, grp_b,n=20)
  dev.off()
  write.csv(res[res$panel=="metabolism",],sprintf("out/metabolismpanel_%s.csv",tocomp[i,1]))
  
  pdf(sprintf("out/cartpanel_%s.pdf",tocomp[i,1]))
  plot_topgenes_xy(res[res$panel=="cart",],grp_a, grp_b,n=20)
  dev.off()
  write.csv(res[res$panel=="cart",],sprintf("out/cartpanel_%s.csv",tocomp[i,1]))
  
  pdf(sprintf("out/merged_%s.pdf",tocomp[i,1]))
  plot_topgenes_xy(res[res$panel=="merged",],grp_a, grp_b,n=20)
  dev.off()
  write.csv(res[res$panel=="merged",],sprintf("out/merged_%s.csv",tocomp[i,1]))
  
  allres[[i]] <- res
}
#tocomp



#sort(rowMeans(log_cnt[cnt_lab$gene_class =="Negative",]))
#sort(rowMeans(log_cnt[cnt_lab$gene_class =="Positive",]))



#############################################################################
# Heatmap of all comparisons
#############################################################################

m <- matrix(NA, ncol=length(allres), nrow=nrow(allres[[1]]))
rownames(m) <- allres[[1]]$gene
colnames(m) <- tocomp$name_of_comparison
for(i in 1:nrow(tocomp)){
  m[allres[[i]]$gene,i] <- allres[[i]]$z
}
m <- na.omit(m)

#z = 3 is 0.1%
redm <- m[apply(abs(m)>3,1,any),]
redm <- redm[rownames(redm) %in% allres[[i]]$gene[allres[[i]]$gene_class=="Endogenous"],]
dim(redm)
pdf("zheatmap_de.pdf",width = 4)
my_palette <- colorRampPalette(c("red", "white", "green"))(n = 299)
heatmap.2(redm,cexRow=0.3,cexCol = 0.5,
          density.info="none",  # turns off density plot inside color legend
          col=my_palette, 
          trace="none"         # turns off trace lines inside the heat map
)
dev.off()


#############################################################################
# Heatmap of all vs D0
#############################################################################

### Compare all samples to day0
allres <- list()
allres_vs <- list()
j<-1
for(i in 1:nrow(cond)){
  if(cond$category[i]!="day0"){
    grp_a <- rownames(cond)[cond$category=="day0"]
    grp_b <- rownames(cond)[i]
    print(i)
    res <- fun_comp(grp_a, grp_b)
    res <- res[res$panel=="merged" & res$gene_class=="Endogenous",]
    allres[[j]] <- res
    allres_vs[[j]] <- grp_b
    j <- j+1
  }
}


m <- matrix(NA, ncol=length(allres), nrow=nrow(allres[[1]]))
rownames(m) <- allres[[1]]$gene
colnames(m) <- unlist(allres_vs)
for(i in 1:ncol(m)){
  m[allres[[i]]$gene,i] <- allres[[i]]$z
}
m <- na.omit(m)


#redm <- m[apply(abs(m)>3,1,any),]  #z = 3 is 0.1%
redm <- m[apply(abs(m)>1,1,any),]
redm <- redm[rownames(redm) %in% allres[[i]]$gene[allres[[i]]$gene_class=="Endogenous"],]
dim(redm)
pdf("zheatmap_vsday0.pdf",width = 4, height = 150)
my_palette <- colorRampPalette(c("red", "white", "green"))(n = 299)
heatmap.2(redm,cexRow=0.3,cexCol = 0.5,
          density.info="none",  # turns off density plot inside color legend
          col=my_palette, 
          trace="none"         # turns off trace lines inside the heat map
)
dev.off()




#############################################################################
# Dimensional reduction to compare all samples
#############################################################################



library(umap)



cond$color <- "black"
cond$color[cond$category=="day0"] <- "darkgreen"
cond$color[cond$category=="interim"] <- "darkgray"
cond$color[cond$category=="startsick"] <- "black"
cond$color[cond$category=="sick"] <- "black"
um <- umap(t(log_cnt[cnt_lab$gene_class=="Endogenous" & rowMeans(log_cnt)>3,]),colvec=c('skyblue'),n_neighbors=10)
pdf("out/umap.pdf")
plot(um$layout,pch=19,cex=0.5,xlim=c(-2,2),ylim=c(-2,2), xlab="UMAP1",ylab="UMAP2",  col=cond$color)
text(x=um$layout[,1], y=um$layout[,2]+.1, labels = rownames(um$layout),cex=0.8, col=cond$color)
dev.off()





#############################################################################
# IVIS & SMRP
#############################################################################

library(limma)

## Best on log scale
plot(log10(na.exclude(cond$SMRP)))
## Best on linear scale
plot((na.exclude(cond$ivis_sac)))
plot((na.exclude(cond$ivis_pre)))


#### ivis sac
keep <- !is.na(cond$ivis_sac)
cond_red <- cond[keep,]
cnt_red <- log_cnt[startsWith(rownames(log_cnt),"merged_"),keep]
phenoData <- new("AnnotatedDataFrame",data=cond_red)
eset <- ExpressionSet(
  assayData=cnt_red,
  phenoData=phenoData)

design <- model.matrix(~1 + ivis_sac, cond_red)
limma_fit <- lmFit(eset, design)
fit <- eBayes(limma_fit)
res <- topTable(fit, number = 100000)
res_ivis_sac <- res

gene <- "merged_cart_DDIT4"
gene <- "merged_cart_CXCL13"
gene <- "merged_cart_CXCR4"

plot(cnt_red[gene,], cond_red$ivis_sac,cex=0)
text(cnt_red[gene,], cond_red$ivis_sac, rownames(cond_red),cex=0.5)



#### ivis pre
keep <- !is.na(cond$ivis_pre)
cond_red <- cond[keep,]
cnt_red <- log_cnt[startsWith(rownames(log_cnt),"merged_"),keep]
phenoData <- new("AnnotatedDataFrame",data=cond_red)
eset <- ExpressionSet(
  assayData=cnt_red,
  phenoData=phenoData)

design <- model.matrix(~1 + ivis_pre, cond_red)
limma_fit <- lmFit(eset, design)
fit <- eBayes(limma_fit)
res <- topTable(fit, number = 100000)
res_ivis_pre <- res

gene <- "merged_cart_DDIT4"
gene <- "merged_cart_CXCL13"
gene <- "merged_cart_CXCR4"

plot(cnt_red[gene,], cond_red$ivis_pre,cex=0)
text(cnt_red[gene,], cond_red$ivis_pre, rownames(cond_red),cex=0.5)
#HIF1A comes up high
#NEG_H high :(



#### ivis pre
cond$log_SMRP <- log10(cond$SMRP)
keep <- !is.na(cond$log_SMRP)
cond_red <- cond[keep,]
cnt_red <- log_cnt[startsWith(rownames(log_cnt),"merged_"),keep]
phenoData <- new("AnnotatedDataFrame",data=cond_red)
eset <- ExpressionSet(
  assayData=cnt_red,
  phenoData=phenoData)

design <- model.matrix(~1 + log_SMRP, cond_red)
limma_fit <- lmFit(eset, design)
fit <- eBayes(limma_fit)
res <- topTable(fit, number = 100000)
res_SMRP <- res

gene <- "merged_metabolism_FLT1"
gene <- "merged_cart_DDIT4"

plot(cnt_red[gene,], cond_red$log_SMRP,cex=0)
text(cnt_red[gene,], cond_red$log_SMRP, rownames(cond_red),cex=0.5)
#FLT1 high here too


#####

limmares <- merge(
  merge(
    data.frame(
      measurement=rownames(res_SMRP),
      res_SMRP$logFC,
      res_SMRP$P.Value),
    data.frame(
      measurement=rownames(res_ivis_sac),
      res_ivis_sac$logFC,
      res_ivis_sac$P.Value)
    ),
  data.frame(
    measurement=rownames(res_ivis_pre),
    res_ivis_pre$logFC,
    res_ivis_pre$P.Value)
  )
limmares$measurement <- str_sub(limmares$measurement,8)
limmares <- limmares[limmares$measurement %in% cnt_lab$measurement[cnt_lab$gene_class!="ignore"],]

write.csv(limmares, "limma/limmalist.csv")

pdf("limma/limma_1.pdf")
plot(limmares$res_SMRP.logFC, limmares$res_ivis_sac.logFC,cex=0)
text(limmares$res_SMRP.logFC, limmares$res_ivis_sac.logFC, labels = limmares$measurement,cex=0.3)
dev.off()
fig <- plot_ly(data = limmares, x = ~res_SMRP.logFC, y = ~res_ivis_sac.logFC, text = ~measurement)
htmlwidgets::saveWidget(partial_bundle(fig), "limma/limma_1.html",selfcontained = T)

pdf("limma/limma_2.pdf")
plot(limmares$res_SMRP.logFC, limmares$res_ivis_pre.logFC,cex=0)
text(limmares$res_SMRP.logFC, limmares$res_ivis_pre.logFC, labels = limmares$measurement,cex=0.3)
dev.off()
fig <- plot_ly(data = limmares, x = ~res_SMRP.logFC, y = ~res_ivis_pre.logFC, text = ~measurement)
htmlwidgets::saveWidget(partial_bundle(fig), "limma/limma_2.html")

pdf("limma/limma_3.pdf")
plot(limmares$res_ivis_pre.logFC, limmares$res_ivis_sac.logFC,cex=0)
text(limmares$res_ivis_pre.logFC, limmares$res_ivis_sac.logFC, labels = limmares$measurement,cex=0.3)
dev.off()
fig <- plot_ly(data = limmares, x = ~res_ivis_pre.logFC, y = ~res_ivis_sac.logFC, text = ~measurement)
htmlwidgets::saveWidget(partial_bundle(fig), "limma/limma_3.html")




