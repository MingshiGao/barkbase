#!/usr/bin/env R
library(ggplot2)
library(reshape2)

file='signal.matrix'
df=read.table(file)
rownames(df)=df[,1]
l=df[,1]
df=df[,-1]
logdf=log10(df+0.1)
norm=as.data.frame(t(apply(logdf,1,scale)))
head(norm)
df=norm


set.seed(1)
Summary=data.frame("k"=double(),'tot.withinss'=double())
i=1
for( k in seq(5,20,1)){
  print(k)
  cluster=kmeans(df ,centers = k, iter.max = 1000, algorithm="Lloyd")
  Summary[i,1]=k
  Summary[i,2]=cluster$tot.withinss
  i=i+1
}
colnames(Summary)[2]="totwithinss"
p1=ggplot(Summary,aes(x=k,y=totwithinss))+geom_line(stat='identity')+geom_point()+
  scale_x_continuous(limits = c(5, 30),breaks=seq(5,30,1))+
  ggtitle('rOCRs accessibility elbow plot')+ theme_minimal()
png('kmeans.elbow.png',height=5,width=5,unit='in',res=150)
plot(p1)
dev.off()

k=13
cluster=kmeans(df ,centers = k, iter.max = 1000, algorithm="Lloyd")
df$cluster=cluster$cluster
sorted=df[order(df$cluster),]
out=as.data.frame(cbind(rownames(df),df$cluster))
write.table(out,'kmeans.cluster.list' ,col.names=F,row.names=F,sep="\t",quote=F)
