library(ggplot2)
library(readr)
library(reshape2)
library(tidyr)

files <- Sys.glob("*.coverage.csv") #code needs to be in directory with files
SEGMENTS = c('MR766')
segment_df = data.frame(segment=as.character(),totalcount=as.numeric())
segment_name = data.frame(name=as.character(),segment=as.character(),ref = as.character(), totalcount=as.numeric())
metadata = '../MetadataZika/Zika_Metadata_v3.csv'
meta_df = read.csv(file=metadata,header=T,sep=",",na.strings = c('nan'))
for (filename in files) {
  print(filename)
  mydata1 <- read.csv(file=filename,header=T,sep=",",na.strings = c(''))
  print(head(mydata1))
  mydata=merge(mydata1,meta_df,by.x='name',by.y='name')
  mydata$segment = factor(mydata$segment,levels=SEGMENTS)
  
  x = aggregate(totalcount~segment, data=mydata, mean)
  segment_df = rbind(x,segment_df)
  y = aggregate(totalcount~ name + segment + ref, data=mydata,mean)
  print(y)
  segment_name = rbind(y,segment_name)

  x <- ggplot(data = mydata, aes(x=ntpos, y=totalcount, colour=(totalcount>=200))) + 
    geom_line(aes(group=1)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10), 
          axis.text.y = element_text(size=10)) +
    scale_color_manual(values=c('red','black')) + 
    facet_wrap(.~id2)
    #print(x)
    ggsave(x, filename = paste(filename, ".pdf", sep=""), width = 40, height = 40, limitsize=FALSE)
}
head(segment_df)
write.csv(segment_df, file='Segment.Average.csv',row.names = F)

head(segment_name)
write.csv(segment_name, file='Sample.Average.csv', row.names=F)

#Plots for the average coverage per segment and sample: 
head(segment_name)
plot_segment = ggplot(segment_name,aes(x=segment,y=totalcount)) +
  geom_boxplot() +
  geom_point() + 
  #facet_grid(ref~.) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=-90, hjust=1))
print(plot_segment)
ggsave(plot_segment, filename = "AvgCoverage.Segment.pdf", width = 15, height = 8, limitsize=FALSE)
