d <- read.table("snp_index_SNPdensity.bed", header=F, stringsAsFactor=T, col.names=c("chr", "start", "stop", "id", "Density"))
d$chr <- factor(d$chr, levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'))
d <- with(d, d[order(chr),])
library(scales)
library(ggplot2)
library(RColorBrewer)

display.brewer.all()
cbp1 <- c("#f3f3f0", "#0f6207", "#a9c904", "#e3eb12","#fca102","#f76f03","#f32502")

p <- ggplot(data=d, aes(x=start, y=1)) +
  facet_grid(chr ~ ., switch='y',scales='free_y', space='free_y') +
  geom_tile(aes(fill=Density)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 180)) +
  scale_fill_gradientn(colours = cbp1, breaks=seq(min(d$Density),max(d$Density),(max(d$Density)-min(d$Density))/8))+
  scale_x_continuous(breaks = c(0,30000000,60000000,90000000,120000000,150000000,180000000,210000000,240000000,270000000,300000000,330000000,360000000),
                     label = c("0Mb", "30Mb", "60Mb","90Mb", "120Mb", "150Mb","180Mb","210Mb","240Mb","270Mb","300Mb","330Mb","360Mb"))
print(p)
#####################
