options(stringsAsFactors=FALSE)
require(sqldf)
gt = read.table('~/d/j/cureffilab/media/2014/08/nuvolone-ko-mice-genotypes.txt',sep='\t',header=TRUE)

gt$chr[which(gt$edbg==1)]

color_b6 = '#222222'
color_129 = '#B28647' # cappucino - http://www.december.com/html/spec/color1.html
color_het = '#685434' # halfway between above
color_lines = '#FC1501' # gummi red - http://www.december.com/html/spec/color1.html
color_bg = '#AAAAAA' # background color

value_b6 = 1
value_129 = 0
value_het = 0.5

# remove SNPs that are not B6/129 informative
gt = gt[gt$differentiate=="",]

# get lengths of each chromosome (highest SNP pos is fine)
chrlen = sqldf("
select   chr, max(pos) maxpos
from     gt
group by 1
order by 1
;")
# properly number and sort the chromosomes
chrlen$chrno = 0
chrlen$chrno[chrlen$chr=='X'] = 23
chrlen$chrno[chrlen$chr!='X'] = as.integer(chrlen$chr[chrlen$chr!='X'])
# calculate the "absolute position" where each chromosome starts
chrlen = chrlen[with(chrlen, order(chrno)),]
chrlen$startpos = cumsum(as.numeric(chrlen$maxpos)) - chrlen$maxpos + 1
# assign "absolute position" to each position
gt$abspos = gt$pos + chrlen$startpos[match(gt$chr,chrlen$chr)]
# double check that absolute positions are monotonically increasing
any(gt$abspos != cummax(gt$abspos)) # FALSE is good
# find which SNPs are the highest on their chromosome - for plotting vertical lines at chr breaks
gt$chrmax = gt$pos == chrlen$maxpos[match(gt$chr,chrlen$chr)]
# find places to plot vertical lines between chromosomes
chrbreaks = c(0,gt$abspos[gt$chrmax])
# find places to label chromosomes - midpoint between each line
chrmids = (chrbreaks[1:(length(chrbreaks)-1)] + chrbreaks[2:length(chrbreaks)]) / 2
# test:
plot(NA,NA,xlim=range(gt$abspos),ylim=c(0,1))
abline(v=chrbreaks)
mtext(side=1,at=chrmids,text=chrlen$chr)

plot(gt$abspos,gt$zh1,type='h',col='#000000')

plot_order = c(8,9,10,11,12,13)
display_names = colnames(gt)
display_names[plot_order] = c("C57BL/6J","B6 WT","ZrchIII","ZrchI","129/Ola","Edbg")

pdf('~/d/j/cureffilab/media/2014/08/ko-mice-genotypes-20140806.pdf',width=8,height=5)
par(mfrow=c(7,1),mar=c(.2,5,.2,1),oma=c(2,3,3,1))
iteration = 1
for (colno in plot_order) {
  plot(NA,NA,xlim=range(gt$abspos),ylim=c(0,1),yaxs='i',xaxs='i',
       yaxt='n',xaxt='n',ylab='',xlab='',cex.lab=2)
  mtext(side=2,text=paste(display_names[colno],"  "),las=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=color_bg)
  m129 = gt[,colno]==value_129 & !is.na(gt[,colno])
  if (sum(m129) > 0) {
      points(gt$abspos[m129],rep(1,sum(m129)),col=color_129,type='h',lwd=2,lend=2)
  }
  het = gt[,colno]==value_het & !is.na(gt[,colno])
  if (sum(het) > 0) {
      points(gt$abspos[het],rep(1,sum(het)),col=color_het,type='h',lwd=2,lend=2)
  }
  b6 = gt[,colno]==value_b6 & !is.na(gt[,colno])
  if (sum(b6) > 0) {
    points(gt$abspos[b6],rep(1,sum(b6)),col=color_b6,type='h',lwd=2,lend=2)
  }
  axis(side=1,at=c(1,chrbreaks),labels=NA)
  if (iteration == 6) {
      mtext(side=1,at=chrmids,text=chrlen$chr,cex=.6)  
  } 
  iteration = iteration + 1
}
# seventh plot is borderless, to hold legend and scale bar
plot(NA,NA,type='l',lend=2,lwd=3,axes=FALSE,
     xlim=range(gt$abspos),ylim=c(0,1),yaxs='i',xaxs='i',
     ylab='',xlab='')
arrows(x0=1e7,x1=6e7,y0=.5,y1=.5,angle=90,code=3,lwd=3,length=.05,lend=1)
text(6e7,.5,pos=4,labels='50 Mb')
legend('right',ncol=4,c("B6","Het","129","No coverage"),
       col=c(color_b6,color_het,color_129,color_bg),
       pch=NA,lwd=5,bty='n')
mtext(side=3,outer=TRUE,text='SNP genotyping of PrP knockout mice and controls',padj=-1)

dev.off()

