options(stringsAsFactors=FALSE)
require(sqldf)
gt = read.table('~/R/nuvolone-ko-mice-genotypes.txt',sep='\t',header=TRUE)

colnames(gt) = tolower(colnames(gt)) # lower case column names

color_b6    = '#222222'
color_129   = '#975D2F' # 
color_het   = '#714623' # halfway between above
color_lines = '#FC1501' # gummi red - http://www.december.com/html/spec/color1.html
color_bg    = '#555555' # background color

value_b6  = 1
value_129 = 0
value_het = 0.5

# get lengths of each chromosome (highest SNP pos is fine)
chrlen = sqldf("
               select   chr, max(position) maxpos
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

# test plotting the chromosome boundaries
par(mfrow=c(1,1))
plot(NA,NA,xlim=c(0,max(gt$abspos)),ylim=c(0,1),xaxs='i',axes=FALSE)
abline(v=chrbreaks)
mtext(side=1,at=chrmids,text=chrlen$chr)

colnames(gt)

plot_order = c(5,6,7,10,9,8,11,12,13,14,15) # which columns to plot, and in what order
display_names = colnames(gt)
# placeholder - can add nicer-looking names later

pdf('~/R/ko-mice-genotypes.pdf',width=10,height=5)
par(mfrow=c(length(plot_order)+1,1),mar=c(.2,5,.2,1),oma=c(2,3,3,1))
iteration = 1
for (colno in plot_order) {
  plot(NA,NA,xlim=range(gt$abspos),ylim=c(0,1),yaxs='i',xaxs='i',
       yaxt='n',xaxt='n',ylab='',xlab='',cex.lab=2)
  mtext(side=2,text=paste(display_names[colno],"  "),las=1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=color_bg)
  
  het = gt[,colno]==value_het & !is.na(gt[,colno])
  if (sum(het) > 0) {
    points(gt$abspos[het],rep(1,sum(het)),col=color_het,type='h',lwd=2,lend=2)
  }
  b6 = gt[,colno]==value_b6 & !is.na(gt[,colno])
  if (sum(b6) > 0) {
    points(gt$abspos[b6],rep(1,sum(b6)),col=color_b6,type='h',lwd=2,lend=2)
  }
  m129 = gt[,colno]==value_129 & !is.na(gt[,colno])
  if (sum(m129) > 0) {
    points(gt$abspos[m129],rep(1,sum(m129)),col=color_129,type='h',lwd=2,lend=2)
  }
  axis(side=1,at=c(1,chrbreaks),labels=NA)
  if (iteration == length(plot_order)) {
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
# mtext(side=3,outer=TRUE,text='SNP genotyping of PrP knockout mice and controls',padj=-1)

dev.off()
