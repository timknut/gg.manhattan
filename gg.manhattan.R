gg.manhattan = function(dataframe, title=NULL, max.y="max", suggestiveline=1, genomewideline=-log10(5e-8), 
								size.x.labels=9, size.y.labels=10, annotate=F, SNPlist=NULL, logp = TRUE) {
	if(!require(ggplot2)) install.packages("ggplot2")
	if(!require(dplyr)) install.packages("dplyr")
	if (annotate & is.null(SNPlist)) stop("You requested annotation but provided no SNPlist!")
	d=arrange(dataframe, CHR, BP)
	#limit to chroms provided with function argument? if (bla bla)
	#d=d[d$CHR %in% 1:29, ] # limit to bovine autosomes?
	if ("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d) ) {
		d=na.omit(d)
		#d=dplyr::filter(d, P > 0 & P <= 1)
		if (logp) d = dplyr::mutate(d,logp = -log10(P)) else d = dplyr::mutate(d, logp = P)
		d$pos=NA
		ticks=NULL
		lastbase=0
		numchroms=nrow(dplyr::distinct(d, CHR))
		if (numchroms==1) {
			d$pos=d$BP
		} else {
			for (i in unique(d$CHR)) {
				if (i==1) {
					d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
				}   else {
					lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
					d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
				}
				ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
			}
			ticklim=c(min(d$pos),max(d$pos))
			
		}
		mycols=rep(c("gray10","gray60"),(max(d$CHR)/2)+1)
		#mycols=rep(brewer.pal(3, "Dark2"), max(numchroms))
		if (max.y=="max") maxy=ceiling(max(d$logp)) else maxy=max.y
		#if (maxy<8) maxy=8 # disable to work with F-values
		# if (annotate) d.annotate=d[as.numeric(substr(d$SNP,3,100)) %in% SNPlist, ]
		if (numchroms==1) {
			plot=qplot(pos/1000000,logp,data=d,ylab=expression(-log[10](italic(p))), 
						  xlab=paste("Chromosome",unique(d$CHR),"position", "Mbp"))
		}   else {
			plot=ggplot(d, aes(x = pos,y = logp)) + geom_point(aes(colour=factor(CHR))) + 
				ylab(expression(-log[10](italic(p))))
			plot=plot+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(d$CHR)))
			if (logp) plot=plot+scale_y_continuous(breaks=1:maxy+1, labels=1:maxy+1)
			plot=plot+scale_colour_manual(values=mycols)
		}
		#if (annotate)   plot=plot + geom_point(data=d.annotate, colour=I("green3")) 
		plot=plot + theme(legend.position = "none") 
		if (is.character(title)) plot=plot + ggtitle(title)
		plot=plot + theme(
			panel.background=element_blank(), 
			panel.grid.minor=element_blank(),
			axis.line=element_line(),
			axis.line.x=element_blank()
			
		)
		if (suggestiveline) plot=plot+geom_hline(yintercept=suggestiveline,colour="blue", alpha=I(1/3))
		if (genomewideline) plot=plot+geom_hline(yintercept=genomewideline,colour="red")
		plot
	}   else {
		stop("Make sure your data frame contains columns CHR, BP, and P")
	}
}

## QQ plot using ggplot2
gg.qq = function(pvector, title=NULL, spartan=F) {
	library(ggplot2)
	o = -log10(sort(pvector,decreasing=F))
	#e = -log10( 1:length(o)/length(o) )
	e = -log10( ppoints(length(pvector) ))
	plot=qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + stat_abline(intercept=0,slope=1, col="red")
	plot=plot+ggtitle(title)
	plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
	plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
	if (spartan) plot=plot+opts(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
	plot
}

## Regional plot to PDF. Unfinished code
plot_reg <- function(dataframe, chrom, saveplot = FALSE){
	p1 <- gg.manhattan(filter(dataframe, CHR == chrom), genomewideline = F) + 
		geom_point(aes(color = trait)) + theme(legend.position = "right") + 
		scale_colour_brewer(palette="Set1")
	if (saveplot) ggsave(filename = sprintf("QTL_plot_chrom_%s.pdf", chrom), 
								plot = p1, paper = "USr")
	p1
}
