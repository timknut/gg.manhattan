gg.manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP",
							 cols=c("gray10", "gray60"), plot_title=NULL, chrlabs=NULL,
							 suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8),
							 logp=TRUE) {

## base plot function. Must test and add several features
# manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP",
# 							 cols=c("gray10", "gray60"), chrlabs=NULL, title= NULL,
# 							 suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8),
# 							 highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {


	# Require ggplot2, install if not present in labrary.
	if (!require(ggplot2)) install.packages("ggpplot2")

	# Not sure why, but package check will warn without this.
	CHR=BP=P=index=NULL

	# Check for sensible dataset
	## Make sure you have chr, bp and p columns.
	if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
	if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
	if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
	## warn if you don't have a snp column
	if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
	## make sure chr, bp, and p columns are numeric.
	if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
	if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
	if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))

	# Create a new data.frame with columns called CHR, BP, and P.
	d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])

	# If the input data frame has a SNP column, add it to the new data frame you're creating.
#	if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])

	# Set positions, ticks, and labels for plotting
	## Sort and keep only values where is numeric.
	#d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
	d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
	d <- d[order(d$CHR, d$BP), ]
	#d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
	if (logp) {
		d$logp <- -log10(d$P)
	} else {
		d$logp <- d$P
	}
	d$pos=NA


	# Fixes the bug where one chromosome is missing by adding a sequential index column.
	d$index=NA
	ind = 0
	for (i in unique(d$CHR)){
		ind = ind + 1
		d[d$CHR==i,]$index = ind
	}

	# This section sets up positions and ticks. Ticks should be placed in the
	# middle of a chromosome. The a new pos column is added that keeps a running
	# sum of the positions of each successive chromsome. For example:
	# chr bp pos
	# 1   1  1
	# 1   2  2
	# 2   1  3
	# 2   2  4
	# 3   1  5
	nchr = length(unique(d$CHR))
	if (nchr==1) { ## For a single chromosome
		## Uncomment the next line to plot single chr results in Mb
		d$pos=d$BP/1e6
		ticks=floor(length(d$pos))/2+1
		xlabel = paste('Chromosome',unique(d$CHR),'position(Mb)')
		labs = ticks
	} else { ## For multiple chromosomes
		lastbase=0
		ticks=NULL
		for (i in unique(d$index)) {
			if (i==1) {
				d[d$index==i, ]$pos=d[d$index==i, ]$BP
			} else {
				lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
				d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
			}
			# Old way: assumes SNPs evenly distributed
			# ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
			# New way: doesn't make that assumption
			ticks = c(ticks, (min(d[d$CHR == i,]$pos) + max(d[d$CHR == i,]$pos))/2 + 1)
		}
		xlabel = 'Chromosome'
		## Set x labels to Chromsomes, or provided custom list from chrlabs argument.
		if (!is.null(chrlabs)) {
			labs <- chrlabs
		} else {
			labs <- unique(d$CHR)
			}
	}

	# Initialize plot
	xmax = ceiling(max(d$pos) * 1.03)
	xmin = floor(max(d$pos) * -0.03)
	ymin = floor(min(d$logp))
	ymax = ceiling(max(d$logp))
	# Set ymax to even number.
	if ((ymax %% 2) != 0 ) {
		ymax = ymax + 1
	}
	mycols = rep(cols, nchr/2+1)
		if (nchr==1) {
			plot=qplot(pos,logp,data=d,ylab=expression(-log[10](italic(p))),
						  xlab=xlabel) +
				scale_y_continuous(breaks=seq(2,ymax,2), labels=seq(2,ymax,2),
										 limits = c(ymin-0.5, ymax), expand = c(0,0))
		}   else {
			plot=ggplot(d, aes(x = pos,y = logp))
			plot=plot + geom_point(aes(colour=factor(CHR))) + ylab(expression(-log[10](italic(p))))
			plot=plot+scale_x_continuous(name=xlabel, breaks=ticks, labels=labs)
			if (logp) plot=plot+scale_y_continuous(breaks=seq(2,ymax,2), labels=seq(2,ymax,2),
																limits = c(ymin-0.5, ymax), expand=c(0,0)) # expand ensures pretty y-axis
			plot=plot+scale_colour_manual(values=mycols)
		}
		#if (annotate)   plot=plot + geom_point(data=d.annotate, colour=I("green3"))
		if (!is.null(plot_title)) plot=plot + ggtitle(plot_title)
		plot=plot + theme(
			panel.background=element_blank(),
			panel.grid.minor=element_blank(),
			axis.line=element_line(),
			axis.line.x=element_blank(),
			legend.position = "none"
		)
		if (suggestiveline) plot=plot+geom_hline(yintercept=suggestiveline,colour="blue", alpha=I(1/3))
		if (genomewideline) plot=plot+geom_hline(yintercept=genomewideline,colour="red")
		plot
	}
