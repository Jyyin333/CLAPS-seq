#!/usr/bin/env Rscript

# import request libraries 
library(ggplot2)
library(argparse)
library(stringr)

parser <- ArgumentParser(description=
	'Plot OG profile in promoter with or without G4 and test statistical signifiance for their differences')

parser$add_argument("--treat", required=T, nargs=2, help="2 Treat sample matrices, [G4+ G4-]")
parser$add_argument("--control", required=T, nargs=2, help="2 Input sample matrices, [G4+ G4-]")
parser$add_argument("--outfig", required=T, help="file name to save image.")
parser$add_argument("--interval", required=F, type="integer", default=1000, help="specify an interval around TSS for comparison, default:1000.")

args <- parser$parse_args()


# parse header line
process_header <- function(infile){
  fl <- gzfile(infile, open='r')
  header <- readLines(fl, n=1)
  close(fl)
  # check if header line exists
  if (!startsWith(header,"#")){
    stop("Do not detect the header line, please check input file\n
         note:the header line should starts with '#'.\n")
  }
  
  # match the header line
  upstream <- as.integer(str_split(
    str_extract(header, pattern = "Upstream:\\d{1,}"), ':'
  )[[1]][2]
  )
  
  downstream <- as.integer(str_split(
    str_extract(header, pattern = "Downstream:\\d{1,}"), ':'
  )[[1]][2]
  )
  
  binsize <- as.integer(str_split(
    str_extract(header, pattern = "BinSize:\\d{1,}"), ':'
  )[[1]][2]
  )
  
  if(is.na(upstream) | is.na(downstream) | is.na(binsize)){
    stop('Do not detect necessary info in header line.\n\t
         [Upstream, Downstream, BinSize] should be included in Header.\n')
  }
  
  res <- c(upstream, downstream, binsize)
  return (res)
}


cal_means <- function(df, info, region=1000){
  up <- info[1]
  down <- info[2]
  bs <- info[3]

  if (up<region/2 | down<region/2){
  	stop('Error:\n\tThe interval length out of range.\n')
  }

  sub_bins <- region / bs / 2
  nbins <- ncol(df) / 2
  center <- nbins / 2
  
  sub_ts_df <- df[,seq(center-sub_bins-1, center+sub_bins)]
  sub_nts_df <- df[,seq(nbins+center-sub_bins-1, nbins+center+sub_bins)]
  ts_res <- rowMeans(sub_ts_df)
  nts_res <- rowMeans(sub_nts_df)
  res <- as.matrix(data.frame(ts=ts_res,nts=nts_res))
  
  return (res)
}



# name params
treat_a <- args$treat[1]
treat_b <- args$treat[2]
ctrl_a <- args$control[1]
ctrl_b <- args$control[2]
outfig <- args$outfig
interval <- args$interval


# check treat/control header
t_a_header <- process_header(treat_a)
t_b_header <- process_header(treat_b)
c_a_headr <- process_header(ctrl_a)
c_b_headr <- process_header(ctrl_b)

# if not same in header than echo warning message
flag <- assertthat::are_equal(t_a_header,t_b_header) & assertthat::are_equal(t_a_header,c_a_headr) & assertthat::are_equal(t_a_header,c_b_headr)
if (flag == FALSE){
	warning("Input files are not equal in Header infos.\n")
}


# load control samples and calculate mean RPKM
input_og4 <- read.table(ctrl_a,
                        header=F, comment.char = '#', sep='\t')

input_wog4 <- read.table(ctrl_b,
                        header=F, comment.char = '#', sep='\t')

# default 1kb around TSS ( megre TS & NTS)
inp_og4_meanRPKM <- mean(rowSums(cal_means(input_og4, c_a_headr, interval)))
inp_wog4_meanRPKM <- mean(rowSums(cal_means(input_wog4, c_b_headr, interval)))



# load treat samples
og4 <- read.table(treat_a,
                  header=F, comment.char = '#', sep='\t')

wog4 <- read.table(treat_b,
                  header=F, comment.char = '#', sep='\t')

# normalize to input
og4_1k_norm_values <- rowSums(cal_means(og4, t_a_header, interval)) / inp_og4_meanRPKM
wog4_1k_norm_values <- rowSums(cal_means(wog4, t_b_header, interval)) / inp_wog4_meanRPKM


# wilcoxon test
w_p <- wilcox.test(og4_1k_norm_values,
                   wog4_1k_norm_values,
                   alternative = "less",
                   conf.level = 0.95)$p.value

# print P-value into stdout
sample_name <- strsplit(basename(treat_a),'_',fixed=TRUE)[[1]][1]
cat(sample_name,'\n')
cat('Wilcox test P-value:\t', w_test_p, '\n')


# prepare data for plotting
df <- data.frame(counts=c(wog4_1k_norm_values,og4_1k_norm_values),
                 sample=c(rep('G4-',length(wog4_1k_norm_values)),
                          rep('G4+',length(og4_1k_norm_values)))
)

#adjust axis-x orders
df$sample <-factor(df$sample,levels=c('G4+','G4-'))
ylim <- boxplot.stats(df$counts)$stats[c(1, 5)]

#pdf(args$out,width=6,height=4.5)
ggplot(data=df) + geom_boxplot(aes(x=sample,y=counts,fill=sample),outlier.shape = NA) +
  labs(x='',title='',y='Normalized OG signal') + theme(plot.title=element_text(size=20,face="bold",hjust=0.5)) +
  theme(axis.title.y=element_text(size=15,face="bold",vjust=0.5))  +
  theme(axis.text.x=element_text(size=12,vjust=0.5,face='italic')) +
  theme(legend.title = element_blank(), legend.text = element_text(size=12)) + coord_cartesian(ylim = ylim*1.5) +
  ggsave(outfig,width=8,height=5.5,dpi=300)

dev.off()

cat('Done.\n')
