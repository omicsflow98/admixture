#!/usr/bin/env Rscript

correlate_components <- function(k, k_min) {
  start_this_k = 1+sum(1:k-1)-(sum(1:k_min-1))+1
  end_this_k = start_this_k+k-1
  end_prev_k = start_this_k-1
  start_prev_k = end_prev_k-(k-2)
  # print (paste0("This K: ",start_this_k,":",end_this_k))
  # print (paste0("Prev K: ",start_prev_k,":",end_prev_k))
  cor(raw_data[start_prev_k:end_prev_k], raw_data[start_this_k:end_this_k])
}

## A function that returns the column name of the best correlated component in 
##   the K-1 run, for each component in the K run.
fix_colours <- function(k, k_min) {
  ## If the K being processed is the minimum K in the data, keep columns 
  ##  unchanged (since there is nothing to compare to).
  if (k == k_min) {
    return(paste0(k, ":", 1:k))
  }
  
  ## Calculate a correlation matrix for each component of this K to the 
  ##   components of K-1.
  cor_mat <- correlate_components(k, k_min)
  
  ## Find the most correlated component from the last K for each component in 
  ##   this K.
  component_order <- c()
  for (x in 1:k-1) {
    component_order <- append(component_order, which.max(cor_mat[x, ]))
  }
  
  ## If one component in the last K is the most correlated with two components 
  ##   on this K, find the next best correlated component, and if that is 
  ##   unique, assign that as the correct component.
  ## If it is not unique, repeat the process until a unique component is found.
  if (any(duplicated(component_order)) == T && 
      sum(duplicated(component_order)) == 1) {
    duplicate <- which(duplicated(component_order))
    condition <- T
    top_correlates <- c()
    while (condition == T) {
      ## 
      ## until one that isn't already crrelated with another component is found.
      top_correlates <- c(top_correlates, which.max(cor_mat[duplicate, ]))
      cor_mat[duplicate, top_correlates] <- NA
      component_order[duplicate] <- which.max(cor_mat[duplicate, ])
      if (any(duplicated(component_order)) == F) {condition = F}
    }
  } else if (sum(duplicated(component_order)) > 1) {
    stop (paste0("Correlation of components failed. Usually this is caused by high CV errors for some of the components you are trying to plot. 
Please consider limiting your input dataset to K=",k_min," to ",k-1,".

You can use this command to extract the suggested columns from the input file:
    cut -d ' ' -f 1-",sum(seq(k_min,k-1))+2,"
    "), call.=FALSE)
  }
  ## If a component hasn't been resolved yet, add it as the newest component.
  missing_component = setdiff(1:k, component_order)
  component_order<-append(component_order, missing_component)
  return(paste0(k,":", component_order))
}


suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

option_list = list(
  make_option(c("-Q", "--Q_files"), default=NULL,
              help="properly formatted file with Q values", metavar = "Q value file"),
  make_option(c("-V", "--Eigen"), default=NULL,
              help="directory that contains your eigen files for PCA", metavar = "Eigen files"),
  make_option(c("-E", "--cv_file"), default=NULL,
              help="file containing your CV error values", metavar = "CV file")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

if(is.null(opt$Q_files) | is.null(opt$Eigen) | is.null(opt$cv_file)) {
  print_help(opt_parser)
  stop("Please provide the necessary files ", call.=FALSE)
  
}

cv_plot <- read_delim(opt$cv_file, delim=" ", col_names=c("K", "CV_error"))

ggplot(cv_plot, aes(x=K, y=CV_error)) +
  geom_point() +
  geom_line() +
  theme_light()
ggsave("cv_plot.png")

## read data
raw_data <- read_delim(opt$Q_files, " ", col_types = cols())

## Sort components of each K according to correlation with components of K-1. 
##   This needs to happen per K, otherwise the correlations will not match 
##   beyond the first pair of Ks.
header <- names(raw_data) ## Take column names from original data

## 'Ind' and 'Pop' should always be at the start of the reformatted data.
refcols <- c("Samples")

othercols <- header[header != refcols]

split_cols <- do.call(rbind, strsplit(othercols, ":", fixed = TRUE))
x <- as.numeric(split_cols[,1])
y <- as.numeric(split_cols[,2])

order_idx <- order(x, y)

raw_data <- raw_data[, c(refcols, othercols[order_idx])]

header <- names(raw_data) ## Take column names from original data

## Infer min and max K values.
k_min <- as.numeric(str_split_fixed(names(raw_data[,2]),":",2)[1])
k_max <- as.numeric(str_split_fixed(names(raw_data[,ncol(raw_data)]),":",2)[1])

## For each K in the data, use fix_colours to extract the vector of most correlated
##   column names for each Component in the K. Then sort the components of this 
##   K in the raw data.
for (k in k_min:k_max) { 
  refcols <- c(refcols, fix_colours(k,k_min)) ##    
  raw_data <- raw_data[, c(refcols, setdiff(names(raw_data), refcols))] ## 
}

## Finally, fix the column names so that inference of component numbers is correct
names(raw_data) <- header 

## Flatten data to long format
long_data <- gather(raw_data, temp, value, 2:ncol(raw_data)) %>%
  separate_wider_delim(cols = temp, delim = ":", names = c("k", "Q"))

long_data$k <- as.integer(long_data$k)

long_data %>%
  ggplot(.,aes(x=Samples,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set3",name="K", labels=seq(1:max(cv_plot$K))) +
  facet_wrap(~k,ncol=1)
ggsave("admixture.png")

k_min <- cv_plot$K[which(cv_plot$CV_error == min(cv_plot$CV_error))]

long_data %>%
  filter(k == k_min) %>%
  ggplot(.,aes(x=Samples,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set3",name="K", labels=1:k_min)
ggsave("admixture_bestk.png")

pca <- read_table(paste0(basename(opt$Eigen),".eigenvec"), col_names = FALSE)
pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

eigenval <- scan(paste0(basename(opt$Eigen),".eigenval"))

pve <- data.frame(PC=1:length(eigenval), pve = eigenval/sum(eigenval)*100)

ggplot(pve, aes(PC, pve)) +
  geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_light()
ggsave("variance.png")

ggplot(pca, aes(PC1, PC2)) +
  geom_point(size = 3) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("pca_plot.png")