#Title: Normalization
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-07

###Import_library###
library(data.table)

###Prepare_file_infor###
hour <- c(0,1,2,4,8,12)
time_point <- length(hour)
unit <- "hr"
normalization_factor_path <- paste('C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/genes_RefSeq_result_mRNA_normalization_factor.fpkm_table',sep="")
input_path <- paste('C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/genes_RefSeq_result_mRNA_rel.fpkm_table',sep="")
output_path <- paste('C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/genes_RefSeq_result_mRNA_rel_normalized.fpkm_table',sep="")

###Load_input_file###
data_size <- time_point+1
normalization_factor <- fread(normalization_factor_path, header=T)[,2:data_size,with=F]
nf_99 <- as.vector(as.matrix(normalization_factor[1,]))
nf_95 <- as.vector(as.matrix(normalization_factor[2,]))
input_file <- fread(input_path, header=T)

###print_header###
cat("gr_id","symbol","accession_id","locus", sep="\t", file=output_path)
cat("\t", file=output_path, append=T)
hour_label <- NULL
for(x in hour){
    hour_label <- append(hour_label, paste(x, unit, sep=""))
}
cat(hour_label, sep="\t", file=output_path, append=T)
cat("\n", sep="", file=output_path, append=T)

###calc_normalized_Relative_RPKM###
sample_size <- length(input_file[[1]])
for (x in 1:sample_size){
    exp <- as.vector(as.matrix(input_file[x,]))
}

###Extract_Stable_genes_from_dataset###
stable_genes_exp <- input_file[A >= 1 & B >= 1 & C >= 1 & D >= 1 & E >= 1 & F >= 1,] ##Should be changed!!
for (x in 1:length(stable_genes_exp[[1]])){
    data <- as.vector(as.matrix(stable_genes_exp[x,]))
    cat(data, sep="\t", file=log_path, append=T)
    cat("\n", file=log_path, append=T)
}

###Calc_normalization_factor###
data_length <- length(stable_genes_exp[1,])
quantile_99_data <- NULL
quantile_95_data <- NULL
for (x in 5:data_length){
    each_time_exp <- stable_genes_exp[[x]]
    each_time_quantile_99 <- quantile(each_time_exp, prob=0.99, na.rm=T) #99th_percentile
    each_time_quantile_95 <- quantile(each_time_exp, prob=0.95, na.rm=T) #95th_percentile
    quantile_99_data <- append(quantile_99_data, each_time_quantile_99)
    quantile_95_data <- append(quantile_95_data, each_time_quantile_95)
}
quantile_99_data <- as.vector(quantile_99_data)
quantile_95_data <- as.vector(quantile_95_data)
plot_data_quantile_99 <- data.frame(hour,quantile_99_data)
plot_data_quantile_95 <- data.frame(hour,quantile_95_data)
cat('Percentile',hour_label, sep="\t", file=output_path)
cat("\n", file=output_path, append=T)
cat('99%_percentale',quantile_99_data, sep="\t", file=output_path, append=T)
cat("\n", file=output_path, append=T)
cat('95%_percentale',quantile_95_data, sep="\t", file=output_path, append=T)
cat("\n", file=output_path, append=T)


###Plot_stable_genes_exp###
stable_genes_number <- length(stable_genes_exp[[1]])
exp_data <- as.matrix(stable_genes_exp[,5:data_length, with=F])
exp_data <- t(exp_data) #Inverse
exp_data <- factor(exp_data)
exp_data <- as.numeric(as.character(exp_data))
time_data <- as.numeric(rep(hour,stable_genes_number))
class_data <- NULL
for (x in 1:stable_genes_number){
    class_data <- append(class_data, rep(x, time_point))
}
plot_data <- data.frame(class_data,time_data,exp_data)
#plot_data <- plot_data[1:6000,] ##test

png(filename=filename,width = 1200, height = 1200)

p.scatter <- ggplot()
p.scatter <- p.scatter + layer(data=plot_data, 
                                    mapping=aes(x=time_data, y=exp_data, group=class_data), 
                                    geom="line",
                                    colour="black",
                                    size=0.02,
                                    alpha=0.05)
p.scatter <- p.scatter + layer(data=plot_data_quantile_99, 
                               mapping=aes(x=hour, y=quantile_99_data), 
                               geom="line",
                               colour="red",
                               size=0.5,
                               alpha=1)
p.scatter <- p.scatter + layer(data=plot_data_quantile_95, 
                               mapping=aes(x=hour, y=quantile_95_data), 
                               geom="line",
                               colour="blue",
                               size=0.5,
                               alpha=1)
p.scatter <- p.scatter + xlim(0,max(plot_data$time_data)) + ylim(1,max(plot_data$exp_data))
p.scatter <- p.scatter + ggtitle("Stable genes distribution")
p.scatter <- p.scatter + xlab("Relative RPKM (0hr = 1)")
p.scatter <- p.scatter + ylab("Time course(hr)")
plot(p.scatter)

dev.off() #close_fig


