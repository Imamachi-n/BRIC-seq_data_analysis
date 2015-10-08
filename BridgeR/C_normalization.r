#Title: Normalization
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-08

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
nf_data_size <- time_point+1
normalization_factor <- fread(normalization_factor_path, header=T)[,2:nf_data_size,with=F]
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
data_size <- length(input_file[1,])
for(x in 1:sample_size){
    data <- as.vector(as.matrix(input_file[x,]))
    gene_infor <- data[1:4]
    exp <- as.numeric(data[5:data_size])
    normalized_exp <- NULL
    if(exp[2] >= nf_95[2] & exp[3] >= nf_99[3]){ #should be changed!!
        normalized_exp <- exp/nf_99
    }else{
        normalized_exp <- exp/nf_95
    }
    cat(gene_infor, sep="\t", file=output_path, append=T)
    cat("\t", sep="\t", file=output_path, append=T)
    cat(normalized_exp, sep="\t", file=output_path, append=T)
    cat("\n", sep="\t", file=output_path, append=T)
}

