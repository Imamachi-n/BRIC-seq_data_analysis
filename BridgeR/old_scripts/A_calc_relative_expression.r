#Title: Calc_relative_expression
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-07

###Prepare_file_infor###
hour <- c(0,1,2,4,8,12)
time_point <- length(hour)
unit <- "hr"
cutoff <- 0.1
input_path <- paste('C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/genes_RefSeq_result_mRNA.fpkm_table',sep="")
output_path <- paste('C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/genes_RefSeq_result_mRNA_rel.fpkm_table',sep="")

###Load_input_file###
input_file = read.table(input_path, header=T)

###print_header###
cat("gr_id","symbol","accession_id","locus", sep="\t", file=output_path)
cat("\t", file=output_path, append=T)
hour_label <- NULL
for(x in hour){
    hour_label <- append(hour_label, paste(x, unit, sep=""))
}
cat(hour_label, sep="\t", file=output_path, append=T)
cat("\n", sep="", file=output_path, append=T)

###Read_each_line###
for (x in 1:length(input_file[,1])){
    data <- input_file[x,]
    data_length <- length(data)
    gene_infor <- as.vector(as.matrix(data[1:4]))
    cat(gene_infor, sep="\t", file=output_path, append=T)
    cat("\t", sep="", file=output_path, append=T)
    exp <- data[5:data_length]
    exp <- as.numeric(exp)
    start_time <- exp[1]
    if(start_time <= cutoff){
        cat(rep(0,time_point), sep="\t", file=output_path, append=T)
        cat("\n", sep="", file=output_path, append=T)
        next
    }
    rel_exp <- NULL
    for (y in 1:time_point){
        rel_exp <- append(rel_exp, exp[y]/start_time)
    }
    cat(rel_exp, sep="\t", file=output_path, append=T)
    cat("\n", sep="", file=output_path, append=T)
}

