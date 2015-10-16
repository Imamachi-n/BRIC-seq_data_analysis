#Title: Normalization
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-08

###current_directory###
setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data")

###input_file_infor###
hour <- c(0,1,2,4,8,12)
#group <- c("siCTRL","siStealth")
group <- c("siStealth","siPUM1")
#files <- c("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/siCTRL_genes_RefSeq_result_mRNA.fpkm_table",
#           "C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/siStealth_genes_RefSeq_result_mRNA.fpkm_table")
files <- c("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/siStealth_genes_RefSeq_result_mRNA.fpkm_table",
           "C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/siPUM1_genes_RefSeq_result_mRNA.fpkm_table")

###Estimate_normalization_factor_function###
BridgeRNormalization <- function(filename = "BridgeR_1_Relative_expression_data.txt", group, hour, InforColumn = 4, SelectNormFactor=T, NormFactor = "BridgeR_2_Normalizaion_factor_dataset", OutputFile = "BridgeR_3_Normalized_expression_data.txt"){
    ###Import_library###
    library(data.table)
    
    ###Load_NormFactor_file###
    group_number <- length(group)
    time_points <- length(hour)
    nf_st <- 2
    nf_ed <- time_points + 1
    nf_99 <- NULL
    nf_95 <- NULL
    for(a in 1:group_number){
        NormFactor_name <- paste(NormFactor,a,"_",group[a],".txt", sep="")
        normalization_factor <- fread(NormFactor_name, header=T)[,nf_st:nf_ed,with=F]
        if(is.null(nf_99) && is.null(nf_95)){
            nf_99 <- as.vector(as.matrix(normalization_factor[1,]))
            nf_95 <- as.vector(as.matrix(normalization_factor[2,]))
        }else{
            nf_99 <- append(nf_99, as.vector(as.matrix(normalization_factor[1,])))
            nf_95 <- append(nf_95, as.vector(as.matrix(normalization_factor[2,])))
        }
    }
    
    ###Load_files###
    input_file <- fread(filename, header=T)
    output_file <- OutputFile
    
    ###print_header###
    cat("",file=output_file)
    hour_label <- NULL
    for(a in 1:group_number){
        if(!is.null(hour_label)){
            cat("\t", file=output_file, append=T)
        }
        hour_label <- NULL
        for(x in hour){
            hour_label <- append(hour_label, paste("T", x, "_", a, sep=""))
        }
        infor_st <- 1 + (a - 1)*(time_points + InforColumn)
        infor_ed <- (InforColumn)*a + (a - 1)*time_points
        infor <- colnames(input_file)[infor_st:infor_ed]
        cat(infor,hour_label, sep="\t", file=output_file, append=T)
    }
    cat("\n", sep="", file=output_file, append=T)
    
    ###calc_normalized_Relative_RPKM###
    gene_number <- length(input_file[[1]])
    for(x in 1:gene_number){
        data <- as.vector(as.matrix(input_file[x,]))
        for(a in 1:group_number){
            if(a != 1){
                cat("\t", sep="", file=output_file, append=T)
            }
            infor_st <- 1 + (a - 1)*(time_points + InforColumn)
            infor_ed <- (InforColumn)*a + (a - 1)*time_points
            exp_st <- infor_ed + 1
            exp_ed <- infor_ed + time_points
            gene_infor <- data[infor_st:infor_ed]
            exp <- as.numeric(data[exp_st:exp_ed])
            
            nf_99_st <- (time_points)*(a - 1) + 1
            nf_99_ed <- (time_points)*a
            nf_99_exp <- nf_99[nf_99_st:nf_99_ed]
            nf_95_st <- (time_points)*(a - 1) + 1
            nf_95_ed <- (time_points)*a
            nf_95_exp <- nf_95[nf_95_st:nf_95_ed]
            
            normalized_exp <- NULL
            if(SelectNormFactor == TRUE){
                if(exp[2] >= nf_95_exp[2] & exp[3] >= nf_95_exp[3]){
                    normalized_exp <- exp/nf_99_exp
                }else{
                    normalized_exp <- exp/nf_95_exp
                }
            }else{
                normalized_exp <- exp/nf_95_exp
            }
            cat(gene_infor, sep="\t", file=output_file, append=T)
            cat("\t", sep="\t", file=output_file, append=T)
            cat(normalized_exp, sep="\t", file=output_file, append=T)
        }
        cat("\n", sep="\t", file=output_file, append=T)
    }
}

###Test###
#BridgeRNormalization(group=group, hour=hour)
#BridgeRNormalization(filename = "BridgeR_1_Relative_expression_data_siStealth_siPUM1.txt", group=group, hour=hour, OutputFile = "BridgeR_3_Normalized_expression_data_siStealth_siPUM1.txt")

#BridgeRNormalization(filename = "BridgeR_0_Simulation_dataset_rpkm.txt", group=c("Simulation_A"), hour=hour, InforColumn = 1, OutputFile = "BridgeR_3_Normalized_simulation_data.txt")
#BridgeRNormalization(filename = "BridgeR_0_Simulation_B_dataset_rpkm.txt", group=c("Simulation_B"), hour=hour, InforColumn = 1, OutputFile = "BridgeR_3_Normalized_simulation_B_data.txt")
BridgeRNormalization(filename = "BridgeR_0_Simulation_C_dataset_rpkm.txt", group=c("Simulation_C"), hour=hour, InforColumn = 1, OutputFile = "BridgeR_3_Normalized_simulation_C_data.txt")

