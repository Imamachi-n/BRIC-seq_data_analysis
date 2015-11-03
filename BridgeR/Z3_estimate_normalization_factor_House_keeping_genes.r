#Title: estimate_normalization_factor
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-07

setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_siStealth_siPUM1_ver1/time_course_0_1_2_4_8_12h")

InputFile <- "BridgeR_1_Relative_expression_data_cuffnorm_total_siStealth_siPUM1.txt"

hour <- c(0,1,2,4,8,12)
group <- c("A_siStealth","B_siPUM1")


###Estimate_normalization_factor_function###
BridgeRNormalizationFactorsHK <- function(InputFile, group, hour, InforColumn = 4, InforHKGenes = 2, HKGenes = c("GAPDH","PGK1","PPIA","ENO1","ATP5B","ALDOA"), OutputFile = "BridgeR_3B_Normalizaion_factor_from_House_keeping_genes.txt"){
    ###Import_library###
    library(data.table)

    ###Calc_normalization_factor###
    group_number <- length(group)
    time_points <- length(hour)
    input_file <- fread(InputFile, header=T)
    
    ###Output_file_infor###nfname
    hour_label <- NULL
    for(x in hour){
        label <- x
        if(x < 10){
            label <- paste("0",x,sep="")
        }
        hour_label <- append(hour_label, paste("T", label, "_", a, sep=""))
    }
    cat('Sample',hour_label, sep="\t", file=OutputFile)
    cat("\n", file=OutputFile, append=T)
    
    for(a in 1:group_number){
        ###Search_House-keeping_genes###
        HKgenes_infor <- input_file[[InforHKGenes]]
        HKgenes_infor_index <- NULL
        for(x in 1:length(HKGenes)){
            HKgenes_infor_index <- append(HKgenes_infor_index, which(HKgenes_infor == HKGenes[x]))
        }
        
        ###Information&exp_data column###
        infor_st <- 1 + (a - 1)*(time_points + InforColumn)
        infor_ed <- (InforColumn)*a + (a - 1)*time_points
        exp_st <- infor_ed + 1
        exp_ed <- infor_ed + time_points

        HKGenes_raw_data <- NULL
        for(x in 1:length(HKgenes_infor_index)){
            if(x == 1){
                HKGenes_raw_data <- input_file[HKgenes_infor_index[x],exp_st:exp_ed,with=F]
            }else{
                HKGenes_raw_data <- rbind(HKGenes_raw_data,input_file[HKgenes_infor_index[x],exp_st:exp_ed,with=F])
            }
        }
        
        nf <- NULL
        for(x in 1:length(HKGenes_raw_data)){
            if(x == 1){
                nf <- as.vector(as.matrix(HKGenes_raw_data[x,]))
            }else{
                nf <- HKGenes_raw_data[x,] * nf
            }
        }
        nf <- nf^(1/length(HKGenes))
        
        cat(group[a],nf, sep="\t", file=OutputFile, append=T)
        cat("\n", file=OutputFile, append=T)

    }
}
