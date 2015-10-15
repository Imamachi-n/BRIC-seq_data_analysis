#Title: BRIC-seq simulation dataset generator
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-15

###current_directory###
setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data")

###input_file_infor###
hour <- c(0,1,2,4,8,12)

###Draw_fitting_curve_function###
BridgeRSimulationGenerator <- function(GeneNumber = 1000, hour = c(0,1,2,4,8,12)){
    ###Import_library###
    #library(numbers)
    
    ###Parameter###
    gene_number = GeneNumber
    
    ###Make_artifical_half-life/decay_rate###
    art_half_life <- runif(gene_number,min=0.5,max=24)
    decay_model <- function(t){
        a <- log(2)/t
    }
    decay_rate <- decay_model(art_half_life) #Calc_decay_rate
    
    ###Make_relative_RNA_remaining_list###
    relative_rna_remaining_table <- NULL
    flg <- 0
    for(x in hour){
        exp_data <- exp(-decay_rate*x)
        if(flg == 0){
            relative_rna_remaining_table <- data.frame(exp_data)
        }else{
            relative_rna_remaining_table <- data.frame(relative_rna_remaining_table,exp_data)
        }
        flg <- 1
    }
    
    ###Make_raw_expression_dataset(the number of RNA molecule)###
    raw_rna_expression_table <- ceiling(relative_rna_remaining_table*1000000000)
    
    sum_rna_expression <- NULL
    for(x in 1:length(hour)){
        exp <- sum(raw_rna_expression_table[,x])
        sum_rna_expression <- append(sum_rna_expression, exp)
    }
    
    #LCM <- mLCM(sum_rna_expression)
}

###Test###
BridgeRSimulationGenerator()
