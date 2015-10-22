#Title: BRIC-seq simulation dataset generator
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-15

###current_directory###
setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data")

###input_file_infor###
hour <- c(0,1,2,4,8,12)

###Draw_fitting_curve_function###
BridgeRSimulationGenerator <- function(test = "none", FileName = "BridgeR_0_Simulation_dataset", minHalf = 0.5, maxHalf = 48, GeneNumber = 10000, hour = c(0,1,2,4,8,12)){
    ###Import_library###
    #library(numbers)
    
    ###Parameter###
    gene_number = GeneNumber
    
    ###Make_artifical_half-life/decay_rate###
    art_half_life_original <- runif(gene_number,min=minHalf,max=maxHalf)
    art_half_life_0_24 <- art_half_life[art_half_life<=24]
    art_half_life_24_48 <- art_half_life[art_half_life>=24]
        r_36_48 <- runif(1000,min=36,max=48)
        sample_0_24 <- sample(art_half_life_0_24,size=1000)
        for(x in 1:1000){
            index_number <- which(art_half_life_0_24 == sample_0_24[x])
            art_half_life_0_24[index_number] = r_36_48[x]

                    r_0_12 <- runif(1000,min=0.5,max=12)
        sample_24_48 <- sample(art_half_life_24_48,size=1000)
        for(x in 1:1000){
            index_number <- which(art_half_life_24_48 == sample_24_48[x])
            art_half_life_24_48[index_number] = r_0_12[x]
    }
    art_half_life <- append(art_half_life_0_24,art_half_life_24_48)
    #art_half_life <- rgamma(gene_number,shape=2,scale=4)
    #art_half_life <- append(art_half_life, runif(3500,min=0.5,max=6)) #B/2000, C/3500
    #art_half_life <- append(art_half_life, runif(500,min=24,max=100)) #B/2000, C/500
    
    
    
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
    #raw_rna_expression_table <- ceiling(relative_rna_remaining_table*1000000000)
    
    #sum_rna_expression <- NULL
    #for(x in 1:length(hour)){
    #    exp <- sum(raw_rna_expression_table[,x])
    #    sum_rna_expression <- append(sum_rna_expression, exp)
    #}
    
    #LCM <- mLCM(sum_rna_expression)
    
    ###Relative_expression_dataset###
    sum_rna_expression <- NULL
    time_point <- length(hour)
    for(x in 1:time_point){
        exp <- sum(relative_rna_remaining_table[,x])
        sum_rna_expression <- append(sum_rna_expression, exp)
    }
    relative_expression_table <- NULL
    flg <- 0
    for(x in 1:time_point){
        if(flg == 0){
            relative_expression_table <- data.frame(relative_rna_remaining_table[,x])
        }else{
            relative_expression_table <- data.frame(relative_expression_table,relative_rna_remaining_table[,x]/sum_rna_expression[x]*gene_number)
        }
        flg <- 1
    }
    
    hour_label <- NULL
    for(x in hour){
        hour_label <- append(hour_label, paste("T", x, sep=""))
    }
    colnames(relative_rna_remaining_table) <- hour_label
    colnames(relative_expression_table) <- hour_label
    
    ###Write_data###
    relative_rna_remaining_table <- data.frame(relative_rna_remaining_table, half_life=art_half_life)
    filename_rpkm <- paste(FileName,"_rpkm.txt",sep="")
    filename_raw <- paste(FileName,"_raw.txt",sep="")
    
    cat("gr_id",hour_label, sep="\t",file=filename_rpkm)
    cat("\n", sep="\t",file=filename_rpkm, append=T)
    cat("gr_id",hour_label,"half_life", sep="\t",file=filename_raw)
    cat("\n", sep="\t",file=filename_raw, append=T)
    
    for(x in 1:length(relative_rna_remaining_table[,1])){
        cat(x,as.vector(as.matrix(relative_expression_table[x,])), sep="\t",file=filename_rpkm, append=T)
        cat("\n", sep="\t",file=filename_rpkm, append=T)
        cat(x,as.vector(as.matrix(relative_rna_remaining_table[x,])), sep="\t",file=filename_raw, append=T)
        cat("\n", sep="\t",file=filename_raw, append=T)
    }
    #write.table(relative_expression_table,file=filename_rpkm, sep="\t")
    #write.table(relative_rna_remaining_table,file=filename_raw, sep="\t")
}

###Test###
BridgeRSimulationGenerator(FileName = "BridgeR_0_Simulation_Test1_dataset", minHalf = 0.5, maxHalf = 48, GeneNumber = 10000, hour = c(0,1,2,4,8,12))
