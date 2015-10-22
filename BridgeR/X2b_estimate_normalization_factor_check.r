#Title: estimate_normalization_factor_check
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-07

###current_directory###
setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data")

###input_file_infor###
hour <- c(0,1,2,4,8,12)

###Estimate_normalization_factor_function###
BridgeRNormalizationFactors <- function(InputFile = "BridgeR_2_Stable_genes_dataset", group, hour, InforColumn = 4, figname = "BridgeR_2b_Normalizaion_factor_check_fig"){
    ###Import_library###
    library(data.table)
    library(ggplot2)

    ###load_input_file###
    group_number <- length(group)
    time_points <- length(hour)
    
    input_filename <- paste(InputFile, sep="")
    input_file <- fread(input_filename, header=T)
    
    input_file <- input_file[T0_1 == 1,]
    input_file <- input_file[T0_2 == 1,]
    
    house_keeping_genes_siStealth <- NULL
    house_keeping_genes_siPUM1 <- NULL
    for(x in c("GAPDH","PGK1","ENO1","ALDOA","ATP5B","PPIA")){
        if(is.null(house_keeping_genes_siPUM1)){
            house_keeping_genes_siStealth <- data.frame(hour=hour, exp=log10(as.numeric(as.vector(as.matrix(input_file[symbol == x,]))[5:10])),factor=rep(x,time_points))
            house_keeping_genes_siPUM1 <- data.frame(hour=hour, exp=log10(as.numeric(as.vector(as.matrix(input_file[symbol == x,]))[15:20])),factor=rep(x,time_points))
        }else{
            Stealth_data <- data.frame(hour=hour, exp=log10(as.numeric(as.vector(as.matrix(input_file[symbol == x,]))[5:10])),factor=rep(x,time_points))
            PUM1_data <- data.frame(hour=hour, exp=log10(as.numeric(as.vector(as.matrix(input_file[symbol == x,]))[15:20])),factor=rep(x,time_points))
            house_keeping_genes_siStealth <- rbind(house_keeping_genes_siStealth, Stealth_data)
            house_keeping_genes_siPUM1 <- rbind(house_keeping_genes_siPUM1, PUM1_data)
        }
    }
    
    all_genes_siStealth <- as.matrix(input_file[,5:10,with=F])
    all_genes_siStealth <- t(all_genes_siStealth) #Inverse
    all_genes_siStealth <- factor(all_genes_siStealth)
    all_genes_siStealth <- log10(as.numeric(as.character(all_genes_siStealth)))
    
    all_genes_siPUM1 <- as.matrix(input_file[,15:20,with=F])
    all_genes_siPUM1 <- t(all_genes_siPUM1)
    all_genes_siPUM1 <- factor(all_genes_siPUM1)
    all_genes_siPUM1 <- log10(as.numeric(as.character(all_genes_siPUM1)))
    
    hour_infor <- rep(hour, length(all_genes_siStealth)/time_points)
    class_infor <- NULL
    for (x in 1:length(input_file[[1]])){
        class_infor <- append(class_infor, rep(x, time_points))
    }
    
    all_data_siStealth <- data.frame(hour=hour, exp=all_genes_siStealth, class=class_infor)
    all_data_siPUM1 <- data.frame(hour=hour, exp=all_genes_siPUM1, class=class_infor)
    
    ###########
    png(paste(figname,"_siStealth.png",sep=""), width = 1200, height = 1200)
    
    p <- ggplot()
    p <- p + layer(data = all_data_siStealth,
                   mapping = aes(x=hour, y=exp, class=factor(class)),
                   geom = "line",
                   colour="black",
                   size=0.02,
                   alpha=0.05)
    p <- p + layer(data = house_keeping_genes_siStealth,
                   mapping = aes(x=hour, y=exp, colour=factor(factor)),
                   geom = "line",
                   size = 0.5)
    p <- p + xlim(0,12) + ylim(-2, 1)
    plot(p)
    
    dev.off() #close_fig
    plot.new()
    
    ############
    png(paste(figname,"_siPUM1.png",sep=""), width = 1200, height = 1200)
    
    p <- ggplot()
    p <- p + layer(data = all_data_siPUM1,
                   mapping = aes(x=hour, y=exp, class=factor(class)),
                   geom = "line",
                   colour="black",
                   size=0.02,
                   alpha=0.05)
    p <- p + layer(data = house_keeping_genes_siPUM1,
                   mapping = aes(x=hour, y=exp, colour=factor(factor)),
                   geom = "line",
                   size = 0.5)
    p <- p + xlim(0,12) + ylim(-2, 1)
    plot(p)
    
    dev.off() #close_fig
    plot.new()
}

###Test###
setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_Luc2_normalization")
BridgeRNormalizationFactors(InputFile = "BridgeR_3Luc2_Normalized_expression_data_siStealth_siPUM1_compatible.txt", group=c("siStealth","siPUM1"), hour=hour, figname = "BridgeR_2b_Normalizaion_factor_check_fig")

