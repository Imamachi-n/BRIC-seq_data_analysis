#Title: estimate_normalization_factor
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-07

###current_directory###
setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data")

###input_file_infor###
hour <- c(0,1,2,4,8,12)
group <- c("siCTRL","siStealth")
files <- "C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_1_Relative_expression_data.txt"

###Estimate_normalization_factor_function###
BridgeRStableGenesDataSet <- function(filename, group, hour, InforColumn = 4, OutputFile = "BridgeR_2_Stable_genes_dataset"){
    ###Import_library###
    library(data.table)

    ###Prepare_files###
    time_points <- length(hour)
    input_file <- suppressWarnings(fread(filename, header=T))
    
    ###Extract_Stable_genes_from_dataset###
    group_number <- length(group)
    for(a in 1:group_number){
        output_filename <- paste(OutputFile,a,"_",group[a],".txt",sep="")
        
        ###print_header###
        cat("",file=output_filename)
        hour_label <- NULL
        for(x in hour){
            hour_label <- append(hour_label, paste("T", x, "_", a, sep=""))
        }
        infor_st <- 1 + (a - 1)*(time_points + InforColumn)
        infor_ed <- (InforColumn)*a + (a - 1)*time_points
        infor <- colnames(input_file)[infor_st:infor_ed]
        cat(infor,hour_label, sep="\t", file=output_filename, append=T)
        cat("\n", sep="", file=output_filename, append=T) 
        
        ###Extract_Stable_genes_from_dataset###
        gene_number <- length(input_file[[1]]) #Total number of genes
        for(y in 1:gene_number){
            data <- as.vector(as.matrix(input_file[y,]))
            ###Exp_data###
            exp_st <- infor_ed + 1
            exp_ed <- infor_ed + time_points
            exp <- data[exp_st:exp_ed]
            exp <- as.numeric(exp)
            if(all(exp >= 1)){
                ###Infor_data###
                gene_infor <- data[infor_st:infor_ed]
                cat(gene_infor, sep="\t", file=output_filename, append=T)
                cat("\t", sep="", file=output_filename, append=T)
                ###Exp_data###
                cat(exp, sep="\t", file=output_filename, append=T)
                cat("\n", sep="", file=output_filename, append=T)
            }
        }
    }
    
    
}

BridgeRNormalizationFactors <- function(InputFile = "BridgeR_2_Stable_genes_dataset", group, hour, InforColumn = 4, NormFactor = "BridgeR_2_Normalizaion_factor_dataset", figname = "BridgeR_2_Normalizaion_factor_fig_dataset"){
    ###Import_library###
    library(data.table)
    library(ggplot2)

    ###Calc_normalization_factor###
    group_number <- length(group)
    time_points <- length(hour)
    
    for(a in 1:group_number){
        ###Load_files###
        input_filename <- paste(InputFile,a,"_",group[a],".txt", sep="")
        input_file <- fread(input_filename, header=T)
        
        ###Output_files###
        output_filename <- paste(NormFactor,a,"_",group[a],".txt", sep="")
        
        hour_label <- NULL
        for(x in hour){
            hour_label <- append(hour_label, paste("T", x, "_", a, sep=""))
        }
        
        quantile_99_data <- NULL
        quantile_95_data <- NULL
        infor_st <- 1
        infor_ed <- InforColumn
        exp_st <- infor_ed + 1
        exp_ed <- infor_ed + time_points
        for(x in exp_st:exp_ed){
            each_time_exp <- input_file[[x]]
            each_time_quantile_99 <- quantile(each_time_exp, prob=0.99, na.rm=T) #99th_percentile
            each_time_quantile_95 <- quantile(each_time_exp, prob=0.95, na.rm=T) #95th_percentile
            quantile_99_data <- append(quantile_99_data, each_time_quantile_99)
            quantile_95_data <- append(quantile_95_data, each_time_quantile_95)
        }
        quantile_99_data <- as.vector(quantile_99_data)
        quantile_95_data <- as.vector(quantile_95_data)
        plot_data_quantile_99 <- data.frame(hour,quantile_99_data)
        plot_data_quantile_95 <- data.frame(hour,quantile_95_data)
        cat('Percentile',hour_label, sep="\t", file=output_filename)
        cat("\n", file=output_filename, append=T)
        cat('99%_percentale',quantile_99_data, sep="\t", file=output_filename, append=T)
        cat("\n", file=output_filename, append=T)
        cat('95%_percentale',quantile_95_data, sep="\t", file=output_filename, append=T)
        cat("\n", file=output_filename, append=T)
        
        
        ###Plot_stable_genes_exp###
        figfile <- paste(figname,a,"_",group[a],".png", sep="")
        
        stable_genes_number <- length(input_file[[1]])
        exp_data <- as.matrix(input_file[,exp_st:exp_ed, with=F])
        exp_data <- t(exp_data) #Inverse
        exp_data <- factor(exp_data)
        exp_data <- as.numeric(as.character(exp_data))
        time_data <- as.numeric(rep(hour,stable_genes_number))
        class_data <- NULL
        for (x in 1:stable_genes_number){
            class_data <- append(class_data, rep(x, time_points))
        }
        plot_data <- data.frame(class_data,time_data,exp_data)
        #plot_data <- plot_data[1:6000,] ##test
        
        png(filename=figfile,width = 1200, height = 1200)
        
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
        p.scatter <- p.scatter + xlab("Relative RPKM (Time0 = 1)")
        p.scatter <- p.scatter + ylab("Time course")
        plot(p.scatter)
        
        dev.off() #close_fig
        plot.new()
    }
}

###Test###
BridgeRStableGenesDataSet(filename=files, group=group, hour=hour)
BridgeRNormalizationFactors(group=group, hour=hour)
