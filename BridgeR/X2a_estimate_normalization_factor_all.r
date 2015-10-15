#Title: estimate_normalization_factor
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-07

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
BridgeRNormalizationFactorsAll <- function(InputFile = "BridgeR_2_Stable_genes_dataset", group, hour, InforColumn = 4, YMin = -2, YMax = 2, figname = "BridgeR_2_Normalizaion_factor_fig_all_dataset"){
    ###Import_library###
    library(data.table)
    library(ggplot2)

    ###Calc_normalization_factor###
    group_number <- length(group)
    time_points <- length(hour)
    
    for(a in 1:group_number){
        ###Load_files###
        #input_filename <- paste(InputFile,a,"_",group[a],".txt", sep="")
        input_filename <- InputFile
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
        quantile_99_data <- log10(as.vector(quantile_99_data)) #
        quantile_95_data <- log10(as.vector(quantile_95_data)) #
        plot_data_quantile_99 <- data.frame(hour,quantile_99_data)
        plot_data_quantile_95 <- data.frame(hour,quantile_95_data)
        
        ###Plot_stable_genes_exp###
        figfile <- paste(figname,a,"_",group[a],".png", sep="")
        
        stable_genes_number <- length(input_file[[1]])
        exp_data <- as.matrix(input_file[,exp_st:exp_ed, with=F])
        exp_data <- t(exp_data) #Inverse
        exp_data <- factor(exp_data)
        exp_data <- as.numeric(as.character(exp_data))
        exp_data <- log10(exp_data) #
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
        p.scatter <- p.scatter + xlim(0,max(plot_data$time_data)) + ylim(YMin, YMax)
        p.scatter <- p.scatter + ggtitle("Stable genes distribution")
        p.scatter <- p.scatter + xlab("Time course")
        p.scatter <- p.scatter + ylab("Relative RPKM (Time0 = 1)")
        plot(p.scatter)
        
        dev.off() #close_fig
        plot.new()
    }
}

###Test###
#BridgeRNormalizationFactors(group=group, hour=hour)
#BridgeRNormalizationFactors(InputFile = "BridgeR_2_All_genes_dataset",group=group, hour=hour, figname = "BridgeR_2_Normalizaion_factor_fig_dataset_test")
#BridgeRNormalizationFactors(InputFile = "siStealth_genes_RefSeq_result_mRNA.fpkm_table",group=group, hour=hour, figname = "BridgeR_2_Normalizaion_factor_fig_dataset_test_siStealth")

BridgeRNormalizationFactorsAll(InputFile = "BridgeR_0_Simulation_dataset_rpkm.txt", YMin=-7, YMax = 0.5, group=c("Simulation_A"), hour=hour, InforColumn = 1, figname = "BridgeR_2_Normalizaion_factor_fig_all_dataset")

