#Title: Calc_relative_expression
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-23

#InputFile <- "BridgeR_1_Relative_expression_data_siCTRL_siPUM2_compatible.txt"
#group <- c("siCTRL","siPUM2")
#hour <- c(0,1,2,4,8,12)
#OutputFile <- "test.txt"

###Check_BRIC-seq_dataset###
test_q <- function(x,y){
    q_99 <- as.vector(quantile(x, prob=0.99, na.rm=T))
    q_95 <- as.vector(quantile(x, prob=0.95, na.rm=T))
    q_90 <- as.vector(quantile(x, prob=0.90, na.rm=T))
    q_80 <- as.vector(quantile(x, prob=0.80, na.rm=T))
    q_70 <- as.vector(quantile(x, prob=0.70, na.rm=T))
    q_60 <- as.vector(quantile(x, prob=0.60, na.rm=T))
    q_50 <- as.vector(quantile(x, prob=0.50, na.rm=T))
    q_40 <- as.vector(quantile(x, prob=0.40, na.rm=T))
    q_30 <- as.vector(quantile(x, prob=0.30, na.rm=T))
    q_20 <- as.vector(quantile(x, prob=0.20, na.rm=T))
    q_10 <- as.vector(quantile(x, prob=0.10, na.rm=T))
    q_5 <- as.vector(quantile(x, prob=0.05, na.rm=T))
    q_1 <- as.vector(quantile(x, prob=0.01, na.rm=T))
    vec <- c(q_99,q_95,q_90,q_80,q_70,q_60,q_50,q_40,q_30,q_20,q_10,q_5,q_1)
    factor_label <- c("99%","95%","90%","80%","70%","60%","50%","40%","30%","20%","10%","05%","01%")
    label <- rep(y,13)
    q_table <- data.frame(name=label,q=vec,factor=factor_label)
    return(q_table)
}

BridgeRDatasetChecker <- function(InputFile, group, hour, InforColumn=4, OutputFile="BridgeR_2_Relative_RPKM_distribution"){
    ###Import_library###
    library(data.table)
    library(ggplot2)
    
    ###Prepare_files###
    time_points <- length(hour)
    
    input_file <- fread(InputFile, header=T)
    
    sample_size <- length(group)
    test_data <- NULL #input data for fig
    if(sample_size == 1){
        test_data <- input_file[T00_1 == 1,]
    }else if(sample_size == 2){
        test_data <- input_file[T00_1 == 1 & T00_2 == 1,]
    }else if(sample_size == 3){
        test_data <- input_file[T00_1 == 1 & T00_2 == 1 & T00_3 == 1,]
    }else if(sample_size == 4){
        test_data <- input_file[T00_1 == 1 & T00_2 == 1 & T00_3 == 1 & T00_4 == 1,]
    }
    
    ###Boxplot_for_each_sample###
    merge_fig_data <- NULL
    for(a in 1:sample_size){
        ###Information&exp_data column###
        infor_st <- 1 + (a - 1)*(time_points + InforColumn)
        infor_ed <- (InforColumn)*a + (a - 1)*time_points
        exp_st <- infor_ed + 1
        exp_ed <- infor_ed + time_points
        
        ###Hour_label###
        hour_label <- NULL
        for(x in hour){
            if(x == 0){
                next
            }
            label <- x
            if(x < 10){
                label <- paste("0",x,sep="")
            }
            hour_label <- append(hour_label, paste(label,"hr_",group[a], sep=""))
        }

        ###Prepare_exp_data###
        exp_st <- exp_st + 1
        exp_data <- test_data[,exp_st:exp_ed,with=F] #except time0
        exp_data <- t(exp_data) #Inverse
        exp_data <- factor(exp_data)
        exp_data <- as.numeric(as.character(exp_data)) #exp_data for fig_data
        exp_data <- log10(exp_data)
        
        ###Label_data###
        gene_number <- length(test_data[[1]])
        label_data <- rep(hour_label,gene_number) #label_data for fig_data
        
        ###Fig_data###
        fig_data <- data.frame(exp=exp_data, label=factor(label_data))
        if(a == 1){
            merge_fig_data <- fig_data
        }else{
            merge_fig_data <- rbind(merge_fig_data, fig_data)
        }
        
        ###fig_name###
        fig_name <- paste(OutputFile,"_Boxplot_",group[a],".png",sep="")
        fig_width <- 150*(time_points-1)
        png(filename=fig_name,width = fig_width, height = 1200)
        
        ###fig_plot###
        p <- ggplot()
        p <- p + layer(data=fig_data,
                       mapping=aes(x=label,y=exp),
                       geom="boxplot")
        p <- p + ylim(-2,2)
        plot(p)
        dev.off() #close_fig
        plot.new()
    }
    
    ###Boxplot_for_all_sample###
    ###fig_name###
    fig_name <- NULL
    for(a in 1:sample_size){
        if(a == 1){
            fig_name <- paste(OutputFile,"_Boxplot_",group[a],sep="")
        }else{
            fig_name <- paste(fig_name,"_",group[a],sep="")
        }
    }
    fig_name <- paste(fig_name,".png",sep="")
    fig_width <- 150*(time_points-1)*sample_size
    png(filename=fig_name,width = fig_width, height = 1200)
    
    ###fig_plot###
    merge_fig_data$label <- factor(merge_fig_data$label, levels=sort(unique(as.character(merge_fig_data$label))))
    p <- ggplot()
    p <- p + layer(data=merge_fig_data,
                   mapping=aes(x=label,y=exp),
                   geom="boxplot")
    p <- p + ylim(-2,2)
    plot(p)
    dev.off() #close_fig
    plot.new()
}

###Test###
setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_siStealth_siPUM2_ver1/time_course_0_1_2_4_8_12h")
BridgeRDatasetChecker(InputFile="BridgeR_1_Relative_expression_data_siCTRL_siPUM2_compatible.txt", group=c("siCTRL","siPUM2"), hour=c(0,1,2,4,8,12))
