#Title: Draw fitting curve
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-08

###current_directory###
setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data")

###input_file_infor###
hour <- c(0,1,2,4,8,12)
group <- c("siCTRL","siStealth")
files <- "C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_3_Normalized_expression_data.txt"

###Draw_fitting_curve_function###
BridgeRDrawFittingCurve <- function(filename = "BridgeR_3_Normalized_expression_data.txt", group, hour, ComparisonFile, CutoffRelExp = 0.1, CutoffDataPoint = 3, InforColumn = 4, OutputDir = "BridgeR_fig"){
    ###Import_library###
    library(data.table)
    library(ggplot2)
    
    ###Make_stored_directory###
    ComparisonFile_name = paste(ComparisonFile,collapse="_")
    output_dir_name <- paste(OutputDir,ComparisonFile_name,sep="_")
    dir.create(output_dir_name)
    
    ###Prepare_file_infor###
    time_points <- length(hour)
    group_number <- length(group)
    input_file <- fread(filename, header=T)
    comp_file_number <- NULL
    for(a in 1:length(ComparisonFile)){
        comp_file_number <- append(comp_file_number, which(group == ComparisonFile[a]))
    }
    
    ###Draw_fitting_curve###
    setwd(paste("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data","/",output_dir_name,sep=""))
    gene_number <- length(input_file[[1]])
    for(x in 1:gene_number){
        data <- as.vector(as.matrix(input_file[x,]))
        
        #Save_fig
        gene_name <- as.character(data[2])
        file_name <- sprintf("%1$s.png",gene_name)
        paste(output_dir_name,"/",file_name,sep="")
        png(filename=file_name,width = 640, height = 640)
        
        #Prepare_ggplot2
        p.fitting <- ggplot()

        flg <- 0
        fig_color <- NULL
        for(a in comp_file_number){
            if(flg == 0){
                fig_color <- "black"
            }else{
                fig_color <- "red"
            }
            infor_st <- 1 + (a - 1)*(time_points + InforColumn)
            infor_ed <- (InforColumn)*a + (a - 1)*time_points
            exp_st <- infor_ed + 1
            exp_ed <- infor_ed + time_points
            
            exp <- as.numeric(data[exp_st:exp_ed])
            #hour <- c(0,1,2,4,8,12) ###TEST###
            #exp <- c(1,0.8350178,0.7806458,0.6386572,0.3946119,0.2603135)
            #exp <- c(1,0.9637603,1.131551,1.025119,1.246695,1.437431)
            #gene_name <- "test"
            #CutoffRelExp <- 0.1
            #CutoffDataPoint <- 3
            time_point_exp_original <- data.frame(hour,exp)
            
            p.fitting <- p.fitting + layer(data=time_point_exp_original, 
                                           mapping=aes(x=hour, y=exp), 
                                           geom="point",
                                           size=4,
                                           shape=19,
                                           colour=fig_color)
            
            time_point_exp <- time_point_exp_original[time_point_exp_original$exp >= CutoffRelExp, ] #>=0.1
            data_point <- length(time_point_exp$exp)
            if(!is.null(time_point_exp)){
                if(data_point >= CutoffDataPoint){ #>=3
                    model <- lm(log(time_point_exp$exp) ~ time_point_exp$hour - 1)
                    model_summary <- summary(model)
                    coef <- -model_summary$coefficients[1]
                    coef_error <- model_summary$coefficients[2]
                    coef_minus <- coef - coef_error
                    coef_plus <- coef + coef_error
                    
                    xmin <- min(hour[1])
                    xmax <- max(hour[length(hour)])
                    predicted_hour <- seq(xmin,xmax,length.out=100)
                    predicted <- data.frame(hour=predicted_hour)
                    predicted_minus <- data.frame(hour=predicted_hour)
                    predicted_plus <- data.frame(hour=predicted_hour)
                    
                    model1_pred <- function(x){
                        mRNA_exp <- exp(-x * predicted$hour)
                    }
                    predicted$exp <- model1_pred(coef)
                    predicted_minus$exp <- model1_pred(coef_minus)
                    predicted_plus$exp <- model1_pred(coef_plus)
                    
                    #p.fitting <- ggplot()
                    
                    p.fitting <- p.fitting + layer(data=predicted,
                                                   mapping=(aes(x=hour, y=exp)),
                                                   geom="line",
                                                   size=1.2,
                                                   colour=fig_color)
                    p.fitting <- p.fitting + layer(data=predicted_minus,
                                                   mapping=(aes(x=hour, y=exp)),
                                                   geom="line",
                                                   size=0.5,
                                                   linetype="dashed",
                                                   colour=fig_color)
                    p.fitting <- p.fitting + layer(data=predicted_plus,
                                                   mapping=(aes(x=hour, y=exp)),
                                                   geom="line",
                                                   size=0.5,
                                                   linetype="dashed",
                                                   colour=fig_color)
                    p.fitting <- p.fitting + ylab("Relative RPKM (Time0 = 1)")
                    p.fitting <- p.fitting + xlim(0,12)
                    ybreaks <- seq(0,10,0.1)[2:101]
                    p.fitting <- p.fitting + scale_y_log10(breaks=ybreaks,labels=ybreaks)
                    #p.fitting <- p.fitting + scale_y_log10()
                    plot(p.fitting)
                    #model_summary <- summary(model)
                    #coef <- -model_summary$coefficients[1]
                    #coef_error <- model_summary$coefficients[2]
                    #coef_p <- model_summary$coefficients[4]
                    #r_squared <- model_summary$r.squared
                    #adj_r_squared <- model_summary$adj.r.squared
                    #residual_standard_err <- model_summary$sigma
                    #half_life <- log(2)/coef
                    #if(coef < 0 || half_life >= 24){
                    #    half_life <- 24
                    #}
                    #cat("Exponential_Decay_Model",coef,coef_error,coef_p,r_squared,adj_r_squared,residual_standard_err,half_life, sep="\t", file=output_file, append=T)
                }else{
                    #cat("few_data","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                }
            }else{
                #cat("low_expresion","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
            } ###TEST###
            flg = 1
        }
        #cat("\n", file=output_file, append=T)
        dev.off() #close_fig
        plot.new()
    }
}

BridgeRDrawFittingCurve(group=group, hour=hour, ComparisonFile=group)
