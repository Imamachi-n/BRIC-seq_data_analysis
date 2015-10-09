#Title: Decay rate infor
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-09

###current_directory###
setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data")

###input_file_infor###
hour <- c(0,1,2,4,8,12)
group <- c("siCTRL","siStealth")
files <- "C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_4_half-life_calculation.txt"

###Decay_rate_Infor_function###
BridgeRHalfLifeComparison <- function(filename = "BridgeR_4_half-life_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, OutputFig = "BridgeR_5_Half-life_comparison"){
    #ComparisonFile: The length of vector must be 2 => c("Control","Knockdown")
    ###Import_library###
    library(data.table)
    library(ggplot2)
    
    ###Prepare_file_infor###
    time_points <- length(hour)
    group_number <- length(group)
    comp_file_number <- NULL
    for(a in 1:length(ComparisonFile)){
        comp_file_number <- append(comp_file_number, which(group == ComparisonFile[a]))
    }
    
    input_file <- fread(filename, header=T)
    figfile <- paste(OutputFig,"_",group[comp_file_number[1]],"_vs_",group[comp_file_number[2]],".png", sep="")
    png(filename=figfile,width = 1200, height = 1200)
    
    ###Plot_Half-life_comparison###
    half_life_column_1 <- comp_file_number[1]*(time_points + InforColumn + 8) #number
    half_1 <- input_file[[half_life_column_1]]
    half_life_column_2 <- comp_file_number[2]*(time_points + InforColumn + 8) #number
    half_2 <- input_file[[half_life_column_2]]
        
    gene_number <- length(half_1)
    half_1_fig <- NULL
    half_2_fig <- NULL
    factor_fig <- NULL
    for(x in 1:gene_number){
        if(is.na(half_1[x]) || is.na(half_2[x])){
            next
        }
        half_1_fig <- append(half_1_fig, half_1[x])
        half_2_fig <- append(half_2_fig, half_2[x])
        div <- half_1[x]/half_2[x]
        if(div <= 0.5 || div >= 2){
            factor_fig <- append(factor_fig, 1)
        }else{
            factor_fig <- append(factor_fig, 0)
        }
    }
    half_1_fig <- log2(half_1_fig)
    half_2_fig <- log2(half_2_fig)
    counter <- length(which(factor_fig == 1))

    plot_data <- data.frame(half_1_fig,half_2_fig,factor_fig)
    print_out <- paste("Plotted: ",length(plot_data[,1])," genes", sep="")
    print_out2 <- paste("At least 2-fold: ",counter," genes", sep="")
    print(print_out)
    print(print_out2)

    p.scatter <- ggplot()
    p.scatter <- p.scatter + layer(data=plot_data, 
                                   mapping=aes(x=half_1_fig, y=half_2_fig,colour=factor(factor_fig)), 
                                   geom="point",
                                   #colour="black",
                                   size=2.5,
                                   alpha=0.3)
    p.scatter <- p.scatter + xlim(0,max(plot_data$half_1_fig)) + ylim(0,max(plot_data$half_2_fig))
    p.scatter <- p.scatter + ggtitle("Half-life comparison")
    name_xlab <- paste(group[comp_file_number[1]]," (Time)", sep="")
    name_ylab <- paste(group[comp_file_number[2]]," (Time)", sep="")
    p.scatter <- p.scatter + xlab(name_xlab)
    p.scatter <- p.scatter + ylab(name_ylab)
    p.scatter <- p.scatter + theme(legend.position="none") #Remove Guides
    p.scatter <- p.scatter + scale_colour_manual(values=c("black","red")) #Change color
    plot(p.scatter)
        
    dev.off() #close_fig
    plot.new()
}
    
BridgeRHalfLifeComparison(group=group, hour=hour, ComparisonFile=group)
