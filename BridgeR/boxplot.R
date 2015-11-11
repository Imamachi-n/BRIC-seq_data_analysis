setwd("C:/Users/Naoto/Documents/github/test_R_libraries/BridgeR_test/test_1hr_RPKM_halflife/R2_comp")

library(ggplot2)

HKgenes_MTR4 <- as.numeric(as.vector(as.matrix(read.table("R2_HKGenes_MTR4_study_siCTRL.txt"))))
HKgenes_PUM2 <- as.numeric(as.vector(as.matrix(read.table("R2_HKGenes_PUM2_study_siCTRL.txt"))))
Percentile_MTR4 <- as.numeric(as.vector(as.matrix(read.table("R2_Percentile_MTR4_study_siCTRL.txt"))))
Percentile_PUM2 <- as.numeric(as.vector(as.matrix(read.table("R2_Percentile_PUM2_study_siCTRL.txt"))))

R2_data <- c(HKgenes_MTR4,HKgenes_PUM2,Percentile_MTR4,Percentile_PUM2)

R2_HKGenes_label <- rep("House-keeping_genes",length(R2_HKGenes_data))
R2_Percentile_label <- rep("97.5-99th_Percentile",length(R2_Percentile_data))
R2_label <- c(R2_HKGenes_label,R2_Percentile_label)

facet1 <- rep("MTR4_study",length(HKgenes_MTR4))
facet2 <- rep("PUM2_study",length(HKgenes_PUM2))
facet3 <- rep("MTR4_study",length(Percentile_MTR4))
facet4 <- rep("PUM2_study",length(Percentile_PUM2))

facet_data <- c(facet1,facet2,facet3,facet4)

fig_data <- data.frame(data=R2_data,label=R2_label,facet=facet_data)
fig_data$label <- factor(fig_data$label, levels=c("House-keeping_genes","97.5-99th_Percentile"))
                                             
p <- ggplot()
p <- p + layer(data=fig_data, 
               mapping=aes(x=label, y=data, factor=label), 
               geom="boxplot")
p <- p + ylim(0,1) + facet_grid(~facet)
p <- p + ggtitle("")
p <- p + ylab("R2") + xlab("")
plot(p)
