setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_Luc2_normalization")

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

library(ggplot2)

############################

input_data <- read.table("BridgeR_3Luc2_Normalized_expression_data_siStealth_siPUM1_compatible.txt",header=T)
test_data <- input_data[input_data$T0_1 == 1,]
test_data <- test_data[test_data$T0_2 == 1,]

T1_1 <- test_q(log10(test_data$T1_1),"1h_Stealth")
T2_1 <- test_q(log10(test_data$T2_1),"2h_Stealth")
T4_1 <- test_q(log10(test_data$T4_1),"4h_Stealth")
T8_1 <- test_q(log10(test_data$T8_1),"8h_Stealth")
T12_1 <- test_q(log10(test_data$T12_1),"12h_Stealth")

T1_2 <- test_q(log10(test_data$T1_2),"1h_siPUM1")
T2_2 <- test_q(log10(test_data$T2_2),"2h_siPUM1")
T4_2 <- test_q(log10(test_data$T4_2),"4h_siPUM1")
T8_2 <- test_q(log10(test_data$T8_2),"8h_siPUM1")
T12_2 <- test_q(log10(test_data$T12_2),"12h_siPUM1")

fig_data <- rbind(T1_1,T1_2,T2_1,T2_2,T4_1,T4_2,T8_1,T8_2,T12_1,T12_2)

png(filename='Pointplot_Comparison_of_Rel_RPKM_siStealth_vs_siPUM1_Luc2_norm.png',width = 900, height = 900)

p <- ggplot()
p <- p + layer(data=fig_data, 
               mapping=aes(x=name, y=q, colour=factor(factor)), 
               geom="point",
               size=5,
               shape=19)
p <- p + xlab("") + ylab("Relative RPKM (Time0 = 1)")
plot(p)

dev.off() #close_fig
plot.new()

###################

input_data <- read.table("BridgeR_3Luc2_Normalized_expression_data_siStealth_siCTRL_compatible.txt",header=T)
test_data <- input_data[input_data$T0_1 == 1,]
test_data <- test_data[test_data$T0_2 == 1,]

T1_1 <- test_q(log10(test_data$T1_1),"1h_Stealth")
T2_1 <- test_q(log10(test_data$T2_1),"2h_Stealth")
T4_1 <- test_q(log10(test_data$T4_1),"4h_Stealth")
T8_1 <- test_q(log10(test_data$T8_1),"8h_Stealth")
T12_1 <- test_q(log10(test_data$T12_1),"12h_Stealth")

T1_2 <- test_q(log10(test_data$T1_2),"1h_siCTRL")
T2_2 <- test_q(log10(test_data$T2_2),"2h_siCTRL")
T4_2 <- test_q(log10(test_data$T4_2),"4h_siCTRL")
T8_2 <- test_q(log10(test_data$T8_2),"8h_siCTRL")
T12_2 <- test_q(log10(test_data$T12_2),"12h_siCTRL")

fig_data <- rbind(T1_1,T1_2,T2_1,T2_2,T4_1,T4_2,T8_1,T8_2,T12_1,T12_2)

png(filename='Pointplot_Comparison_of_Rel_RPKM_siStealth_vs_siCTRL_Luc2_norm.png',width = 900, height = 900)

p <- ggplot()
p <- p + layer(data=fig_data, 
               mapping=aes(x=name, y=q, colour=factor(factor)), 
               geom="point",
               size=5,
               shape=19)
p <- p + xlab("") + ylab("Relative RPKM (Time0 = 1)")
plot(p)

dev.off() #close_fig
plot.new()
