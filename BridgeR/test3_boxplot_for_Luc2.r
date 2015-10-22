setwd("C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_Luc2_normalization")

input_data <- read.table("BridgeR_3Luc2_Normalized_expression_data_siStealth_siPUM1_compatible.txt",header=T)
test_data <- input_data[input_data$T0_1 == 1,]
test_data <- test_data[test_data$T0_2 == 1,]

png(filename='Boxplot_Comparison_of_Rel_RPKM_siStealth_vs_siPUM1_Luc2_norm.png',width = 1300, height = 1000)

boxplot(log10(test_data$T1_1),log10(test_data$T1_2),
        log10(test_data$T2_1),log10(test_data$T2_2),
        log10(test_data$T4_1),log10(test_data$T4_2),
        log10(test_data$T8_1),log10(test_data$T8_2),
        log10(test_data$T12_1),log10(test_data$T12_2),
        ylim=c(-3,1),
        names=c("1h_siStealth","1h_siPUM1",
                "2h_siStealth","2h_siPUM1",
                "4h_siStealth","4h_siPUM1",
                "8h_siStealth","8h_siPUM1",
                "12h_siStealth","12h_siPUM1"))

dev.off() #close_fig
plot.new()

####################

input_data <- read.table("BridgeR_3Luc2_Normalized_expression_data_siStealth_siCTRL_compatible.txt",header=T)
test_data <- input_data[input_data$T0_1 == 1,]
test_data <- test_data[test_data$T0_2 == 1,]

png(filename='Boxplot_Comparison_of_Rel_RPKM_siStealth_vs_siCTRL_Luc2_norm.png',width = 1300, height = 1000)

boxplot(log10(test_data$T1_1),log10(test_data$T1_2),
        log10(test_data$T2_1),log10(test_data$T2_2),
        log10(test_data$T4_1),log10(test_data$T4_2),
        log10(test_data$T8_1),log10(test_data$T8_2),
        log10(test_data$T12_1),log10(test_data$T12_2),
        ylim=c(-3,1),
        names=c("1h_siStealth","1h_siCTRL",
                "2h_siStealth","2h_siCTRL",
                "4h_siStealth","4h_siCTRL",
                "8h_siStealth","8h_siCTRL",
                "12h_siStealth","12h_siCTRL"))

dev.off() #close_fig
plot.new()
