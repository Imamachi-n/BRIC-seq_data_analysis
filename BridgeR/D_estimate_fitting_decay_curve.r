#Title: Estimate fitting decay curve
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-08

###Import_library###
library(data.table)

###Prepare_file_infor###
hour <- c(0,1,2,4,8,12)
time_point <- length(hour)
unit <- "hr"
data_path <- paste('C:/Users/Naoto/Documents/github/BRIC-seq_data_analysis/BridgeR/data/genes_RefSeq_result_mRNA',sep="")

###Load_input_file###
input_path <- paste(data_path,"_rel_normalized.fpkm_table",sep="")
output_path <- paste(data_path,"_rel_normalized_cal.fpkm_table",sep="")
input_data <- fread(input_path, header=T)
sample_size <- length(input_data[[1]])

###print_header###
cat("gr_id","symbol","accession_id","locus", sep="\t", file=output_path)
cat("\t", file=output_path, append=T)
hour_label <- NULL
for(x in hour){
    hour_label <- append(hour_label, paste(x, unit, sep=""))
}
cat(hour_label, sep="\t", file=output_path, append=T)
cat("\t", sep="", file=output_path, append=T)
cat("Model","Decay_rate_coef","coef_error","coef_p-value","R2","Adjusted_R2","Residual_standard_error","half_life", sep="\t", file=output_path, append=T)
cat("\n", sep="", file=output_path, append=T)

###calc_RNA_half_lives###
cutoff_rel_exp = 0.1
cutoff_data_point = 3

data_size <- length(as.vector(as.matrix(input_data[1,])))
for(x in 1:sample_size){
    data <- as.vector(as.matrix(input_data[x,])) ###Should be changed after!!### x -> 6160
    gene_infor <- data[1:4]
    cat(gene_infor, sep="\t", file=output_path, append=T)
    cat("\t", file=output_path, append=T)
    exp <- as.numeric(data[5:data_size])
    cat(exp, sep="\t", file=output_path, append=T)
    cat("\t", file=output_path, append=T)
    time_point_exp <- data.frame(hour,exp)
    time_point_exp <- time_point_exp[time_point_exp$exp >= cutoff_rel_exp, ]
    data_point <- length(time_point_exp$exp)
    if(!is.null(time_point_exp)){
        if(data_point >= cutoff_data_point){
            #Exponential Decay Model / f(t)=exp(-at)
            model <- lm(log(time_point_exp$exp) ~ time_point_exp$hour - 1)
            model_summary <- summary(model)
            coef <- -model_summary$coefficients[1]
            coef_error <- model_summary$coefficients[2]
            coef_p <- model_summary$coefficients[4]
            r_squared <- model_summary$r.squared
            adj_r_squared <- model_summary$adj.r.squared
            residual_standard_err <- model_summary$sigma
            half_life <- log(2)/coef
            if(coef < 0 || half_life >= 24){
                half_life <- 24
            }
            cat("Exponential_Decay_Model",coef,coef_error,coef_p,r_squared,adj_r_squared,residual_standard_err,half_life, sep="\t", file=output_path, append=T)
            cat("\t", file=output_path, append=T)
            #if(coef_p > 0.05){
            #    model_linear <- lm((time_point_exp$exp-1) ~ time_point_exp$hour - 1)
            #    model_linear_summary <- summary(model_linear)
            #    coef_linear <- model_linear_summary$coefficients[1]
            #    coef_error_linear <- model_linear_summary$coefficients[2]
            #    coef_p_linear <- model_linear_summary$coefficients[4]
            #    r_squared_linear <- model_linear_summary$r.squared
            #    adj_r_squared_linear <- model_linear_summary$adj.r.squared
            #    residual_standard_err_linear <- model_linear_summary$sigma
            #    half_life_linear <- -0.5/coef_linear
            #    cat("linear_Model",coef_linear,coef_error_linear,coef_p_linear,r_squared_linear,adj_r_squared_linear,residual_standard_err_linear,half_life_linear, sep="\t", file=output_path, append=T)
            #}
            cat("\n", file=output_path, append=T)
        }else{
            cat("few_data","\n", sep="", file=output_path, append=T)
        }
    }else{
        cat("low_expresion","\n", sep="", file=output_path, append=T)
    }
}

