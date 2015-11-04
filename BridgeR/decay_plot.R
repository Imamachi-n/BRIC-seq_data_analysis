setwd(paste("C:/Users/Naoto/Documents/github/akimitsu_sensei/fig",sep=""))

input_file <- read.table("C:/Users/Naoto/Documents/github/akimitsu_sensei/norm_figdata_for_R.txt")

hour_norm <- c(0,0.25,0.67,1,2,3,4,6,8,10,12)
hour_hypo <- c(0,0.25,0.50,1,2,3,4,6,8,10,12)


for(z in 1:length(input_file[,1])){
gene_symbol = as.character(as.vector(as.matrix(input_file[z,]))[1])
file_name <- sprintf("%1$s.png",gene_symbol)
png(filename=file_name,width = 640, height = 640)

model_norm = as.vector(as.matrix(input_file[z,]))[3]
model_hypo = as.vector(as.matrix(input_file[z,]))[46]

exp_norm <- as.vector(as.matrix(input_file[z,]))[7:17]
exp_hypo <- as.vector(as.matrix(input_file[z,]))[50:60]

fig_color <- c("black","red")

time_point_norm <- data.frame(hour=as.numeric(hour_norm),exp=as.numeric(exp_norm))
time_point_hypo <- data.frame(hour=as.numeric(hour_hypo),exp=as.numeric(exp_hypo))

library(ggplot2)
p.fitting <- ggplot()
p.fitting <- p.fitting + layer(data=time_point_norm, 
                               mapping=aes(x=hour, y=exp), 
                               geom="point",
                               size=4,
                               shape=19,
                               colour=fig_color[1])
p.fitting <- p.fitting + layer(data=time_point_hypo, 
                               mapping=aes(x=hour, y=exp), 
                               geom="point",
                               size=4,
                               shape=19,
                               colour=fig_color[2])

#################################
t <- model_norm
model1_a_norm <- as.numeric(as.vector(as.matrix(input_file[z,]))[27])
model2_a_norm <- as.numeric(as.vector(as.matrix(input_file[z,]))[33])
model2_b_norm <- as.numeric(as.vector(as.matrix(input_file[z,]))[34])
model3_a_norm <- as.numeric(as.vector(as.matrix(input_file[z,]))[39])
model3_b_norm <- as.numeric(as.vector(as.matrix(input_file[z,]))[40])
model3_c_norm <- as.numeric(as.vector(as.matrix(input_file[z,]))[41])

if(model_norm == "model1"){
    model1_curve<-function(t){exp(-model1_a_norm*t)}
    p.fitting <- p.fitting + layer(geom = "path",        # Default. Can be omitted.
                                   stat = "function",
                                   fun = model1_curve,
                                   mapping = aes(color = "Normoxia"),
                                   size=1.2)
}else if(model_norm == "model2"){
    model2_curve<-function(t){(1-model2_b_norm)*exp(-model2_a_norm*t)+model2_b_norm}
    p.fitting <- p.fitting + layer(geom = "path",        # Default. Can be omitted.
                                   stat = "function",
                                   fun = model2_curve,
                                   mapping = aes(color = "Normoxia"),
                                   size=1.2)
}else if(model_norm == "model3"){
    model3_curve<-function(t){model3_c_norm*exp(-model3_a_norm*t)+(1-model3_c_norm)*exp(-model3_b_norm*t)}
    p.fitting <- p.fitting + layer(geom = "path",        # Default. Can be omitted.
                                   stat = "function",
                                   fun = model3_curve,
                                   mapping = aes(color = "Normoxia"),
                                   size=1.2)
}

#################################
t <- model_hypo
model1_a_hypo <- as.numeric(as.vector(as.matrix(input_file[z,]))[70])
model2_a_hypo <- as.numeric(as.vector(as.matrix(input_file[z,]))[76])
model2_b_hypo <- as.numeric(as.vector(as.matrix(input_file[z,]))[77])
model3_a_hypo <- as.numeric(as.vector(as.matrix(input_file[z,]))[82])
model3_b_hypo <- as.numeric(as.vector(as.matrix(input_file[z,]))[83])
model3_c_hypo <- as.numeric(as.vector(as.matrix(input_file[z,]))[84])

if(model_hypo == "model1"){
    model1_curve<-function(t){exp(-model1_a_hypo*t)}
    p.fitting <- p.fitting + layer(geom = "path",        # Default. Can be omitted.
                                   stat = "function",
                                   fun = model1_curve,
                                   mapping = aes(color = "Hypoxia"),
                                   size=1.2)
}else if(model_hypo == "model2"){
    model2_curve<-function(t){(1-model2_b_hypo)*exp(-model2_a_hypo*t)+model2_b_hypo}
    p.fitting <- p.fitting + layer(geom = "path",        # Default. Can be omitted.
                                   stat = "function",
                                   fun = model2_curve,
                                   mapping = aes(color = "Hypoxia"),
                                   size=1.2)
}else if(model_hypo == "model3"){
    model3_curve<-function(t){model3_c_hypo*exp(-model3_a_hypo*t)+(1-model3_c_hypo)*exp(-model3_b_hypo*t)}
    p.fitting <- p.fitting + layer(geom = "path",        # Default. Can be omitted.
                                   stat = "function",
                                   fun = model3_curve,
                                   mapping = aes(color = "Hypoxia"),
                                   size=1.2)
}

p.fitting <- p.fitting + ggtitle(gene_symbol)
p.fitting <- p.fitting + xlab("Time")
p.fitting <- p.fitting + ylab("Relative RPKM (Time0 = 1)")
p.fitting <- p.fitting + scale_x_continuous(limits = c(0, 12))
p.fitting <- p.fitting + scale_color_manual(name = "Sample",
                                            values = c("red","black"), # Color specification
                                            labels = c("Hypoxia","Normoxia"))
#ybreaks <- seq(0,10,0.1)[2:101]
#p.fitting <- p.fitting + xlim(0,12) 
p.fitting <- p.fitting + scale_y_log10()
plot(p.fitting)
dev.off()
plot.new()
}
