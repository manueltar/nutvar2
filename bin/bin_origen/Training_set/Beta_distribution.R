setwd("~/Escritorio/Proyecto_clasificador/Raw_Data/1_Antonio_Global_Set")

table.initial<-read.table("tabla_AFs_Beta_distribution.txt")
wt<-table.initial[,8]
HET<-table.initial[,9]
HOM<-table.initial[,10]
total.vector<-2*(wt+HET+HOM)
observed.vector<-HET+2*HOM
lower<-qbeta(0.025,1+observed.vector,1+total.vector-observed.vector, ncp=0)
upper<-qbeta(0.975,1+observed.vector,1+total.vector-observed.vector, ncp=0)

new.file<-cbind(table.initial,lower,upper)
write.table(new.file,"OUT_plus_Intervals_of_confidence.txt", sep="\t")