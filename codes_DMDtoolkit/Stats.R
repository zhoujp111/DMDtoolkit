# read input file
#file_name<-"testing data_LVEDD"
indicator<-unlist(strsplit(file_name,"_"))[2]
data<-read.table(paste0(file_name,".rdata"),head=T,sep="\t")

# Stat 1. summary
sink(paste0(file_name,".sum"))
print(summary(data))
sink()

# Stat 2. correlation coefficient
data2<-data[,-c(1:2,9)] # remove columns 1,2 & 9
for (i in 1:ncol(data2)){
	data2[,i]<-as.numeric(data2[,i])
}

sink(paste0(file_name,".cor"))
for (i in 1:ncol(data2)){
	print(paste(indicator,names(data2)[i],sep=" and "))
	print(cor.test(data2[,indicator],data2[,i]))
}
sink()

# Stat 3. regression coefficient
y<-data[,indicator]
x<-data[,-c(1:2,9,which(names(data)%in%indicator))]

sink(paste0(file_name,".reg"))
	lm_1<-lm(y~.,data=x)
	print(names(data)[which(names(data)%in%indicator)])
	print(summary(lm_1))
	for (i in 1:ncol(x)){
		x_i<-x[i]
		lm_2<-lm(y~.,data=x_i)
		print(summary(lm_2))
	}
sink()

# Stat 4. t test of indicator subgroups according to indicator threshold
#indicator_threshold<-40 # set threshold for indicator
sink(paste0(file_name," ",indicator_threshold,".t-test"))
	for (i in 1:ncol(data2)){
		x1<-data2[data2[,indicator]<indicator_threshold & !is.na(data2[,i]),i]
		x2<-data2[data2[,indicator]>=indicator_threshold & !is.na(data2[,i]),i]
		print(paste0("Comparison for ",names(data2)[i]," between subgroups of ",indicator,"<",indicator_threshold," (n=",length(x1),") and ",indicator,">=",indicator_threshold," (n=",length(x2),")"))
		print(t.test(x1[!is.na(x1)],x2[!is.na(x2)]))
	}
sink()
