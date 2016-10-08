# read input file
#file_name<-"testing data_LVEDD"
indicator<-unlist(strsplit(file_name,"_"))[2]
data<-read.table(paste0(file_name,".rdata"),head=T,sep="\t")

# Plot 2. plot frequency histogram of mutation types
plot.freq<-function(type,num){
	if (type == "all"){
		barplot(summary(data$Mutation)[1:num],las=2,cex.names=0.7)
	} else {
		barplot(summary(data$Mutation[grep(paste0("exon.*",type),data$Mutation)])[1:num],las=2,cex.names=0.7)
	}
}
# plot.freq("del",10) # an example of plot command

# Plot 3. scatter plot with trend line
plot.trend<-function(column_no1,column_no2){
	plot(data[,column_no1],data[,column_no2],xlab=names(data)[column_no1],ylab=names(data)[column_no2])
	abline(lm(data[,column_no2]~data[,column_no1]),col="black")
}
# plot.trend(3,6) # the 3rd and 6th column of input file

# Plot 4. stem and leaf plot
plot.stem<-function(column_no){
	sink(paste0(file_name,".stem"))
		#for (i in 1:ncol(data2)){
		#	print(names(data2)[i])
		#	stem(data2[,i])
		#}
		for (i in column_no){
			print(names(data)[i])
			stem(data[,i])
		}
	sink()
}
# plot.stem(c(3,6)) #  the 3rd and 6th column of input file

# Plot 5. cluster dendrogram
plot.clust<-function(column_no,cex_no){
	row.names(data2)<-data[,1]
	hc<-hclust(dist(data2[,column_no])^2,"cen")
	plot(hc,cex=cex_no,xlab="samples",main="Cluster of Samples")
#	memb<-cutree(hc,k)
#	cent <- NULL
#	for(i in 1:k){
#		cent <- rbind(cent, colMeans(data2[memb == i, , drop = FALSE]))
#	}
#	hc2<-hclust(dist(cent)^2,"cen",members=table(memb))
#	plot(hc2,cex=.5,xlab="samples",main="Cluster of Samples")
	savePlot("cluster of samples.pdf",type="pdf")
}
# plot.clust(1:6,0.4) # an example of cluster for samples
