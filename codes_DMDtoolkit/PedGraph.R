# Plot 1. plot pedigree
plot.ped<-function(filename){
	dataPed<-read.table(paste0(filename,".txt"),head=T,sep="\t")
	
	setwd("./kinship/")
	source("plot.pedigree.R")
	source("pedigree.R")
	source("autohint.R")
	source("align.pedigree.R")
	source("kindepth.R")
	source("check.hint.R")
	source("alignped1.R")
	source("alignped2.R")
	source("alignped3.R")
	source("alignped4.R")
	
	setwd("../")
	dataPed$famid<-as.factor(dataPed$famid)
	
	for (i in levels(dataPed$famid)){
		dataPed2<-dataPed[dataPed$famid==i,]
		#attach(dataPed2)
		pdf(paste0("family_",i,".pdf"),width=4,height=3) # create and open a graph file
		ped<-pedigree(dataPed2$id,dataPed2$fid,dataPed2$mid,dataPed2$sex,dataPed2$aff)
		#par(xpd=T)
		#plot.pedigree(ped)
		
		strid<-dataPed2$id
		strid<-paste(dataPed2$id,dataPed2$mutation,sep="\n")
		par(xpd=T)
		#plot(ped,id=strid)
		plot(ped,id=strid,mar=c(3,3,3,3),cex=0.5,lwd=5)
		
		dev.off() # close the graph file
		#detach(data2)
	}
}
# plot.ped("pedigree") # input filename is "pedigree"
