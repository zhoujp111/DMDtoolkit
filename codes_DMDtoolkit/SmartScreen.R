pkgLoad <- function(x){
	if (!require(x,character.only = TRUE)){
		install.packages(x,dep=TRUE, repos='http://star-www.st-andrews.ac.uk/cran/')
		if(!require(x,character.only = TRUE)) stop("Package not found")
	}
	#now load library and suppress warnings
	#suppressPackageStartupMessages(library(x))
	library(x,character.only=TRUE)
}

Impute<-function(filename,method){

	data<-read.table(paste0(filename,".txt"),head=T,sep="\t")
	if (method == "rf"){
	pkgLoad("randomForest")
	data.imputed<-rfImpute(Mutation ~ .,data)
	write.table(data.imputed, paste0(filename,"_imputed.txt"), quote=F, sep="\t", col.names=T, row.names=F)
	}else if (method == "mean" | method == "median"){
	pkgLoad("e1071")
	data.imputed<-impute(data,method)
	write.table(data.imputed, paste0(filename,".imputed.txt"), quote=F, sep="\t", col.names=T, row.names=F)
	}
}

Weight<-function(filename,method){

	data<-read.table(paste0(filename,".txt"),head=T,sep="\t")
	weights<-lm(Mutation ~ ., data=data)$coefficients[2:ncol(data)]
	sink("weights_estimated.txt")
	print(abs(weights)/mean(abs(weights)))
	sink()
}

# Impute("testing data","rf")
# Weight("testing data_imputed","lm")
