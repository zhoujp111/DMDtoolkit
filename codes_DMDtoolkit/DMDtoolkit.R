################################
# 
#       by Jiapeng Zhou
#      zhoujp111@126.com
#
################################

########## Drawing mutated DMD proteins ##########
Args<-commandArgs()
data3<-read.table(paste0(Args[1].".Dp427m.rdata"),head=T,sep="\n") #reading clinical data
#data4<-read.table(paste0(Args[1].".Dp427m.pros.diag",head=T,sep="\t") #reading predition 1
data4<-read.table(paste0(Args[1].".Dp427m.pros.diag2",head=T,sep="\t") #reading predition 2

for (i in 1:nrow(data3)){
	str_split<-unlist(strsplit(as.character(unlist(data3[i,])),"\t"))
	pdf(paste0(str_split[1],".pdf"),width=26,height=13) #creating graph file
	par(bg="white",mar=c(0.2,0.2,0.2,0.2)) #seting background
	plot(x=c(0),y=c(0),xlim=c(0,26),ylim=c(0,13),col="white",axes=F,ann=F,pch=20) #seting coordinate
	text(13,13,paste0(data4[i,1],": length of potential protein ",data4[i,6],"; number of potential stop-gains ",data4[i,7],"; prediction of reading-frame rule ",data4[i,8]),cex=2,col="black") #seting title

### function of drawing introns ###
	segments(-2,1,10,1,col="grey",lty=1,lwd=2)
	for (i in 2:6){
		abline(h=2*i-1,col="grey",lwd=2)
	}

### function of drawing exons ###

	draw.exon<-function(exon_region,pos_region1,pos_region2,col_region){
	
	#drawing polygons/exons
		draw.polygon<-function(x1,x2,y,z1,z2,n,l,c){
			polygon(c(x1+0.3*exon_region,x2+0.3*exon_region,x2+0.3*exon_region+0.07*z2,x2+0.3*exon_region+0.07*2*z2,x2+0.3*exon_region+0.07*z2,x2+0.3*exon_region,x1+0.3*exon_region,x1+0.3*exon_region+0.07*z1,x1+0.3*exon_region+0.07*2*z1,x1+0.3*exon_region+0.07*z1),c(y+0.3,y+0.3,y+0.15,y,y-0.15,y-0.3,y-0.3,y-0.15,y,y+0.15), col=c, border="blue", density=c(1000))
			text((x1+x2)/2+0.3*exon_region,y+0.15,n,col="black")
			text((x1+x2)/2+0.3*exon_region,y-0.15,l,col="red",cex=0.6)
		}
		draw.polygon2<-function(x1,x2,y,z1,z2,n,l,c){
			polygon(c(x1+0.3*exon_region,x2+0.3*exon_region,x2+0.3*exon_region+0.03*z2,x2+0.3*exon_region+0.03*2*z2,x2+0.3*exon_region+0.03*z2,x2+0.3*exon_region,x1+0.3*exon_region,x1+0.3*exon_region+0.03*z1,x1+0.3*exon_region+0.03*2*z1,x1+0.3*exon_region+0.03*z1),c(y+0.3,y+0.3,y+0.15,y,y-0.15,y-0.3,y-0.3,y-0.15,y,y+0.15), col=c, border="blue", density=c(1000))
			text((x1+x2)/2+0.3*exon_region,y+0.15,n,col="black")
			text((x1+x2)/2+0.3*exon_region,y-0.15,l,col="red",cex=0.6)
		}
		draw.polygon3<-function(x1,x2,y,z1,z2,n,l,c){
			polygon(c(x1+0.3*exon_region,x2+0.3*exon_region,x2+0.3*exon_region+0.07*z2,x2+0.3*exon_region+0.07*2*z2,x2+0.3*exon_region+0.07*z2,x2+0.3*exon_region,x1+0.3*exon_region,x1+0.3*exon_region+0.03*z1,x1+0.3*exon_region+0.03*2*z1,x1+0.3*exon_region+0.03*z1),c(y+0.3,y+0.3,y+0.15,y,y-0.15,y-0.3,y-0.3,y-0.15,y,y+0.15), col=c, border="blue", density=c(1000))
			text((x1+x2)/2+0.3*exon_region,y+0.15,n,col="black")
			text((x1+x2)/2+0.3*exon_region,y-0.15,l,col="red",cex=0.6)
		}
		draw.polygon4<-function(x1,x2,y,z1,z2,n,l,c){
			polygon(c(x1+0.3*exon_region,x2+0.3*exon_region,x2+0.3*exon_region+0.03*z2,x2+0.3*exon_region+0.03*2*z2,x2+0.3*exon_region+0.03*z2,x2+0.3*exon_region,x1+0.3*exon_region,x1+0.3*exon_region+0.07*z1,x1+0.3*exon_region+0.07*2*z1,x1+0.3*exon_region+0.07*z1),c(y+0.3,y+0.3,y+0.15,y,y-0.15,y-0.3,y-0.3,y-0.15,y,y+0.15), col=c, border="blue", density=c(1000))
			text((x1+x2)/2+0.3*exon_region,y+0.15,n,col="black")
			text((x1+x2)/2+0.3*exon_region,y-0.15,l,col="red",cex=0.6)
		}
		
	#cycling
		pos_exon_1<-round(pos_region1/100,2)
		pos_exon_2<-round(pos_region2/100,2)
		len_exon<-pos_region2-pos_region1+1
		index_1<-(pos_region1-1)%%3
		index_2<-pos_region2%%3
		if (exon_region <= 16){
			if (index_1 == 1 & index_2 == 1){
				draw.polygon2(pos_exon_1,pos_exon_2,11,index_1,index_2,exon_region,len_exon,col_region)
			}else if (index_1 == 1 & (index_2 == 0 | index_2 == 2)){
				draw.polygon3(pos_exon_1,pos_exon_2,11,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & index_2 == 1){
				draw.polygon4(pos_exon_1,pos_exon_2,11,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & (index_2 == 0 | index_2 == 2)){
				draw.polygon(pos_exon_1,pos_exon_2,11,index_1,index_2,exon_region,len_exon,col_region)
			}
		}
		if (exon_region >= 17 & exon_region <= 30){
			if (index_1 == 1 & index_2 == 1){
				draw.polygon2(pos_exon_1-25*1,pos_exon_2-25*1,9,index_1,index_2,exon_region,len_exon,col_region)
			}else if (index_1 == 1 & (index_2 == 0 | index_2 == 2)){
				draw.polygon3(pos_exon_1-25*1,pos_exon_2-25*1,9,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & index_2 == 1){
				draw.polygon4(pos_exon_1-25*1,pos_exon_2-25*1,9,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & (index_2 == 0 | index_2 == 2)){
				draw.polygon(pos_exon_1-25*1,pos_exon_2-25*1,9,index_1,index_2,exon_region,len_exon,col_region)
			}
		}
		if (exon_region >= 31 & exon_region <= 43){
			if (index_1 == 1 & index_2 == 1){
				draw.polygon2(pos_exon_1-25*2,pos_exon_2-25*2,7,index_1,index_2,exon_region,len_exon,col_region)
			}else if (index_1 == 1 & (index_2 == 0 | index_2 == 2)){
				draw.polygon3(pos_exon_1-25*2,pos_exon_2-25*2,7,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & index_2 == 1){
				draw.polygon4(pos_exon_1-25*2,pos_exon_2-25*2,7,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & (index_2 == 0 | index_2 == 2)){
				draw.polygon(pos_exon_1-25*2,pos_exon_2-25*2,7,index_1,index_2,exon_region,len_exon,col_region)
			}
		}
		if (exon_region >= 44 & exon_region <= 56){
			if (index_1 == 1 & index_2 == 1){
				draw.polygon2(pos_exon_1-25*3,pos_exon_2-25*3,5,index_1,index_2,exon_region,len_exon,col_region)
			}else if (index_1 == 1 & (index_2 == 0 | index_2 == 2)){
				draw.polygon3(pos_exon_1-25*3,pos_exon_2-25*3,5,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & index_2 == 1){
				draw.polygon4(pos_exon_1-25*3,pos_exon_2-25*3,5,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & (index_2 == 0 | index_2 == 2)){
				draw.polygon(pos_exon_1-25*3,pos_exon_2-25*3,5,index_1,index_2,exon_region,len_exon,col_region)
			}
		}
		if (exon_region >= 57 & exon_region <= 73){
			if (index_1 == 1 & index_2 == 1){
				draw.polygon2(pos_exon_1-25*4,pos_exon_2-25*4,3,index_1,index_2,exon_region,len_exon,col_region)
			}else if (index_1 == 1 & (index_2 == 0 | index_2 == 2)){
				draw.polygon3(pos_exon_1-25*4,pos_exon_2-25*4,3,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & index_2 == 1){
				draw.polygon4(pos_exon_1-25*4,pos_exon_2-25*4,3,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & (index_2 == 0 | index_2 == 2)){
				draw.polygon(pos_exon_1-25*4,pos_exon_2-25*4,3,index_1,index_2,exon_region,len_exon,col_region)
			}
		}
		if (exon_region >= 74 & exon_region <= 79){
			if (index_1 == 1 & index_2 == 1){
				draw.polygon2(pos_exon_1-25*5,pos_exon_2-25*5,1,index_1,index_2,exon_region,len_exon,col_region)
			}else if (index_1 == 1 & (index_2 == 0 | index_2 == 2)){
				draw.polygon3(pos_exon_1-25*5,pos_exon_2-25*5,1,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & index_2 == 1){
				draw.polygon4(pos_exon_1-25*5,pos_exon_2-25*5,1,index_1,index_2,exon_region,len_exon,col_region)
			}else if ((index_1 == 0 | index_1 == 2) & (index_2 == 0 | index_2 == 2)){
				draw.polygon(pos_exon_1-25*5,pos_exon_2-25*5,1,index_1,index_2,exon_region,len_exon,col_region)
			}
		}
	}
	
### function of drawing Dp427m protein ###

	pro_isoform<-"Dp427m"
	data<-read.table(paste0(pro_isoform," CDs.txt"), head=T, sep = "\t") #¶ÁÈëDp427m CDsÊý¾Ý
	for (i in data[1,1]:data[nrow(data),1]){
		draw.exon(i,data[i,4],data[i,5],"white")
	}
	
### function of drawing promoter ###

	draw.promoter<-function(pro_isoform,col_isoform){
	
	#Dp427m
		if (pro_isoform == "Dp427" | pro_isoform == "Dp427m"){
			segments(0,11,0,12.5,col=col_isoform,lty=1,lwd=2)
			arrows(0,12.5,0.5,12.5,length=0.1,col=col_isoform,lty=1,lwd=2)
			text(1.0,12.5,"Dp427m",cex=1,col=col_isoform)
			segments(0.1,11,0.1,12.2,col=col_isoform,lty=2,lwd=2)
			arrows(0.1,12.2,0.5,12.2,length=0.1,col=col_isoform,lty=2,lwd=2)
			text(1.0,12.2,"Dp427m",cex=1,col=col_isoform)
		}
	}
	
### function of drawing domains/motifs ###

	data2<-read.table(paste0(pro_isoform," Domains.rdata"),head=T,sep="\t") 
	for (i in 1:nrow(data2)){
		if (data2[i,4] <= 16){
			arrows(0.3*data2[i,4]+data2[i,2]*3/100,11+0.45,0.3*data2[i,5]+data2[i,3]*3/100,11+0.45,length=0.1,col="black",lty=1,lwd=2)
			arrows(0.3*data2[i,5]+data2[i,3]*3/100,11+0.45,0.3*data2[i,4]+data2[i,2]*3/100,11+0.45,length=0.1,col="black",lty=1,lwd=2)
			text((0.3*data2[i,4]+data2[i,2]*3/100+0.3*data2[i,5]+data2[i,3]*3/100)/2,11+0.7,data2[i,1],cex=1,col="black")
		}else if (data2[i,4] >= 17 & data2[i,4] <= 29){
			arrows(0.3*data2[i,4]+data2[i,2]*3/100-25*1,9+0.45,0.3*data2[i,5]+data2[i,3]*3/100-25*1,9+0.45,length=0.1,col="black",lty=1,lwd=2)
			arrows(0.3*data2[i,5]+data2[i,3]*3/100-25*1,9+0.45,0.3*data2[i,4]+data2[i,2]*3/100-25*1,9+0.45,length=0.1,col="black",lty=1,lwd=2)
			text((0.3*data2[i,4]+data2[i,2]*3/100+0.3*data2[i,5]+data2[i,3]*3/100)/2-25*1,9+0.7,data2[i,1],cex=1,col="black")
		}else if (data2[i,4] == 30){
			arrows(round(data[30,5]/100,2)+0.3*30-25*1,9+0.45,0.3*data2[i,4]+data2[i,2]*3/100-25*1,9+0.45,length=0.1,col="black",lty=1,lwd=2)
			text((round(data[30,5]/100,2)+0.3*30+0.3*data2[i,4]+data2[i,2]*3/100)/2-25*1,9+0.7,"repeat",cex=1,col="black")
			arrows(round(data[31,4]/100,2)+0.3*31-25*2,7+0.45,0.3*data2[i,5]+data2[i,3]*3/100-25*2,7+0.45,length=0.1,col="black",lty=1,lwd=2)
			text((round(data[31,4]/100,2)+0.3*31+0.3*data2[i,5]+data2[i,3]*3/100)/2-25*2,7+0.7,"10",cex=1,col="black")
		}else if (data2[i,4] >= 31 & data2[i,4] <= 43){
			arrows(0.3*data2[i,4]+data2[i,2]*3/100-25*2,7+0.45,0.3*data2[i,5]+data2[i,3]*3/100-25*2,7+0.45,length=0.1,col="black",lty=1,lwd=2)
			arrows(0.3*data2[i,5]+data2[i,3]*3/100-25*2,7+0.45,0.3*data2[i,4]+data2[i,2]*3/100-25*2,7+0.45,length=0.1,col="black",lty=1,lwd=2)
			text((0.3*data2[i,4]+data2[i,2]*3/100+0.3*data2[i,5]+data2[i,3]*3/100)/2-25*2,7+0.7,data2[i,1],cex=1,col="black")
		}else if (data2[i,4] >= 44 & data2[i,4] <= 56){
			arrows(0.3*data2[i,4]+data2[i,2]*3/100-25*3,5+0.45,0.3*data2[i,5]+data2[i,3]*3/100-25*3,5+0.45,length=0.1,col="black",lty=1,lwd=2)
			arrows(0.3*data2[i,5]+data2[i,3]*3/100-25*3,5+0.45,0.3*data2[i,4]+data2[i,2]*3/100-25*3,5+0.45,length=0.1,col="black",lty=1,lwd=2)
			text((0.3*data2[i,4]+data2[i,2]*3/100+0.3*data2[i,5]+data2[i,3]*3/100)/2-25*3,5+0.7,data2[i,1],cex=1,col="black")
		}else if (data2[i,4] >= 57 & data2[i,4] <= 73){
			arrows(0.3*data2[i,4]+data2[i,2]*3/100-25*4,3+0.45,0.3*data2[i,5]+data2[i,3]*3/100-25*4,3+0.45,length=0.1,col="black",lty=1,lwd=2)
			arrows(0.3*data2[i,5]+data2[i,3]*3/100-25*4,3+0.45,0.3*data2[i,4]+data2[i,2]*3/100-25*4,3+0.45,length=0.1,col="black",lty=1,lwd=2)
			text((0.3*data2[i,4]+data2[i,2]*3/100+0.3*data2[i,5]+data2[i,3]*3/100)/2-25*4,3+0.7,data2[i,1],cex=1,col="black")
		}else if (data2[i,4] >= 74 & data2[i,4] <= 79){
			arrows(0.3*data2[i,4]+data2[i,2]*3/100-25*5,1+0.45,0.3*data2[i,5]+data2[i,3]*3/100-25*5,1+0.45,length=0.1,col="black",lty=1,lwd=2)
			arrows(0.3*data2[i,5]+data2[i,3]*3/100-25*5,1+0.45,0.3*data2[i,4]+data2[i,2]*3/100-25*5,1+0.45,length=0.1,col="black",lty=1,lwd=2)
			text((0.3*data2[i,4]+data2[i,2]*3/100+0.3*data2[i,5]+data2[i,3]*3/100)/2-25*5,1+0.7,data2[i,1],cex=1,col="black")
		}
	}
	
### function of drawing annotations ###

	draw.annotation<-function(){

		segments(11,0,26,0,col="black",lty=1,lwd=2)
		segments(11,2,26,2,col="black",lty=1,lwd=2)
		segments(11,0,11,2,col="black",lty=1,lwd=2)
		segments(26,0,26,2,col="black",lty=1,lwd=2)
		
		segments(11.3,1.1,11.3,1.7,col="red",lty=1,lwd=2)
		arrows(11.3,1.7,11.8,1.7,length=0.1,col="red",lty=1,lwd=2)
		text(12.2,1.4,"promotor",cex=0.8,col="black")
		segments(13.2,1.1,13.2,1.7,col="red",lty=2,lwd=2)
		arrows(13.2,1.7,13.7,1.7,length=0.1,col="red",lty=2,lwd=2)
		text(14.2,1.5,"traduction",cex=0.8,col="black")
		text(14.3,1.3,"initiation site",cex=0.8,col="black")
		
		arrows(11.3,0.4,12.8,0.4,length=0.1,col="black",lty=1,lwd=2)
		arrows(12.8,0.4,11.3,0.4,length=0.1,col="black",lty=1,lwd=2)
		text(12.05,0.6,"CH",cex=0.8,col="black")
		text(13.8,0.5,"structural domain",cex=0.8,col="black")
		
		draw.polygon2<-function(x1,x2,y,z1,z2,n,l,c){
			polygon(c(x1+0.3*exon_region,x2+0.3*exon_region,x2+0.3*exon_region+0.03*z2,x2+0.3*exon_region+0.03*2*z2,x2+0.3*exon_region+0.03*z2,x2+0.3*exon_region,x1+0.3*exon_region,x1+0.3*exon_region+0.03*z1,x1+0.3*exon_region+0.03*2*z1,x1+0.3*exon_region+0.03*z1),c(y+0.3,y+0.3,y+0.15,y,y-0.15,y-0.3,y-0.3,y-0.15,y,y+0.15), col=c, border="blue", density=c(1000))
			text((x1+x2)/2+0.3*exon_region,y+0.15,n,col="black")
			text((x1+x2)/2+0.3*exon_region,y-0.15,l,col="red",cex=0.6)
		}
		exon_region<-1
		draw.polygon2(15.2,15.5,1,0,1,1,31,"green")
		text(15.9,1.6,"exon number",cex=0.8,col="black")
		text(15.8,0.4,"exon size",cex=0.8,col="red")
		
		exon_region<-1
		draw.polygon2(16.8,17.1,1.5,0,1,1,31,"grey")
		text(17.9,1.5,"deletion",cex=0.8,col="black")

		exon_region<-1
		draw.polygon2(16.8,17.1,0.5,0,1,1,31,"darkgreen")
		text(17.95,0.6,"duplication",cex=0.8,col="black")
		text(17.95,0.4,"/insertion",cex=0.8,col="black")
		
		exon_region<-1
		draw.polygon2(18.55,18.85,1.5,0,1,1,31,"green")
		exon_region<-2
		draw.polygon2(18.65,19.25,1.5,1,0,2,62,"green")
		text(20.3,1.5,"in frame",cex=0.8,col="black")
		
		exon_region<-1
		draw.polygon2(18.55,18.85,0.5,0,1,1,31,"green")
		exon_region<-3
		draw.polygon2(18.4,19.3,0.5,0,0,3,93,"green")
		text(20.75,0.5,"frameshift",cex=0.8,col="black")
		
		text(22.4,1.6,"c.3727_3728delCT",cex=1,col="red")
		text(22.4,1.3,"c.6077_6078insA",cex=1,col="red")
		text(22.4,1,"c.3772dupT",cex=1,col="red")
		text(22.4,0.7,"c.3257dup",cex=1,col="red")
		text(22.4,0.4,"indel",cex=0.8,col="black")
		
		text(24.8,1.6,"c.433C>T p.Arg145*",cex=1,col="red")
		text(24.8,1.3,"point mutation",cex=0.8,col="black")
		text(24.8,0.7,"c.31+1G>T",cex=1,col="red")
		text(24.8,0.4,"splice site",cex=0.8,col="black")
		
	}

########## function of drawing mutated proteins ##########
	
	cds_range<-as.numeric(unlist(strsplit(as.character(unlist(str_split[length(str_split)]))," ")))
	cds_range0<-cds_range[2] - cds_range[1] + 1
	
	draw.promoter(pro_isoform, "green") #drawing promoter
	
	#marking the protein structure
	for (j in 1:((length(str_split)-6)/2)){
		exon_region<-unlist(strsplit(as.character(unlist(str_split[j*2+4]))," "))
		pos_region<-as.numeric(unlist(strsplit(as.character(unlist(str_split[j*2+5]))," ")))
		if (exon_region[3] == "promoter"){
			draw.promoter(as.numeric(exon_region[1]),"grey")
			for (i in data[1,1]:data[nrow(data),1]){
				draw.exon(i,data[i,4],data[i,5],"grey")
			}
			break
		}else if (exon_region[3] == "exon_del" | exon_region[3] == "del" | exon_region[3] == "as"){
			if (pos_region[1] == 1){
				cds_range[1] = cds_range[1] + pos_region[2] - pos_region[1] + 1
				cds_range[2] = cds_range[2] + pos_region[2] - pos_region[1] + 1
			}else{
				cds_range[2] = cds_range[2] + pos_region[2] - pos_region[1] + 1
			}
			if (cds_range[2] - cds_range[1] + 1 >= cds_range0){
				for (i in 1:nrow(data)){
					if (data[i,4] <= cds_range[1] & data[i,5] >= cds_range[1]){
						m<-i
					}
					if (data[i,4] <= cds_range[2] & data[i,5] >= cds_range[2]){
						n<-i
					}
				}
			}
		}else if (exon_region[3] == "exon_dup" | exon_region[3] == "dup" | exon_region[3] == "ins" | exon_region[3] == "point"){
			cds_range[2] = cds_range[2] + pos_region[1] - pos_region[2] - 1
			if (cds_range[2] - cds_range[1] + 1 <= cds_range0){
				for (i in 1:nrow(data)){
					if (data[i,4] <= cds_range[1] & data[i,5] >= cds_range[1]){
						m<-i
					}
					if (data[i,4] <= cds_range[2] & data[i,5] >= cds_range[2]){
						n<-i
					}
				}
			}
		}
	}
	
	#drawing the protein structure
	if (n - m > 1){
		for (i in (m+1):(n-1)){
			draw.exon(i, data[i,4], data[i,5], "green")
		}
	}else if (n - m == 1){
		draw.exon(m, cds_range[1], data[m,5], "green")
		draw.exon(n, data[n,4], cds_range[2], "green")
	}else if (n - m == 0){
		draw.exon(m, cds_range[1], cds_range[2], "green")
	}
	for (j in 1:((length(str_split)-6)/2)){
		mutations<-unlist(strsplit(as.character(unlist(str_split[5])),"; "))
		exon_region<-unlist(strsplit(as.character(unlist(str_split[j*2+4]))," "))
		pos_region<-as.numeric(unlist(strsplit(as.character(unlist(str_split[j*2+5]))," ")))
		if (exon_region[3] == "promoter"){
			break
		}else if (exon_region[3] == "exon_del" | exon_region[3] == "as"){
			if (n - m >= 1){
				draw.exon(m, cds_range[1], data[m,5], "green")
				draw.exon(n, data[n,4], cds_range[2], "green")
			}
			for (i in exon_region[1]:exon_region[2]){
				draw.exon(i, data[i,4], data[i,5], "grey")
			}
		}else if (exon_region[3] == "exon_dup"){
			for (i in exon_region[1]:exon_region[2]){
				draw.exon(i, data[i,4], data[i,5], "darkgreen")
			}
			if (n - m >= 1){
				draw.exon(m, cds_range[1], data[m,5], "green")
				draw.exon(n, data[n,4], cds_range[2], "green")
			}
		}else if (exon_region[3] == "ins" | exon_region[3] == "dup"){
			draw.exon(m, cds_range[1], data[m,5], "green")
			draw.exon(n, data[n,4], cds_range[2], "green")
			draw.exon(as.numeric(exon_region[1]),pos_region[1],pos_region[2],"darkgreen")
		}else if (exon_region[3] == "del"){
			draw.exon(m, cds_range[1], data[m,5], "green")
			draw.exon(n, data[n,4], cds_range[2], "green")
			draw.exon(as.numeric(exon_region[1]),pos_region[1],pos_region[2],"grey")
		}else if (exon_region[3] == "point"){
			draw.exon(m, cds_range[1], data[m,5], "green")
			draw.exon(n, data[n,4], cds_range[2], "green")
		}
		if (exon_region[3] == "as"){
			for (i in 1:nrow(data2)){
				if (as.numeric(exon_region[1]) <= 16){
					if (grepl("-",mutations[j])){
						segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100,11-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100,11+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[1]/100,11-0.3-0.15*j,mutations[j],cex=1,col="red")
					}else{
						segments(0.3*as.numeric(exon_region[1])+pos_region[2]/100,11-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[2]/100,11+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[2]/100,11-0.3-0.15*j,mutations[j],cex=1,col="red")
					}
				}else if (as.numeric(exon_region[1]) >= 17 & as.numeric(exon_region[1]) <= 30){
					if (grepl("-",mutations[j])){
						segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*1,9-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*1,9+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*1,9-0.3-0.15*j,mutations[j],cex=1,col="red")
					}else{
						segments(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*1,9-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*1,9+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*1,9-0.3-0.15*j,mutations[j],cex=1,col="red")
					}
				}else if (as.numeric(exon_region[1]) >= 30 & as.numeric(exon_region[1]) <= 43){
					if (grepl("-",mutations[j])){
						segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*2,7-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*2,7+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*2,7-0.3-0.15*j,mutations[j],cex=1,col="red")
					}else{
						segments(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*2,7-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*2,7+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*2,7-0.3-0.15*j,mutations[j],cex=1,col="red")
					}
				}else if (as.numeric(exon_region[1]) >= 44 & as.numeric(exon_region[1]) <= 56){
					if (grepl("-",mutations[j])){
						segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*3,5-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*3,5+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*3,5-0.3-0.15*j,mutations[j],cex=1,col="red")
					}else{
						segments(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*3,5-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*3,5+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*3,5-0.3-0.15*j,mutations[j],cex=1,col="red")
					}
				}else if (as.numeric(exon_region[1]) >= 57 & as.numeric(exon_region[1]) <= 73){
					if (grepl("-",mutations[j])){
						segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*4,3-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*4,3+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*4,3-0.3-0.15*j,mutations[j],cex=1,col="red")
					}else{
						segments(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*4,3-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*4,3+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*4,3-0.3-0.15*j,mutations[j],cex=1,col="red")
					}
				}else if (as.numeric(exon_region[1]) >= 74 & as.numeric(exon_region[1]) <= 79){
					if (grepl("-",mutations[j])){
						segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*5,1-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*5,1+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*5,1-0.3-0.15*j,mutations[j],cex=1,col="red")
					}else{
						segments(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*5,1-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*5,1+0.3,col="black",lty=1,lwd=2)
						text(0.3*as.numeric(exon_region[1])+pos_region[2]/100-25*5,1-0.3-0.15*j,mutations[j],cex=1,col="red")
					}
				}
			}
		}else if(exon_region[3] == "point"){
			for (i in 1:nrow(data2)){
				if (as.numeric(exon_region[1]) <= 16){
					segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100,11-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100,11+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100,11-0.3-0.15*j,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 17 & as.numeric(exon_region[1]) <= 30){
					segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*1,9-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*1,9+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*1,9-0.3-0.15*j,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 30 & as.numeric(exon_region[1]) <= 43){
					segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*2,7-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*2,7+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*2,7-0.3-0.15*j,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 44 & as.numeric(exon_region[1]) <= 56){
					segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*3,5-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*3,5+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*3,5-0.3-0.15*j,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 57 & as.numeric(exon_region[1]) <= 73){
					segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*4,3-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*4,3+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*4,3-0.3-0.15*j,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 74 & as.numeric(exon_region[1]) <= 79){
					segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*5,1-0.25-0.15*j,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*5,1+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*5,1-0.3-0.15*j,mutations[j],cex=1,col="red")
				}
			}
		}else{
			for (i in 1:nrow(data2)){
				if (as.numeric(exon_region[1]) <= 16){
					#segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100,11-0.4,0.3*as.numeric(exon_region[1])+pos_region[1]/100,11+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100,11-0.6,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 17 & as.numeric(exon_region[1]) <= 30){
					#segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*1,9-0.4,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*1,9+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*1,9-0.6,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 30 & as.numeric(exon_region[1]) <= 43){
					#segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*2,7-0.4,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*2,7+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*2,7-0.6,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 44 & as.numeric(exon_region[1]) <= 56){
					#segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*3,5-0.4,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*3,5+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*3,5-0.6,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 57 & as.numeric(exon_region[1]) <= 73){
					#segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*4,3-0.4,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*4,3+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*4,3-0.6,mutations[j],cex=1,col="red")
				}else if (as.numeric(exon_region[1]) >= 74 & as.numeric(exon_region[1]) <= 79){
					#segments(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*5,1-0.4,0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*5,1+0.3,col="black",lty=1,lwd=2)
					text(0.3*as.numeric(exon_region[1])+pos_region[1]/100-25*5,1-0.6,mutations[j],cex=1,col="red")
				}
			}
		}
	}
	
	#drawing annotations
	draw.annotation()
	
	dev.off() #close the graph file

}
