# Pluritest - skript for use with the open Access tool for classification of pluripotent stem cells 

# Use arraydata points to file that can be read by lumiR (mod IGR)
# wd is the working directory where graphics and tables are stored.
# sample.names is the location of a file containing a vector containing the names of the samples for figure legend purposes. 

#args <- commandArgs(trailingOnly = TRUE)
load("~/Pluritest/workspace")

#sample.names = args[1]
#arraydata = args[2]

# load packages
	
	require(lumi)
	require(xtable)
	require(GO.db)
	
# set the stage
	#  samplenames <- readLines(file(sample.names))

	# #sink(file=paste(arraydata, "pluritest_log.txt", sep="_"))
	
 #    #IMPORT RAW DATA WITH LUMI
 #    working.lumi<-lumiR(arraydata, convertNuID = FALSE, annotationColumn="PROBE_ID")
	# sampleNames(working.lumi) <- samplenames
	
	# #identify probe ID column
 #    probeind <-match("PROBE_ID", colnames(fData(working.lumi)))
	# fData(working.lumi)[,probeind]<-gsub("\"","",fData(working.lumi)[,probeind])
	
	# #BOXPLOT OF RAW VALUES
	# pdf(file=paste(arraydata, "pluritest_image01.pdf", sep="_"))
	# plot(working.lumi, what='boxplot')
	# dev.off()

    #write("[pluritest_image01.pdf,Lumi QC Plot before robust spline normalization]",file="index.txt",append=TRUE)
    #write("Lumi QC Plot before robust spline normalization",file="index.txt",append=TRUE)   
    #write("",file="index.txt",append=TRUE)
	
	#CLUSTERING OF VARIANCE STABILIZED SAMPLES
	# working.lumi<-lumiT(working.lumi) #
	# hc<-hclust(as.dist(1-abs(cor(exprs(working.lumi[,])))))
	# pdf(file=paste(arraydata, "pluritest_image02a.pdf",sep="_"))
	# plot(hc,hang=-1,main="Clustering of vst-transformed samples",sub="distance based on pearson correlations",xlab='')
	# dev.off()
 #    #write("[pluritest_image02.pdf,Lumi QC Plot before robust spline normalization]",file="index.txt",append=TRUE)
    #write("Lumi QC Plot before robust spline normalization",file="index.txt",append=TRUE)   
    #write("",file="index.txt",append=TRUE)
    library(edgeR)
    datas = read.table("RNA-seq_pluri.txt", header = T)
    cpms = cpm(datas[,7:dim(datas)[2]], log = T)
    rownames(cpms) <- datas[,1]
	#RSN NORMALIZATION OF THE DATA
	A <- rownames(cpms)  #matches on ILMN_Ids for lumi/RSN
	B <- fData(H9targetArray)[,1] #for matches on ILMN_Ids for lumi/RSN
	sel.match <- match(B,A)
	#working.lumi <- working.lumi[sel.match[!is.na(sel.match)],] #subsets the exprSet user.lumi down to subsetuser.lumi to match H9 ref array rows
	#working.lumi<-lumiN(working.lumi,method="rsn",target=H9targetArray[is.na(sel.match)==FALSE,])

	#RSN NORMALIZATION OF THE DATA CONTINUED
	# assume that illumina Probe Ids are in fData[,1]
	#A <- fData(working.lumi)[,probeind] 
	#sel.match<-match(rownames(W15),A)
	
	#working.lumi.exprs <- cpms[sel.match[!is.na(sel.match)],]
    #rownames(working.lumi.exprs) <- A	
        
	#CALCULATION OF SCORES 
	try(
			{
				sel<-match(rownames(W15),A)
				coef<-c(-1.267095e+02  ,4.567437e-03 , 4.377068e-03  ,1.043193e-03)

                working.lumi.int  <- cpms[sel[!is.na(sel)],]
				H15.new<-predictH(working.lumi.int, W15[!is.na(sel),])
				H12.new<-predictH(working.lumi.int, W12[!is.na(sel),])

				rss.new<-apply((working.lumi.int - W12[!is.na(sel),]%*%H12.new)^2,2,sum)
				RMSE.new<-sqrt(rss.new/sum(!is.na(sel)))

				novel.new<-apply((working.lumi.int - W12[!is.na(sel),]%*%H12.new)^8,2,sum)
				# Novelty score - ideally less than 2 is good 
				novel.new<-(novel.new/sum(!is.na(sel)))^(1/8)

				# Pluripotency score - ideally over 20 is good
				s.new<-drop(coef[1] +coef[2:4]%*%H15.new[c(1,14,13),])
				print(s.new)
			}
	)
	
    # plot MULTICLASS PLURITEST & overview
	
	table.results<-matrix(,nrow=dim(cpms)[2],ncol=5)
	rownames(table.results)<-colnames(cpms)
	colnames(table.results)<-c("pluri-raw","pluri logit-p","novelty","novelty logit-p","RMSD")
	try(
        {
            print(s.new)
            pdf(file=paste(arraydata, "pluritest_pluripotency_score.pdf",sep="_"))
            par(mar=c(12,4,4,2))
            par(xaxt='n')
            plot(s.new,main="pluripotency",xlab="",ylab="pluripotency",ylim=c(-130,70))
            abline(h=25.25,lty="dashed",col="red")
            abline(h=59.95,lty="dashed",col="red")
            abline(h=-28.92,lty="dashed",col="lightblue")
            abline(h=-130,lty="dashed",col="lightblue")
            par(xaxt='s')
            axis(1,at=c(1:length(s.new)),labels=names(s.new),las=2)
            dev.off()
            
        }
	)

	table.results[,1]<-round(s.new,3)
	table.results[,2]<-round(exp(s.new)/(1+exp(s.new)),3)
	table.results[,3]<-round(novel.new,3)
	table.results[,5]<-round(RMSE.new,3)

	try(
        {
            pdf(file=paste(arraydata, "pluritest_clustering.pdf",sep="_"))
            color.palette = colorRampPalette(c("red","pink1","aliceblue","lightblue","blue"), bias=1)
            filled.contour2(y=c(-129:70),x=c((1:200)/50),background129_70x1_4,col=colram(50),nlevels=35,xlab="novelty",ylab="pluripotency")
            points(s.new~novel.new,cex=.4,main="Overview")
            
            dev.off()
        }			
			
    )
	
	try(
        {
            
        palette(colorRampPalette(c("green", "orange","orange","orange", "red"))(5))
        df.novelty.new<-data.frame(novelty= novel.new)
        pdf(file=paste(arraydata, "pluritest_novelty.pdf",sep="_"))
        par(mar=c(12,4,4,2))
        par(xaxt='n')
        barplot(novel.new,col=pmin(5,10*predict(logit.novelty,type="response",newdata= df.novelty.new)+1),names.arg=c(1:length(novel.new)),xlab="",xlim=c(0,length(novel.new)),width=.9,space=(1/9),ylim=c(0,4))
        
        title(main="novelty")
        par(xaxt='s')
        axis(1,at=c(1:nrow(	table.results))-.4,labels=names(s.new),las=2)
            
        table.results[,4]<-round(predict(logit.novelty,type="response",newdata= df.novelty.new),3)		
           
        dev.off()
            
		}
	)	

	# Save CSV FILE for TABLE
	table.results[,5]<-round(RMSE.new,3)
	write.csv(table.results,file=paste(arraydata, "pluritest.csv",sep="_"))


