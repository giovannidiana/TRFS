
tlist<-c("15","20","25");
flist<-c("sb","27","67","68","29","110");

show_int <- function(f,t){
	layout(matrix(c(1,2),ncol=1))
	mat<-as.matrix(read.table(paste("intervals/interval_","T",tlist[t],"_",flist[f],".dat",sep="")));
	mat_sd<-as.matrix(read.table(paste("intervals/sd_","T",tlist[t],"_",flist[f],".dat",sep="")));
	adf<-barplot(mat[1,],names.arg=c("QL404","QL402","QL435"),main="ADF")
	arrows(adf,mat[1,]-mat_sd[1,],adf,mat[1,]+mat_sd[1,],code=3,length=0)

	asi<-barplot(mat[2,],names.arg=c("QL404","QL402","QL435"),main="ASI")
	arrows(asi,mat[2,]-mat_sd[2,],asi,mat[2,]+mat_sd[2,],code=3,length=0)
}

