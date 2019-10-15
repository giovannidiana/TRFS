require(Hmisc)
NPAR=9
NMOD=512

tlist<-c("15","20","25");
flist<-c("sb","27","67","68","29","110");
parmat<-array(NA,dim=c(18,11,512));
parmat_av<-array(NA,dim=c(18,11));

likelihood_list<-matrix(NA,18,512)
for(t in 1:3){
    for(f in 1:6){
        tab<-read.table(paste("results/T",tlist[t],"_",flist[f],".dat",sep=""))
        names(tab)<-c("t1","t2","t3","t4","t5","t6","t7","t8","t9","model","chi2");
        for(k in 1:NMOD){
			tab2=tab[tab$model+1==k,]			
		    if(NROW(tab2)>0){ 
				ind=which.min(tab2$chi2)
				mat[f+(t-1)*6,k]=(9-tab2[ind,11])/9.;
				parmat[f+(t-1)*6,,k]=as.numeric(tab2[ind,]);
			}
        }
    }
}


