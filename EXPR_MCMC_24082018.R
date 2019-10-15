
library(rjags)
library(mvtnorm)

args <- commandArgs(trailingOnly=TRUE)
chainused=5000
data<-read.table("20170524_FEDt.dat");
names(data)<-c("Batch","date","GT","day","food","T","ADF","ASI","NSM","time");

Diff <- function(x, start) as.numeric(x - as.Date(cut(start, "year")))
flist <- c(1, 2e+07, 6.3e+07, 6.3e+08,2e+09, 1.1e+10);
fcontrol=2e+09;
glist <- c("QL196","QL404","QL402","QL435");
collist <- c("black","blue","red","purple");
tcol<-c("#FFFF00","#FFC800","#FF9600");
neuron_lab <- c("ADF","ASI","NSM")
tlist <-c(15,20,25);
lut<-read.table("combined_batch_list_v2.0.dat")
#lut$V1 <- sapply(as.character(lut$V1),function(x) strsplit(x,"-")[[1]][2])
fcol <- rgb(0,seq(0,1,length=6),seq(1,0,length=6))
beads_red <- read.table("beads_red_13022018.dat")
beads_red$DATE=as.Date(beads_red$DATE,"%Y-%m-%d")
start_date=beads_red$DATE[1]
beads_red$DATE=Diff(beads_red$DATE,beads_red$DATE[1])
beads_red$NVAL=beads_red$NVAL/1e3
#names(beads_red)<-c("index","NVAL","BG","NBEADS","FOLDER","DATE");
beads_green <- read.table("beads_green_13022018.dat")
beads_green$DATE=as.Date(beads_green$DATE,"%Y-%m-%d")
beads_green$DATE=Diff(beads_green$DATE,beads_green$DATE[1])
beads_green$NVAL=beads_green$NVAL/1e3
#names(beads_green)<-c("index","NVAL","BG","NBEADS","FOLDER","DATE");

getdata <- function(Genotype){
		tmpset<-subset(data,GT==Genotype & day==6 & (Batch %in% lut$V1) & (T%in%tlist) & (food%in%flist),select=c("Batch","T","food","ADF","ASI","NSM"));
		g=match(Genotype,glist);

		tabtemp<-table(tmpset$Batch);
        batches<-names(tabtemp)[tabtemp>0]

		yenv=matrix(NA,ncol=2,nrow=length(batches))

		for(i in 1:length(batches)){
			tmpset$bid[tmpset$Batch==batches[i]]=i	
			t=tmpset[tmpset$Batch==batches[i],]$T[1]
			f=tmpset[tmpset$Batch==batches[i],]$food[1]
			yenv[i,] = c(i,6*(match(t,tlist)-1)+match(f,flist))
		}

		tmpset[,c("ADF","ASI","NSM")]=tmpset[,c("ADF","ASI","NSM")]/1e6
		tmpset$DATE=Diff(as.Date(strtrim(tmpset$Batch,8),"%Y%m%d"),start_date)
		meanmat=cbind(by(tmpset$ADF,yenv[tmpset$bid,2],mean),by(tmpset$ASI,yenv[tmpset$bid,2],mean),by(tmpset$NSM,yenv[tmpset$bid,2],mean))

		return(list(y=tmpset[,c("ADF","ASI","NSM")],
					ybid=tmpset$bid,
					ydate=tmpset$DATE,
					yenv=yenv,
					meanmat=meanmat,
					yC=data.frame(beads_red$NVAL,beads_green$NVAL,beads_red$DATE),
					NW=nrow(tmpset),
					NB=length(batches)))
} 

summary_data<-vector("list",4);
for(i in 1:4){
	tmp<-getdata(Genotype=glist[i])
	summary_data[[i]] = cbind(worms=table(tmp$yenv[tmp$ybid,2]),nbatches=table(tmp$yenv[,2]))
}

samples<-vector("list",8)
counter=0
for(i in 1:8) {
	counter=counter+1;
	load(paste("all_samples_28082018_ch",i,".RData",sep="")); 
	samples[[counter]]=all_samples;
}

get_expr_model7<-function(g,f,t,q=.01,ch=1){
	r<-matrix(NA,ncol=3,nrow=3)
	ind=6*(t-1)+f
	chain_len=length(samples[[ch]][[1]][["rat"]][1,1,,1])
	fromchain=(chain_len-chainused):chain_len
	r[1,1]<-mean(samples[[ch]][[g]][["rat"]][ind,1,fromchain,1])
	r[2,1]<-mean(samples[[ch]][[g]][["rat"]][ind,2,fromchain,1])
	r[3,1]<-mean(samples[[ch]][[g]][["rat"]][ind,3,fromchain,1])
	r[1,2:3]=quantile(samples[[ch]][[g]][["rat"]][ind,1,fromchain,1],probs=c(q,1-q))
	r[2,2:3]=quantile(samples[[ch]][[g]][["rat"]][ind,2,fromchain,1],probs=c(q,1-q))
	r[3,2:3]=quantile(samples[[ch]][[g]][["rat"]][ind,3,fromchain,1],probs=c(q,1-q))
		
	return(r)
}

get_cov_model7<-function(g,f,t,ch=1){
	cov_global<-rep(NA,6)
	cov_batch<-rep(NA,6)
	cov_worm<-rep(NA,6)
	cov_tot<-rep(NA,6)

	cov_global_sd<-rep(NA,6)
	cov_batch_sd<-rep(NA,6)
	cov_worm_sd<-rep(NA,6)
	cov_tot_sd<-rep(NA,6)

	ind=6*(t-1)+f

	chain_len=length(samples[[ch]][[1]][["rat"]][1,1,,1])
	fromchain=(chain_len-chainused+1):chain_len
	invT_global <- matrix(NA,ncol=6,nrow=chainused);
	invT_batch <- matrix(NA,ncol=6,nrow=chainused);
	invT_worm <- matrix(NA,ncol=6,nrow=chainused);
	invT_tot <- matrix(NA,ncol=6,nrow=chainused);
	counter=0;
	for(i in fromchain) {
		counter=counter+1;

		invT_global[counter,] = c(1/samples[[ch]][[g]][["TglobADF"]][ind,i,1],
		                          1/samples[[ch]][[g]][["TglobASI"]][ind,i,1],
		                          1/samples[[ch]][[g]][["TglobNSM"]][ind,i,1],
								  0,0,0)

		mat_tmp<-solve(samples[[ch]][[g]][["T.rat"]][,,ind,i,1])
		invT_batch[counter,1]=mat_tmp[1,1]
		invT_batch[counter,2]=mat_tmp[2,2]
		invT_batch[counter,3]=mat_tmp[3,3]
		invT_batch[counter,4]=mat_tmp[1,2]
		invT_batch[counter,5]=mat_tmp[1,3]
		invT_batch[counter,6]=mat_tmp[2,3]
		
		mat_tmp<-solve(samples[[ch]][[g]][["TE"]][,,ind,i,1])
		invT_worm[counter,1]=mat_tmp[1,1]
		invT_worm[counter,2]=mat_tmp[2,2]
		invT_worm[counter,3]=mat_tmp[3,3]
		invT_worm[counter,4]=mat_tmp[1,2]
		invT_worm[counter,5]=mat_tmp[1,3]
		invT_worm[counter,6]=mat_tmp[2,3]

		for(j in 1:6) invT_tot[counter,j]=invT_global[counter,j]+invT_batch[counter,j]+invT_worm[counter,j]

	}

    invT_tot[,1:3]<-sqrt(invT_tot[,1:3])

	cov_global=apply(invT_global,2,mean)
    cov_batch=apply(invT_batch,2,mean)
	cov_worm=apply(invT_worm,2,mean)
	cov_tot=apply(invT_tot,2,mean)
	
	cov_global_sd=apply(invT_global,2,quantile,probs=(c(0.05,.95)))
    cov_batch_sd=apply(invT_batch,2,quantile,probs=(c(0.05,.95)))
	cov_worm_sd=apply(invT_worm,2,quantile,probs=(c(0.05,.95)))
	cov_tot_sd=apply(invT_tot,2,quantile,probs=(c(0.05,.95)))
		
	return(list(cov_global=cov_global,
	            cov_batch=cov_batch,
				cov_worm=cov_worm ,
				cov_tot=cov_tot,
	            cov_global_sd=cov_global_sd, 
				cov_batch_sd=cov_batch_sd, 
				cov_worm_sd=cov_worm_sd,
				cov_tot_sd=cov_tot_sd))
}

decodeF_tph1<-function(g,t,ch=1){
	cov_tot<-vector("list",6)
	means <- vector("list",6)
	CM=matrix(0,6,6);
	dp=rep(NA,50);
	SetToDecode <- subset(data,GT==glist[g] & day==6 & (Batch %in% lut$V1) & (T==tlist[t]) & (food%in%flist),select=c("Batch","food","ADF","ASI","NSM"));
	DATE=Diff(as.Date(strtrim(SetToDecode$Batch,8),"%Y%m%d"),start_date)

	chain_len=length(samples[[ch]][[1]][["rat"]][1,1,,1])
	fromchain=rev((chain_len-50+1):chain_len)
	counter=0;
	for(i in fromchain) {
		counter=counter+1;
		image(CM,col=grey.colors(100))
		cat(counter,'        \r') 
		CMtmp=matrix(0,6,6);

		CC = samples[[ch]][[g]][["C"]][,,1,i,1]
		for(f in 1:6){
			ind=6*(t-1)+f

            means[[f]]=samples[[ch]][[g]][["rat"]][ind,2,i,1]

			cov_tot[[f]]=1/samples[[ch]][[g]][["TglobASI"]][ind,i,1]
			
			# Adding batch contribution
			mat_tmp<-solve(samples[[ch]][[g]][["T.rat"]][,,ind,i,1])[2,2]
			cov_tot[[f]]=cov_tot[[f]]+mat_tmp;
			# adding worm contribution	
			mat_tmp<-solve(samples[[ch]][[g]][["TE"]][,,ind,i,1])[2,2]
			cov_tot[[f]]=cov_tot[[f]]+mat_tmp;
		}

        
		for(w in 1:nrow(SetToDecode)){
			pval=rep(0,6);
			for(fprime in 1:6){
				#print(SetToDecode[w,2:4])
				#print(means[[fprime]])
				#print(cov_tot[[fprime]])
				ratio=c(CC[1,1] + CC[2,1]*DATE[w],
						CC[1,1] + CC[2,1]*DATE[w])
				pval[fprime] = dnorm(as.numeric(SetToDecode[w,4])/1e6/ratio,means[[fprime]],sqrt(cov_tot[[fprime]]),log=TRUE)
			}
			food_decoded=which.max(pval)
			food_orig=match(SetToDecode$food[w],flist)
#			cat(food_orig,' ',food_decoded,'\n')

			CMtmp[food_decoded,food_orig]=CMtmp[food_decoded,food_orig]+1
		}
		CM=CM+t(t(CMtmp)/colSums(CMtmp))
		CMtmp=t(t(CMtmp)/colSums(CMtmp))
		dp[counter]=sum(diag(CMtmp))/sum(CMtmp)


	}

	return(list(CM,dp))
}

decodeF_daf7<-function(g,t,ch=1){
	cov_tot<-vector("list",6)
	means <- vector("list",6)
	CM=matrix(0,6,6);
	dp=rep(NA,50);
	SetToDecode <- subset(data,GT==glist[g] & day==6 & (Batch %in% lut$V1) & (T==tlist[t]) & (food%in%flist),select=c("Batch","food","ADF","ASI","NSM"));
	DATE=Diff(as.Date(strtrim(SetToDecode$Batch,8),"%Y%m%d"),start_date)

	chain_len=length(samples[[ch]][[1]][["rat"]][1,1,,1])
	fromchain=rev((chain_len-50+1):chain_len)
	counter=0;
	for(i in fromchain) {
		counter=counter+1;
		image(CM,col=grey.colors(100))
		cat(counter,'        \r') 
		CMtmp=matrix(0,6,6);

		CC = samples[[ch]][[g]][["C"]][,,1,i,1]
		for(f in 1:6){
			ind=6*(t-1)+f

            means[[f]]=samples[[ch]][[g]][["rat"]][ind,c(1,3),i,1]

			cov_tot[[f]]=matrix( c(1/samples[[ch]][[g]][["TglobADF"]][ind,i,1],0,
		                           0,1/samples[[ch]][[g]][["TglobNSM"]][ind,i,1]),2,2);
			
			# Adding batch contribution
			mat_tmp<-solve(samples[[ch]][[g]][["T.rat"]][,,ind,i,1])[c(1,3),c(1,3)]
			cov_tot[[f]]=cov_tot[[f]]+mat_tmp;
			# adding worm contribution	
			mat_tmp<-solve(samples[[ch]][[g]][["TE"]][,,ind,i,1])[c(1,3),c(1,3)]
			cov_tot[[f]]=cov_tot[[f]]+mat_tmp;
		}

        
		for(w in 1:nrow(SetToDecode)){
			pval=rep(0,6);
			for(fprime in 1:6){
				#print(SetToDecode[w,2:4])
				#print(means[[fprime]])
				#print(cov_tot[[fprime]])
				ratio=c(CC[1,1] + CC[2,1]*DATE[w],
						CC[1,1] + CC[2,1]*DATE[w])
				pval[fprime] = dmvnorm(as.numeric(SetToDecode[w,c(3,5)])/1e6/ratio,means[[fprime]],cov_tot[[fprime]],log=TRUE)
			}
			food_decoded=which.max(pval)
			food_orig=match(SetToDecode$food[w],flist)
#			cat(food_orig,' ',food_decoded,'\n')

			CMtmp[food_decoded,food_orig]=CMtmp[food_decoded,food_orig]+1
		}
		CM=CM+t(t(CMtmp)/colSums(CMtmp))
		CMtmp=t(t(CMtmp)/colSums(CMtmp))
		dp[counter]=sum(diag(CMtmp))/sum(CMtmp)


	}

	return(list(CM,dp))
}

decodeF<-function(g,t,ch=1){
	cov_tot<-vector("list",6)
	means <- vector("list",6)
	CM=matrix(0,6,6);
	dp=rep(NA,50);
	SetToDecode <- subset(data,GT==glist[g] & day==6 & (Batch %in% lut$V1) & (T==tlist[t]) & (food%in%flist),select=c("Batch","food","ADF","ASI","NSM"));
	DATE=Diff(as.Date(strtrim(SetToDecode$Batch,8),"%Y%m%d"),start_date)

	chain_len=length(samples[[ch]][[1]][["rat"]][1,1,,1])
	fromchain=rev((chain_len-50+1):chain_len)
	counter=0;
	for(i in fromchain) {
		counter=counter+1;
		image(CM,col=grey.colors(100))
		cat(counter,'        \r') 
		CMtmp=matrix(0,6,6);

		CC = samples[[ch]][[g]][["C"]][,,1,i,1]
		for(f in 1:6){
			ind=6*(t-1)+f

            means[[f]]=samples[[ch]][[g]][["rat"]][ind,,i,1]

			cov_tot[[f]]=matrix( c(1/samples[[ch]][[g]][["TglobADF"]][ind,i,1],0,0,
		                           0, 1/samples[[ch]][[g]][["TglobASI"]][ind,i,1],0,
		                           0,0,1/samples[[ch]][[g]][["TglobNSM"]][ind,i,1]),3,3);
			
			# Adding batch contribution
			mat_tmp<-solve(samples[[ch]][[g]][["T.rat"]][,,ind,i,1])
			cov_tot[[f]]=cov_tot[[f]]+mat_tmp;
			# adding worm contribution	
			mat_tmp<-solve(samples[[ch]][[g]][["TE"]][,,ind,i,1])
			cov_tot[[f]]=cov_tot[[f]]+mat_tmp;
		}

        
		for(w in 1:nrow(SetToDecode)){
			pval=rep(0,6);
			for(fprime in 1:6){
				#print(SetToDecode[w,2:4])
				#print(means[[fprime]])
				#print(cov_tot[[fprime]])
				ratio=c(CC[1,1] + CC[2,1]*DATE[w],
				        CC[1,2] + CC[2,2]*DATE[w],
						CC[1,1] + CC[2,1]*DATE[w])
				pval[fprime] = dmvnorm(as.numeric(SetToDecode[w,3:5])/1e6/ratio,means[[fprime]],cov_tot[[fprime]],log=TRUE)
			}
			food_decoded=which.max(pval)
			food_orig=match(SetToDecode$food[w],flist)
#			cat(food_orig,' ',food_decoded,'\n')

			CMtmp[food_decoded,food_orig]=CMtmp[food_decoded,food_orig]+1
		}
		CM=CM+t(t(CMtmp)/colSums(CMtmp))
		CMtmp=t(t(CMtmp)/colSums(CMtmp))
		dp[counter]=sum(diag(CMtmp))/sum(CMtmp)


	}

	return(list(CM,dp))
}

decodeT<-function(g,f,ch=1){
	cov_tot<-vector("list",3)
	means <- vector("list",3)
	CM=matrix(0,3,3);
	SetToDecode <- subset(data,GT==glist[g] & day==6 & (Batch %in% lut$V1) & (T%in%tlist) & (food==flist[f]),select=c("Batch","T","ADF","ASI","NSM"));
	DATE=Diff(as.Date(strtrim(SetToDecode$Batch,8),"%Y%m%d"),start_date)

	chain_len=length(samples[[ch]][[1]][["rat"]][1,1,,1])
	fromchain=rev((chain_len-50+1):chain_len)
	counter=0;
	dp=rep(NA,50)
	for(i in fromchain) {
		counter=counter+1;
		image(CM,col=grey.colors(100))
		cat(counter,'        \r') 
		CMtmp=matrix(0,3,3);

		CC = samples[[ch]][[g]][["C"]][,,1,i,1]
		for(t in 1:3){
			ind=6*(t-1)+f

            means[[t]]=samples[[ch]][[g]][["rat"]][ind,,i,1]

			cov_tot[[t]]=matrix( c(1/samples[[ch]][[g]][["TglobADF"]][ind,i,1],0,0,
		                           0, 1/samples[[ch]][[g]][["TglobASI"]][ind,i,1],0,
		                           0,0,1/samples[[ch]][[g]][["TglobNSM"]][ind,i,1]),3,3);
			
			# Adding batch contribution
			mat_tmp<-solve(samples[[ch]][[g]][["T.rat"]][,,ind,i,1])
			cov_tot[[t]]=cov_tot[[t]]+mat_tmp;
			# adding worm contribution	
			mat_tmp<-solve(samples[[ch]][[g]][["TE"]][,,ind,i,1])
			cov_tot[[t]]=cov_tot[[t]]+mat_tmp;
		}

        
		for(w in 1:nrow(SetToDecode)){
			pval=rep(0,3);
			for(tprime in 1:3){
				#print(SetToDecode[w,2:4])
				#print(means[[fprime]])
				#print(cov_tot[[fprime]])
				ratio=c(CC[1,1] + CC[2,1]*DATE[w],
				        CC[1,2] + CC[2,2]*DATE[w],
						CC[1,1] + CC[2,1]*DATE[w])
				pval[tprime] = dmvnorm(as.numeric(SetToDecode[w,3:5])/1e6/ratio,means[[tprime]],cov_tot[[tprime]],log=TRUE)
			}
			t_decoded=which.max(pval)
			t_orig=match(SetToDecode$T[w],tlist)
#			cat(food_orig,' ',food_decoded,'\n')

			CMtmp[t_decoded,t_orig]=CMtmp[t_decoded,t_orig]+1
		}
		CM=CM+t(t(CMtmp)/colSums(CMtmp))
		CMtmp=t(t(CMtmp)/colSums(CMtmp))
		dp[counter]=sum(diag(CMtmp))/sum(CMtmp)

	}

	return(list(CM,dp))
}

gen_decoding_matrices<- function(){
	all_decoding_F<-vector("list",12)
	all_decoding_T<-vector("list",24)
	decoding_F_daf7<-vector("list",3)
	decoding_F_tph1<-vector("list",3)
	for(g in 1:4){
		for(t in 1:3){
			all_decoding_F[[3*(g-1)+t]]=decodeF(g,t)

		}
		for(f in 1:6){
			all_decoding_T[[6*(g-1)+f]]=decodeT(g,f)
		}
	}

	for(t in 1:3){
		decoding_F_daf7[[t]]=decodeF_daf7(3,t)
		decoding_F_tph1[[t]]=decodeF_tph1(2,t)
	}

	return(list(DF=all_decoding_F,DT=all_decoding_T,DF_tph1=decoding_F_tph1,DF_daf7=decoding_F_daf7))
}

#all_decodings=gen_decoding_matrices()

decoding_power.plot <- function(g,print=F){
	gcol=c("black","blue","red","purple")
	tmpF<-matrix(NA,50,12)
	tmpT<-matrix(NA,50,24)
	for(t in 1:3){ tmpF[,(g-1)*3+t]=all_decodings[[1]][[(g-1)*3+t]][[2]]}
	for(f in 1:6){ tmpT[,(g-1)*6+f]=all_decodings[[2]][[(g-1)*6+f]][[2]]}

    if(print) pdf(paste("results_expression_08082018/Decoding_Power_",glist[g],".pdf",sep=""),5,5)
	layout(matrix(1:2,1,2))
	boxplot(names=tlist,tmpF[,(g-1)*3 + 1:3],outline=F,col=gcol[g])
	boxplot(names=flist,tmpT[,(g-1)*6 + 1:6],las=2,outline=F,col=gcol[g])
	if(print) dev.off();
}

decoding_power_fn.plot <- function(print=F){
	gcol=c("black","blue","red","purple")
	tmpF<-array(NA,dim=c(50,3,3))
	for(t in 1:3){ 
		tmpF[,1,t]=all_decodings$DF[[(1-1)*3+t]][[2]]
		tmpF[,2,t]=all_decodings$DF_tph1[[t]][[2]]
		tmpF[,3,t]=all_decodings$DF_daf7[[t]][[2]]
	}

    if(print) pdf(paste("results_expression_04012019/Decoding_Power_fn.pdf",sep=""),8,5)
	layout(matrix(1:3,1,3))
	par(mar=c(10,5,3,2))
	boxplot(names=glist[1:3],tmpF[,1:3,1],outline=F,col=gcol[1:3],las=2,main=tlist[1],cex.main=2,cex.axis=2,cex.lab=2)
	boxplot(names=glist[1:3],tmpF[,1:3,2],outline=F,col=gcol[1:3],las=2,main=tlist[2],cex.main=2,cex.axis=2,cex.lab=2)
	boxplot(names=glist[1:3],tmpF[,1:3,3],outline=F,col=gcol[1:3],las=2,main=tlist[3],cex.axis=2,cex.main=2,cex.lab=2)
	if(print) dev.off();
}

decoding_power_fn.table <- function(print=F){
	tmpF<-array(NA,dim=c(50,3,3))
	for(t in 1:3){ 
		tmpF[,1,t]=all_decodings$DF[[(1-1)*3+t]][[2]]
		tmpF[,2,t]=all_decodings$DF_tph1[[t]][[2]]
		tmpF[,3,t]=all_decodings$DF_daf7[[t]][[2]]
	}
    
	u=matrix(NA,9,2)
	for(g in 1:3){
		for(t in 1:3){
			u[(g-1)*3+t,1]=mean(tmpF[,g,t]); 
			u[(g-1)*3+t,2]=sd(tmpF[,g,t]);
		}
	}

    if(print) write.table(u,file="results_expression_04012019/Decoding_Power_fn.txt",col.names=F,row.names=F)
	return(u)
}

get_FDR_model7<-function(g,t,ch=1){

	chain_len=length(samples[[ch]][[1]][["rat"]][1,1,,1])
	fromchain=(chain_len-chainused+1):chain_len
	dr <- matrix(NA,ncol=3,nrow=chainused);
	counter=0;
	for(i in fromchain) {
		counter=counter+1;
		mat_tmp=matrix(NA,nrow=3,ncol=6);
		for(f in 1:6){
			ind=6*(t-1)+f
			mat_tmp[,f]=samples[[ch]][[g]][["rat"]][ind,,i,1]
		}
		dr[counter,]=apply(mat_tmp,1,function(x) diff(range(x)))
	}
	dr_mean=apply(dr,2,mean)
	dr_sd=apply(dr,2,sd)
		
	return(list(dr_mean,dr_sd))
}

get_TDR_model7<-function(g,f,ch=1){

	chain_len=length(samples[[ch]][[1]][["rat"]][1,1,,1])
	fromchain=(chain_len-chainused+1):chain_len
	dr <- matrix(NA,ncol=3,nrow=chainused);
	counter=0;
	for(i in fromchain) {
		counter=counter+1;
		mat_tmp=matrix(NA,nrow=3,ncol=3);
		for(t in 1:3){
			ind=6*(t-1)+f
			mat_tmp[,t]=samples[[ch]][[g]][["rat"]][ind,,i,1]
		}
		dr[counter,]=apply(mat_tmp,1,function(x) diff(range(x)))
	}
	dr_mean=apply(dr,2,mean)
	dr_sd=apply(dr,2,sd)
		
	return(list(dr_mean,dr_sd))
}

get_expr_GT<-function(g,t,q=0.05,ch=1){
	resp_mean=matrix(NA,nrow=6,ncol=3)
	resp_q1=matrix(NA,nrow=6,ncol=3)
	resp_q2=matrix(NA,nrow=6,ncol=3)
	for(i in 1:6){
		resp=get_expr_model7(g,i,t,q,ch)
		resp_mean[i,]=resp[,1]
		resp_q1[i,]=resp[,2]
		resp_q2[i,]=resp[,3]
	}
	return(list(resp_mean,resp_q1,resp_q2))
}

get_expr_GF<-function(g,f,q=0.05,ch=1){
	resp_mean=matrix(NA,nrow=3,ncol=3)
	resp_q1=matrix(NA,nrow=3,ncol=3)
	resp_q2=matrix(NA,nrow=3,ncol=3)
	for(i in 1:3){
		resp=get_expr_model7(g,f,i,q,ch)
		resp_mean[i,]=resp[,1]
		resp_q1[i,]=resp[,2]
		resp_q2[i,]=resp[,3]
	}
	return(list(resp_mean,resp_q1,resp_q2))
}

get_cov_GT<-function(g,t,q=0.05,ch=1){
	cov_mean=matrix(NA,nrow=6,ncol=6)
	cov_low=matrix(NA,nrow=6,ncol=6)
	cov_hi=matrix(NA,nrow=6,ncol=6)
	for(i in 1:6){
		co=get_cov_model7(g,i,t,ch)
		cov_mean[i,]=co$cov_tot
		cov_low[i,]=co$cov_tot_sd[1,]
		cov_hi[i,]=co$cov_tot_sd[2,]
	}
	return(list(cov_mean,cov_low,cov_hi))
}

get_cov_GF<-function(g,f,q=0.05,ch=1){
	cov_mean=matrix(NA,nrow=3,ncol=6)
	cov_low=matrix(NA,nrow=6,ncol=6)
	cov_hi=matrix(NA,nrow=6,ncol=6)
	for(i in 1:3){
		co=get_cov_model7(g,f,i,ch)
		cov_mean[i,]=co$cov_tot
		cov_low[i,]=co$cov_tot_sd[1,]
		cov_hi[i,]=co$cov_tot_sd[2,]
	}
	return(list(cov_mean,cov_low,cov_hi))
}

get_expr_T.plot<-function(t,q=0.01,print=0,ch=1){
	gcol=c("black","blue","red","purple")
	custom_margins=c(2,2,2);
	if(print==1) pdf(paste("results_expression_02082018/T_",tlist[t],".pdf",sep=""),width=5,height=12)
	layout(matrix(1:4,ncol=1),heights=c(10,10,10,3))
	resp=vector("list",4)
	resp[[1]]<-get_expr_GT(1,t,q,ch);
	resp[[2]]<-get_expr_GT(2,t,q,ch);
	resp[[3]]<-get_expr_GT(3,t,q,ch);
	resp[[4]]<-get_expr_GT(4,t,q,ch);
	lims=matrix(NA,ncol=2,nrow=3)
	lims[1,] = c(2,50);
	lims[2,] = c(2.,50);
	lims[3,] = c(2,100);
	x_axis=flist
	
	for(i in c(3,1,2)){
	    par(mar=c(custom_margins[i],6,2,2),cex.axis=2,mgp=c(4,1.5,0))
		plot(x_axis,resp[[1]][[1]][,i],type='l',ylim=lims[i,],
			 xlab="",ylab=neuron_lab[i],cex=2,xaxt='n',cex.lab=2,lwd=1.5,log="x")
		rect(xleft=10^par("usr")[1]*1, ybottom=par("usr")[3]*1, 
		     xright=10^par("usr")[2],ytop=par("usr")[4]*1, 
			      lwd=3, border=tcol[t], xpd=TRUE)
		la=NA;
		axis(x_axis,at=x_axis,labels=la,las=2,lwd.ticks=1,lwd=0)
		arrows(x0=x_axis,x1=x_axis,y0=resp[[1]][[2]][,i],y1=resp[[1]][[3]][,i],length=0.05,angle=90,code=3)
		for(j in 2:4){
			lines(x_axis,resp[[j]][[1]][,i],col=gcol[j],lwd=1.5)
			#lines(x_axis,rep(1,length(flist)),lty=2)
			arrows(x0=x_axis,x1=x_axis,y0=resp[[j]][[2]][,i],y1=resp[[j]][[3]][,i],length=0.05,code=3,angle=90,lwd=1.5,col=gcol[j])
		}

	}
    
	plot(x_axis,rep(1,6),frame=F,xaxt='n',yaxt='n',type='n',ylab="",xlab="",log="x")
	axis(1,at=x_axis,labels=flist,las=2,lwd.ticks=1,lwd=0,line=-6,tick=F)


	if(print) dev.off()
}

get_expr_T.chaincomp<-function(g,t,q=0.01,print=0,chains=1:8){
	gcol=c("black","blue","red","purple")
	custom_margins=c(2,2,2);
	if(print==1) pdf(paste("results_expression_02082018/T_",tlist[t],".pdf",sep=""),width=5,height=12)
	layout(matrix(1:4,ncol=1),heights=c(10,10,10,3))
	resp=vector("list",length(chains))
	for(ch in 1:length(chains)) resp[[ch]] <- get_expr_GT(g,t,q,chains[ch]);
	lims=matrix(NA,ncol=2,nrow=3)
	lims[1,] = c(2,50);
	lims[2,] = c(2.,50);
	lims[3,] = c(40,80);
	x_axis=flist
	
	for(i in 1:3){
	    par(mar=c(custom_margins[i],6,2,2),cex.axis=2,mgp=c(4,1.5,0))
		plot(x_axis,resp[[1]][[1]][,i],type='l',ylim=lims[i,],
			 xlab="",ylab=neuron_lab[i],cex=2,xaxt='n',cex.lab=2,lwd=1.5,col=gcol[g],log="x")
		rect(xleft=10^par("usr")[1]*1, ybottom=par("usr")[3]*1, 
		     xright=10^par("usr")[2],ytop=par("usr")[4]*1, 
			      lwd=3, border=tcol[t], xpd=TRUE)
		la=NA;
		axis(1,x_axis,at=x_axis,labels=la,las=2,lwd.ticks=1,lwd=0)
		arrows(x0=x_axis,x1=x_axis,y0=resp[[1]][[2]][,i],y1=resp[[1]][[3]][,i],length=0.05,angle=90,code=3,col=gcol[g])
		for(j in 2:length(chains)){
			lines(x_axis,resp[[j]][[1]][,i],col=gcol[g],lwd=1.5)
			#lines(x_axis,rep(1,length(flist)),lty=2)
			arrows(x0=x_axis,x1=x_axis,y0=resp[[j]][[2]][,i],y1=resp[[j]][[3]][,i],length=0.05,code=3,angle=90,lwd=1.5,col=gcol[g])
		}

	}
    
	plot(x_axis,flist,frame=F,xaxt='n',yaxt='n',type='n',ylab="",xlab="",log="x")
	axis(1,at=x_axis,labels=flist,las=2,lwd.ticks=1,lwd=0,line=-6,tick=F)


	if(print) dev.off()
}


get_cov_T.plot<-function(t,print=0,ch=1){	
	
	gcol=c("black","blue","red","purple")
	custom_margins=c(2,2,2,2,2,2,2);
	if(print==1) pdf(paste("results_expression_02082018/noise_T",tlist[t],".pdf",sep=""),width=5,height=18)
	layout(matrix(1:4,ncol=1),heights=c(10,10,10,10))
	resp=vector("list",4)
	resp[[1]]<-get_cov_GT(1,t,ch);
	resp[[2]]<-get_cov_GT(2,t,ch);
	resp[[3]]<-get_cov_GT(3,t,ch);
	resp[[4]]<-get_cov_GT(4,t,ch);
	lims=matrix(NA,ncol=2,nrow=6)
	lims[1,] = c(0,20);
	lims[2,] = c(0,20);
	lims[3,] = c(0,20);
	x_axis=1:length(flist)
	
	for(i in 1:3){
	    par(mar=c(custom_margins[i],6,2,2),cex.axis=2,mgp=c(4,1.5,0))
		plot(x_axis,resp[[1]][[1]][,i],type='l',ylim=lims[i,],
			 xlab="",ylab="",cex=2,xaxt='n',cex.lab=2,lwd=1.5)
		rect(xleft=par("usr")[1]*1, ybottom=par("usr")[3]*1, 
		     xright=par("usr")[2],ytop=par("usr")[4]*1, 
			      lwd=3, border=tcol[t], xpd=TRUE)
		la=NA;
		axis(x_axis,at=x_axis,labels=la,las=2,lwd.ticks=1,lwd=0)
		arrows(x0=x_axis,x1=x_axis,y0=resp[[1]][[2]][,i],y1=resp[[1]][[3]][,i],length=0.05,angle=90,code=3)
		for(j in 2:4){
			lines(x_axis,resp[[j]][[1]][,i],col=gcol[j],lwd=1.5)
			arrows(x0=x_axis,x1=x_axis,y0=resp[[j]][[2]][,i],y1=resp[[j]][[3]][,i],length=0.05,code=3,angle=90,lwd=1.5,col=gcol[j])
		}

	}
    
	plot(x_axis,rep(1,6),frame=F,xaxt='n',yaxt='n',type='n',ylab="",xlab="")
	axis(1,at=x_axis,labels=flist,las=2,lwd.ticks=1,lwd=0,line=-6,tick=F)


	if(print) dev.off()
}

get_FDR_plot<-function(t,print=0,ch=1){
	
	gcol=c("black","blue","red","purple")
	custom_margins=c(2,2,2);
	if(print==1) pdf(paste("results_expression_02082018/EXPR_DR_T",tlist[t],".pdf",sep=""),width=5,height=12)
	layout(matrix(1:4,ncol=1),heights=c(10,10,10,3))
	co=vector("list",4)
	resp_mat=matrix(NA,ncol=4,nrow=3)
	sd_mat=matrix(NA,ncol=4,nrow=3)

	co[[1]]<-get_FDR_model7(1,t,ch);
	co[[2]]<-get_FDR_model7(2,t,ch);
	co[[3]]<-get_FDR_model7(3,t,ch);
	co[[4]]<-get_FDR_model7(4,t,ch);
	resp_mat[,1]<-co[[1]][[1]]
	resp_mat[,2]<-co[[2]][[1]]
	resp_mat[,3]<-co[[3]][[1]]
	resp_mat[,4]<-co[[4]][[1]]
	sd_mat[,1]<-co[[1]][[2]]
	sd_mat[,2]<-co[[2]][[2]]
	sd_mat[,3]<-co[[3]][[2]]
	sd_mat[,4]<-co[[4]][[2]]
	lims=matrix(NA,ncol=2,nrow=3)
	lims[1,] = c(0,30);
	lims[2,] = c(0,30);
	lims[3,] = c(0,50);
	x_axis=1:length(flist)
	
	for(i in 1:3){
	    par(mar=c(custom_margins[i],6,2,2),cex.axis=2,mgp=c(4,1.5,0))
		x_axis<-barplot(resp_mat[i,],col=gcol,ylim=lims[i,])
		rect(xleft=par("usr")[1]*1, ybottom=par("usr")[3]*1, 
		     xright=par("usr")[2],ytop=par("usr")[4]*1, 
			      lwd=3, border=tcol[t], xpd=TRUE)
		la=NA;
		arrows(x0=x_axis,x1=x_axis,y0=resp_mat[i,]-sd_mat[i,],y1=resp_mat[i,]+sd_mat[i,],length=0.05,angle=90,code=3)
			#arrows(x0=x_axis,x1=x_axis,y0=resp[[1]][[1]][,i]-resp[[j]][[2]][,i],y1=resp[[1]][[1]][,i]+resp[[j]][[2]][,i],length=0.05,code=3,angle=90,lwd=1.5,col=gcol[j])

	}
    
	plot(x_axis,rep(1,4),frame=F,xaxt='n',yaxt='n',type='n',ylab="",xlab="")
	axis(1,at=x_axis,labels=glist,las=2,lwd.ticks=1,lwd=0,line=-6,tick=F)


	if(print) dev.off()
}



get_TDR_plot<-function(f,print=0,ch=1){
	
	gcol=c("black","blue","red","purple")
	custom_margins=c(2,2,2);
	if(print==1) pdf(paste("results_expression_02082018/EXPR_DR_F",flist[f],".pdf",sep=""),width=5,height=12)
	layout(matrix(1:4,ncol=1),heights=c(10,10,10,3))
	co=vector("list",4)
	resp_mat=matrix(NA,ncol=4,nrow=3)
	sd_mat=matrix(NA,ncol=4,nrow=3)

	co[[1]]<-get_TDR_model7(1,f,ch);
	co[[2]]<-get_TDR_model7(2,f,ch);
	co[[3]]<-get_TDR_model7(3,f,ch);
	co[[4]]<-get_TDR_model7(4,f,ch);
	resp_mat[,1]<-co[[1]][[1]]
	resp_mat[,2]<-co[[2]][[1]]
	resp_mat[,3]<-co[[3]][[1]]
	resp_mat[,4]<-co[[4]][[1]]
	sd_mat[,1]<-co[[1]][[2]]
	sd_mat[,2]<-co[[2]][[2]]
	sd_mat[,3]<-co[[3]][[2]]
	sd_mat[,4]<-co[[4]][[2]]
	lims=matrix(NA,ncol=2,nrow=3)
	lims[1,] = c(0,30);
	lims[2,] = c(0,30);
	lims[3,] = c(0,70);
	x_axis=1:length(flist)
	
	for(i in 1:3){
	
	    par(mar=c(custom_margins[i],6,2,2),cex.axis=2,mgp=c(4,1.5,0))
		x_axis<-barplot(resp_mat[i,],col=gcol,ylim=lims[i,])
		if(i==1) mtext(side=3,f)
		rect(xleft=par("usr")[1]*1, ybottom=par("usr")[3]*1, 
		     xright=par("usr")[2],ytop=par("usr")[4]*1, 
			      lwd=3, border=fcol[f], xpd=TRUE)
		la=NA;
		arrows(x0=x_axis,x1=x_axis,y0=resp_mat[i,]-sd_mat[i,],y1=resp_mat[i,]+sd_mat[i,],length=0.05,angle=90,code=3)
			#arrows(x0=x_axis,x1=x_axis,y0=resp[[1]][[1]][,i]-resp[[j]][[2]][,i],y1=resp[[1]][[1]][,i]+resp[[j]][[2]][,i],length=0.05,code=3,angle=90,lwd=1.5,col=gcol[j])

	}
    
	x_axis<-barplot(rep(1,4),xaxt='n',yaxt='n',ylab="",xlab="",col="white",border=NA)
	axis(1,at=x_axis,labels=glist,las=2,lwd.ticks=1,lwd=0,line=-6,tick=F)

	if(print) dev.off()
}

get_TDR_plot2<-function(print=0,ch=1){
	
	gcol=c("black","blue","red","purple")
	custom_margins=c(2,2,2);
	if(print==1) pdf(paste("results_expression_28022018/EXPR_TDR2.pdf",sep=""),width=10,height=8)
	layout(matrix(1:12,ncol=4))
	co=vector("list",4*6)
	sd_mat=matrix(NA,ncol=4,nrow=3)

    for(g in 1:4){
		for(f in 1:6){
			co[[6*(g-1)+f]]<-get_TDR_model7(g,f,ch);
		}
	}

	lims=matrix(NA,ncol=2,nrow=3)
	lims[1,] = c(0,2);
	lims[2,] = c(0,2);
	lims[3,] = c(0,2);
	x_axis=1:length(flist)

    for(g in 1:4){
	for(i in 1:3){
		resp_mat=matrix(NA,ncol=6,nrow=3);
		sd_mat=matrix(NA,ncol=6,nrow=3);
		for(f in 1:6) {
			resp_mat[,f]=co[[6*(g-1)+f]][[1]]
			sd_mat[,f]=co[[6*(g-1)+f]][[2]]
		}
		x_axis<-barplot(resp_mat[i,],ylim=lims[i,],main=glist[g],cex.main=2,col=fcol)
		#if(i==1) mtext(side=3,f)
		#rect(xleft=par("usr")[1]*1, ybottom=par("usr")[3]*1, 
		#     xright=par("usr")[2],ytop=par("usr")[4]*1, 
	#		      lwd=3, border=fcol[f], xpd=TRUE)
		la=NA;
		arrows(x0=x_axis,x1=x_axis,y0=resp_mat[i,]-sd_mat[i,],y1=resp_mat[i,]+sd_mat[i,],length=0.05,angle=90,code=3)
			#arrows(x0=x_axis,x1=x_axis,y0=resp[[1]][[1]][,i]-resp[[j]][[2]][,i],y1=resp[[1]][[1]][,i]+resp[[j]][[2]][,i],length=0.05,code=3,angle=90,lwd=1.5,col=gcol[j])

	}}
    
	#x_axis<-barplot(rep(1,4),xaxt='n',yaxt='n',ylab="",xlab="",col="white",border=NA)
	#axis(1,at=x_axis,labels=glist,las=2,lwd.ticks=1,lwd=0,line=-6,tick=F)

	if(print) dev.off()
}

chain.plot <- function(g=4,neur=2, ind=5,chin=1:8){
	plot(samples[[chin[1]]][[g]]$rat[ind,neur,,1],type='l');
	for(i in chin) lines(samples[[i]][[g]]$rat[ind,neur,,1],col=rainbow(8)[i])
	#plot(samples[[chin[1]]][[g]]$TE[1,1,ind,,1],type='l');
	#for(i in chin) lines(samples[[i]][[g]]$TE[1,1,ind,,1],col=rainbow(8)[i])
}

chain.plot.popvariance <- function(g=4,neur=1, ind=5,chin=1:8){
	chain_len=length(samples[[1]][[1]][["rat"]][1,1,,1])
	tmp=vector("numeric",chain_len);
	tmp_mat=matrix(0,3,3)
	for(u in 1:chain_len) tmp[u] = solve(samples[[chin[1]]][[g]]$TE[,,ind,u,1])[c(1,5,9)][neur]
	plot(tmp,type='l');
	for(i in chin){
		for(u in 1:chain_len) tmp[u] = solve(samples[[chin[i]]][[g]]$TE[,,ind,u,1])[c(1,5,9)][neur]		
		lines(tmp,col=rainbow(8)[i])
	}
}
plot_bead_fit<- function(q=0.05,ch=1,print=0){
	chain_len=length(samples[[ch]][[1]][["rat"]][1,1,,1])
	fromchain=(chain_len-chainused+1):chain_len
	if(print==1) pdf("results_expression_28022018/beads_trends.pdf",12,6);
	layout(matrix(1:2,ncol=2))
	par(mar=c(7,7,2,2),cex.lab=2)
	r<-seq(30,3000,5);
	rmat_red=matrix(NA,ncol=length(r),nrow=chainused)
	rmat_green=matrix(NA,ncol=length(r),nrow=chainused)
	counter=0;
	for (i in fromchain){
		counter=counter+1
		rmat_red[counter,]<-samples[[ch]][[1]]$C[1,1,i,ch]+r*samples[[ch]][[1]]$C[2,1,i,ch]
		rmat_green[counter,]<-samples[[ch]][[1]]$C[1,2,i,ch]+r*samples[[ch]][[1]]$C[2,2,i,ch]
	}
	pr<-apply(rmat_red,2,function(x) quantile(x,probs=c(q,1-q)))
	pg<-apply(rmat_green,2,function(x) quantile(x,probs=c(q,1-q)))
    
	plot(NVAL~DATE,beads_red)
	polygon(c(r,rev(r)),c(pr[1,],rev(pr[2,])),col=rgb(1,0,0,.3))
	plot(NVAL~DATE,beads_green)
	polygon(c(r,rev(r)),c(pg[1,],rev(pg[2,])),col=rgb(1,1,0,.3))
	if(print==1) dev.off()
}

plot_bead_fit2<- function(q=0.05,ch=1,print=0,chainused=500){
	u<-gensample_beads()
	chain_len=length(u$C[1,1,,1])
	fromchain=(chain_len-chainused+1):chain_len
	if(print==1) pdf("results_expression_28022018/beads_trends.pdf",12,6);
	layout(matrix(1:2,ncol=2))
	par(mar=c(7,7,2,2),cex.lab=2)
	r<-seq(30,3000,5);
	rmat_red=matrix(NA,ncol=length(r),nrow=chainused)
	rmat_green=matrix(NA,ncol=length(r),nrow=chainused)
	counter=0;
	for (i in fromchain){
		counter=counter+1
		rmat_red[counter,]<-u$C[1,1,i,ch]+r*u$C[2,1,i,ch]
		rmat_green[counter,]<-u$C[1,2,i,ch]+r*u$C[2,2,i,ch]
	}
	pr<-apply(rmat_red,2,function(x) quantile(x,probs=c(q,1-q)))
	pg<-apply(rmat_green,2,function(x) quantile(x,probs=c(q,1-q)))
    
	plot(NVAL~DATE,beads_red)
	polygon(c(r,rev(r)),c(pr[1,],rev(pr[2,])),col=rgb(1,0,0,.3))
	plot(NVAL~DATE,beads_green)
	polygon(c(r,rev(r)),c(pg[1,],rev(pg[2,])),col=rgb(1,1,0,.3))
	if(print==1) dev.off()
}





