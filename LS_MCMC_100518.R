require( survival )
library(rjags)
# Load and adjust the dataset
data <- read.table("20170911_lifespan_curated.csv");
# Give headers to each field. Notice that Maureen dataset is different from Eugeni and Dhaval's.
colnames(data) <- c("Batch","day0","GT","Strain","T2","Food","Age","Freq","Censor","Cause")
# Replicate rows according to the frequency field.
data <- data[rep( 1:nrow( data ), data$Freq ),]
# Remove "Cause" and "Freq" columns from the dataset since we don't need those.
data <- data[,c("Strain","Food","T2","Batch","Age","Censor")]
# Find the list of strains
StrainList <- as.character(names(table(data$Strain)));

# List of food levels to plot
tlist=c(15,17.5,20,25)
flist=c(1,2e7,6.32e7,6.32e8,2e9,1.12e10)
gcol=c("black","blue","red","purple")
gcola=c(rgb(0,0,0,0.2),rgb(0,0,1,.2),rgb(1,0,0,.2),rgb(0.62,0.12,0.94,.2))
tcol=c("blue","green","orange","red")
tcola=c(rgb(0,0,1,.2),rgb(0,1,0,.2),rgb(1,0.65,0,.2),rgb(1,0,0,.2))

GetData <- function(Gset="N2",Fset=1,Tset=20,Use_Previous_Data=F){
	data <- read.table("20170911_lifespan_curated.csv");
	colnames(data) <- c("Batch","day0","GT","Strain","T2","Food","Age","Freq","Censor","Cause");
	data <- data[,c("Strain","Food","T2","Batch","Age","Freq","Censor")];
    data <- droplevels(data[data$Strain==Gset & data$T2==Tset & data$Food==Fset,]);
	tb<-table(data$Batch);
	batch_list<-names(tb[tb>0]);
	NB<-length(batch_list);
	datalist<-vector("list",NB+1);

    Nworm=by(data$Freq,data$Batch,function(x) sum(x))[batch_list]
    N=by(data$Age,data$Batch,max)[batch_list]
	NsurvNB=vector("list",NB);
   
    Ncumul=0;
	for(i in 1:NB){
		NsurvNB[[i]]<-subset(data,Batch==batch_list[i], select=c("Batch", "Age", "Freq","Censor"))
		NsurvNB[[i]]$BatchID=i		
	}

	NsurvCombined=do.call(rbind,NsurvNB);
	NsurvCombined=NsurvCombined[rep(1:nrow(NsurvCombined),NsurvCombined$Freq),]

	return(NsurvCombined);
}

do_rjags <- function(d,B=NA){
    model="model {

	    for(i in 1:N){
			censor[i] ~ dinterval( t[i] , cutoff[i])
			t[i] ~ dweib(nu[g[i]],scaleb[g[i]])
		}
    # dgamma(shape,rate) 
	# nu in [4,20]
	# scale in [40,80]
	    for(i in 1:K){
			nu[i] ~ dnorm(nu0,tnu0)T(0,)
			scale[i] ~ dnorm(scale0,tscale0)T(0,)
			scaleb[i] = scale[i]^(-nu[i])
		}

		nu0 ~ dgamma(2,.1)
		tnu0 ~ dgamma(1,.02)
		scale0 ~ dgamma(2,.01)
		tscale0 ~ dgamma(1,02)

         
	}"
    
	surv_object=survfit( Surv(Age,1-Censor) ~ Batch ,data=d);
	tab1 <- summary( surv_object,
	                 rmean="common")$table;

	d$cutoff=d$Age
	d$Age[d$Censor==1]=NA
	data=list(t=d$Age,cutoff=d$cutoff,censor=d$Censor,N=nrow(d),g=d$BatchID,K=max(d$BatchID))
	
	varnames=c("nu","scaleb","scale","nu0","tnu0","scale0","tscale0")
	burn_in=1000;
	steps=10000;
	thin=1;
	
	fileConn=file("model.tmp")
	writeLines(model,fileConn);
	close(fileConn)
	
	m=jags.model(file="model.tmp",data=data);
	update(m,burn_in)
	draw=jags.samples(m,steps,thin=thin,variable.names=varnames)

	return(list(draw=draw,surv=surv_object))
}

load("all_LS_samples_130518.RData")
#all_samples<-vector("list",24*4)
#for(g in 1:4){
#	for(f in 1:6){
#		for(t in 1:4){
#			all_samples[[24*(g-1)+6*(t-1)+f]]=do_rjags(GetData(glist[g],flist[f],tlist[t]))
#		}
#	}
#}	

GetLS <- function(g, f, t){
	ind=24*(g-1)+6*(t-1)+f
	draw=all_samples[[ind]]$draw
	mean_lifespan = mean(draw$scale0[1,,1]*gamma(1+1/draw$nu0[1,,1]))
	sd_lifespan = sd(draw$scale0[1,,1]*gamma(1+1/draw$nu0[1,,1]))
	return(c(mean_lifespan,sd_lifespan))

}

summary_data <- function(g,t){
	summ<-matrix(NA,6,5)
	colnames(summ)<-c("died","censored","trials","mean LS","sd LS")
	for(f in 1:6){
		dt <- GetData(glist[g],flist[f],tlist[t])
		summ[f,1] = sum(dt$Censor==0)
		summ[f,2] = sum(dt$Censor==1)
		summ[f,3] = max(dt$BatchID)
		summ[f,4:5] = GetLS(g,f,t)
	}
	return(summ)
}


LS_Fresp.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/LS_Fresp.pdf",sep=""),width=15,height=4)
	layout(matrix(1:4,nrow=1))
	par(mar=c(8,3,4,1))
	for(t in 1:4){
		plot(flist,rep(0,6),type='n',xlab="",ylab="",ylim=c(15,66),cex=1.5,cex.axis=1.5,log="x",xaxt='n',main=paste("T =",tlist[t]))
		axis(1,at=flist,labels=flist,cex.axis=1.5,las=2)
		for(g in 1:4){
			resp=matrix(NA,ncol=2,nrow=6);
			for(f in 1:6){
				ind=24*(g-1)+6*(t-1)+f
				resp[f,]=GetLS(g,f,t)
				
			}
			lines(flist,resp[,1],col=gcol[g])
			arrows(x0=flist,x1=flist,y0=resp[,1]-resp[,2],y1=resp[,1]+resp[,2],code=3,length=0,col=gcol[g])

		}
	}
	if(print) dev.off()
}

LS_Fresp_byT.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/LS_Fresp_byT.pdf",sep=""),width=15,height=4)
	layout(matrix(1:4,nrow=1))
	par(mar=c(8,3,4,1))
	for(g in 1:4){
		plot(flist,rep(0,6),type='n',xlab="",ylab="",ylim=c(15,66),cex=1.5,cex.axis=1.5,log="x",xaxt='n',main=paste("g =",glist[g]))
		axis(1,at=flist,labels=flist,cex.axis=1.5,las=2)
		for(t in 1:4){
			resp=matrix(NA,ncol=2,nrow=6);
			for(f in 1:6){
				ind=24*(g-1)+6*(t-1)+f
				resp[f,]=GetLS(g,f,t)
				
			}
			lines(flist,resp[,1],col=tcol[t])
			arrows(x0=flist,x1=flist,y0=resp[,1]-resp[,2],y1=resp[,1]+resp[,2],code=3,length=0,col=tcol[t])

		}
	}
	if(print) dev.off()
}

LS_Tresp.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/LS_Tresp.pdf",sep=""),width=15,height=10)
	layout(matrix(1:6,nrow=2))
	par(mar=c(8,3,4,1))
	for(f in 1:6){
		plot(tlist,rep(0,4),type='n',xlab="",ylab="",ylim=c(15,66),cex=1.5,cex.axis=1.5,xaxt='n',main=paste("T =",flist[f]))
		axis(1,at=tlist,labels=tlist,cex.axis=1.5,las=2)
		for(g in 1:4){
			resp=matrix(NA,ncol=2,nrow=4);
			for(t in 1:4){
				ind=24*(g-1)+6*(t-1)+f
				resp[t,]=GetLS(g,f,t)
				
			}
			lines(tlist,resp[,1],col=gcol[g])
			arrows(x0=tlist,x1=tlist,y0=resp[,1]-resp[,2],y1=resp[,1]+resp[,2],code=3,length=0,col=gcol[g])

		}
	}
	if(print) dev.off()
}

LS_Tresp_byF.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/LS_Tresp_byF.pdf",sep=""),width=15,height=5)
	layout(matrix(1:4,nrow=1))
	par(mar=c(8,3,4,1))
	for(g in 1:4){
		plot(tlist,rep(0,4),type='n',xlab="",ylab="",ylim=c(15,66),cex=1.5,cex.axis=1.5,xaxt='n',main=paste("G =",glist[g]))
		axis(1,at=tlist,labels=tlist,cex.axis=1.5,las=2)
		for(f in 1:6){
			resp=matrix(NA,ncol=2,nrow=4);
			for(t in 1:4){
				ind=24*(g-1)+6*(t-1)+f
				resp[t,]=GetLS(g,f,t)
				
			}
			lines(tlist,resp[,1],col=gray.colors(6)[f])
			arrows(x0=tlist,x1=tlist,y0=resp[,1]-resp[,2],y1=resp[,1]+resp[,2],code=3,length=0,col=gray.colors(6)[f])

		}
	}
	if(print) dev.off()
}

# Plot survival curves across genotypes
Surv_Gcomp.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/Surv_Gcomp.pdf",sep=""),width=10,height=15)
	layout(matrix(1:24,ncol=4,nrow=6))
	par(mar=c(2,3,1,1))

    days<-rep(NA,6*4)

	for(t in 1:4){
		for(g in 1:4){
			for(f in 1:6){
				ind=24*(g-1)+6*(t-1)+f
				days[6*(g-1)+f]=max(summary(all_samples[[ind]]$surv)$time)
			}
		}
		
		maxday=max(days)
		time_axis=seq(0,maxday,1)
		meansurv=rep(NA,length(time_axis))
		surv_band_up=rep(NA,length(time_axis))
		surv_band_down=rep(NA,length(time_axis))
		Nsamples=length(all_samples[[1]]$draw$nu0[1,,1])

		for(f in 1:6){
			for(g in 1:4){
				ind=24*(g-1)+6*(t-1)+f
				add=F; if(g>1) add=T;
				traj=matrix(NA,ncol=length(time_axis),nrow=100)
				draw=all_samples[[ind]]$draw
				for(k in 1:100){
					traj[k,]=exp(-(time_axis/draw$scale0[1,Nsamples-k,1])^draw$nu0[1,Nsamples-k,1])
				}

				band=apply(traj,2,function(x) quantile(x,probs = c(0.2,0.8)))

				if(add) {
					polygon(c(time_axis,rev(time_axis)),c(band[1,],rev(band[2,])),col=gcola[g])
				} else {
					plot(time_axis,band[1,],col=gcol[g],type='l')
					polygon(c(time_axis,rev(time_axis)),c(band[1,],rev(band[2,])),col=gcola[g])
				}

			}
		}
	}
	if(print) dev.off()
}

Surv_Tcomp.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/Surv_Tcomp.pdf",sep=""),width=13,height=15)
	layout(matrix(1:24,ncol=4,nrow=6))
	par(mar=c(2,3,1,1))
	
	days<-rep(NA,96)
	for(g in 1:4){
		for(f in 1:6){
			for(t in 1:4){
				ind=24*(g-1)+6*(t-1)+f
				days[ind]=max(summary(all_samples[[ind]]$surv)$time)
			}
		}
	}
	
	maxday=max(days)
	time_axis=seq(0,maxday,1)
	meansurv=rep(NA,length(time_axis))
	surv_band_up=rep(NA,length(time_axis))
	surv_band_down=rep(NA,length(time_axis))
	Nsamples=length(all_samples[[1]]$draw$nu0[1,,1])

	for(g in 1:4){
		for(f in 1:6){
			for(t in 1:4){
				ind=24*(g-1)+6*(t-1)+f	
				add=F; if(t>1) add=T;
				traj=matrix(NA,ncol=length(time_axis),nrow=100)
				draw=all_samples[[ind]]$draw
				for(k in 1:100){
					traj[k,]=exp(-(time_axis/draw$scale0[1,Nsamples-k,1])^draw$nu0[1,Nsamples-k,1])
				}

				band=apply(traj,2,function(x) quantile(x,probs = c(0.2,0.8)))

				if(add) {
					polygon(c(time_axis,rev(time_axis)),c(band[1,],rev(band[2,])),col=tcola[t])
				} else {
					plot(time_axis,band[1,],col=gcol[t],type='l')
					polygon(c(time_axis,rev(time_axis)),c(band[1,],rev(band[2,])),col=tcola[t])
				}

			}
		}
	}
	if(print) dev.off()
}

LS_FDR.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/LS_FDR.pdf",sep=""),width=12,height=4)
	layout(matrix(1:4,nrow=1))
	par(mar=c(4,6,3,3))
	Nsamples=length(all_samples[[1]]$draw$nu0[1,,1])
	for(g in 1:4){
		dr=matrix(NA,ncol=4,nrow=2000);
	    colnames(dr)=tlist
		for(t in 1:4){
			resp=matrix(NA,ncol=6,nrow=2000);
			for(f in 1:6){
				ind=24*(g-1)+6*(t-1)+f
				draw=all_samples[[ind]]$draw
				resp[,f]=draw$scale0[1,(Nsamples-2000+1):Nsamples,1]*gamma(1+1/draw$nu0[1,(Nsamples-2000+1):Nsamples,1])
			}
			dr[,t]=apply(resp,1,function(x) diff(range(x)))

		}
		ylab="";if(g==1) ylab="Dynamic range";
			boxplot(dr,outline=F,cex.axis=1.5,ylim=c(0,40),main=glist[g],
					cex.main=1.5,col=gcol[g],cex.axis=1.5,cex.lab=1.5,ylab=ylab)
	}
	if(print) dev.off()
}

LS_TDR.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/LS_TDR.pdf",sep=""),width=12,height=4)
	layout(matrix(1:4,nrow=1))
	par(mar=c(4,6,3,3))
	Nsamples=length(all_samples[[1]]$draw$nu0[1,,1])
	for(g in 1:4){
		dr=matrix(NA,ncol=6,nrow=2000);
	    colnames(dr)=flist
		for(f in 1:6){
			resp=matrix(NA,ncol=4,nrow=2000);
			for(t in 1:4){
				ind=24*(g-1)+6*(t-1)+f
				draw=all_samples[[ind]]$draw
				resp[,t]=draw$scale0[1,(Nsamples-2000+1):Nsamples,1]*gamma(1+1/draw$nu0[1,(Nsamples-2000+1):Nsamples,1])
			}
			dr[,f]=apply(resp,1,function(x) diff(range(x)))

		}
		ylab="";if(g==1) ylab="T Dynamic range";
			boxplot(dr,outline=F,cex.axis=1.5,ylim=c(5,50),main=glist[g],
					cex.main=1.5,col=gcol[g],cex.axis=1.5,cex.lab=1.5,ylab=ylab)
	}
	if(print) dev.off()
}

# NEED TO ADD SNR FOR LIFESPAN... or maybe decoding is better.

GetSurv <- function(Strain, Food, Temp){
	temp <- data[data$Strain==Strain & data$Food==Food & data$T2==Temp,];
	return( survfit( Surv(Age,1-Censor) ~ Batch,data=temp) )
	
}

analyse <- function(d){
	res=do_rjags(d);
	s<-seq(0,max(summary(res$surv)$time),1);
	plot(res$surv,col="red")
	inds=sample(2500:3000,size=20); for(i in inds){a=res$draw$nu[1,i,1]; b=res$draw$scaleb[1,i,1];lines(s,exp(-s^a*b))}
}

compare <- function(res,inds){
	s<-seq(0,max(summary(res$surv)$time),1);
	plot(res$surv,col="red");
	for(i in inds){a=res$draw$nu[1,i,1]; b=res$draw$scaleb[1,i,1];lines(s,exp(-s^a*b))}
}

Decoding_F <- function(Gind,Tind){
	compvec=rep(NA,6)
	confmat=matrix(0,ncol=6,nrow=6)
	confmat_tmp=matrix(0,ncol=6,nrow=6)
	subdata<-subset(data,Strain==glist[Gind] & T2==tlist[Tind] & Censor==0 & (Food %in% flist))
	Nsamples=length(all_samples[[1]]$draw$nu0[1,,1])
	sample_used=200
	sampled_MLS=matrix(NA,ncol=6,nrow=sample_used);
	sampled_SDLS=matrix(NA,ncol=6,nrow=sample_used);
	sampled_BVLS=matrix(NA,ncol=6,nrow=sample_used);
	dec_pow=rep(NA,sample_used)
	
	for(f in 1:6){
		draw=all_samples[[(Gind-1)*24+(Tind-1)*6+f]]$draw
		sampled_MLS[,f]=draw$scale0[1,(Nsamples-sample_used+1):Nsamples,1]*gamma(1+1/draw$nu0[1,(Nsamples-sample_used+1):Nsamples,1])
		sampled_SDLS[,f]=draw$scale0[1,(Nsamples-sample_used+1):Nsamples,1]^2*(gamma(1+2/draw$nu0[1,(Nsamples-sample_used+1):Nsamples,1])-gamma(1+1/draw$nu0[1,(Nsamples-sample_used+1):Nsamples,1])^2)
		sampled_BVLS[,f]=1/draw$tnu0[1,(Nsamples-sample_used+1):Nsamples,1]
	}

	for(k in 1:sample_used){
		confmat_tmp=matrix(0,6,6);
		for(worm in 1:nrow(subdata)){
			for(f in 1:6){
				compvec[f]=dnorm(subdata$Age[worm],sampled_MLS[k,f],sqrt(sampled_SDLS[k,f]+sampled_BVLS[k,f]),log=T)
			}
			dec_f=which.max(compvec)
			f_orig=which(subdata$Food[worm]==flist)
			
			confmat_tmp[dec_f,f_orig]=confmat_tmp[dec_f,f_orig]+1;
			confmat[dec_f,f_orig]=confmat[dec_f,f_orig]+1;
		}
		dec_pow[k]=sum(diag(confmat_tmp))/sum(confmat_tmp)
	}

	return(list(CM=confmat,DP=dec_pow))
}	

Decoding_T <- function(Gind,Find){
	compvec=rep(NA,4)
	confmat=matrix(0,ncol=4,nrow=4)
	confmat_tmp=matrix(0,ncol=4,nrow=4)
	subdata<-subset(data,Strain==glist[Gind] & Food==flist[Find] & Censor==0 & (T2 %in% tlist))
	Nsamples=length(all_samples[[1]]$draw$nu0[1,,1])
	sample_used=200
	sampled_MLS=matrix(NA,ncol=4,nrow=sample_used);
	sampled_SDLS=matrix(NA,ncol=4,nrow=sample_used);
	sampled_BVLS=matrix(NA,ncol=4,nrow=sample_used);
	dec_pow=rep(NA,sample_used)
	
	for(t in 1:4){
		draw=all_samples[[(Gind-1)*24+(t-1)*6+Find]]$draw
		sampled_MLS[,t]=draw$scale0[1,(Nsamples-sample_used+1):Nsamples,1]*gamma(1+1/draw$nu0[1,(Nsamples-sample_used+1):Nsamples,1])
		sampled_SDLS[,t]=draw$scale0[1,(Nsamples-sample_used+1):Nsamples,1]^2*(gamma(1+2/draw$nu0[1,(Nsamples-sample_used+1):Nsamples,1])-gamma(1+1/draw$nu0[1,(Nsamples-sample_used+1):Nsamples,1])^2)
		sampled_BVLS[,t]=1/draw$tnu0[1,(Nsamples-sample_used+1):Nsamples,1]
	}

	for(k in 1:sample_used){
		confmat_tmp=matrix(0,4,4)
		for(worm in 1:nrow(subdata)){
			for(t in 1:4){
				compvec[t]=dnorm(subdata$Age[worm],sampled_MLS[k,t],sqrt(sampled_SDLS[k,t]+sampled_BVLS[k,t]),log=T)
			}
			dec_t=which.max(compvec)
			t_orig=which(subdata$T2[worm]==tlist)
			confmat[dec_t,t_orig]=confmat[dec_t,t_orig]+1;
			confmat_tmp[dec_t,t_orig]=confmat[dec_t,t_orig]+1;
		}
		dec_pow[k]=sum(diag(confmat_tmp))/sum(confmat_tmp)
		
	}

	return(list(CM=confmat,DP=dec_pow))
}

Decoding_FT <- function(Gind){
	compvec=rep(NA,4)
	confmat=matrix(0,ncol=24,nrow=24)
	confmat_tmp=matrix(0,ncol=24,nrow=24)
	subdata<-subset(data,Strain==glist[Gind] & (Food %in% flist) & Censor==0 & (T2 %in% tlist))
	Nsamples=length(all_samples[[1]]$draw$nu0[1,,1])
	sample_used=200
	sampled_MLS=matrix(NA,ncol=24,nrow=sample_used);
	sampled_SDLS=matrix(NA,ncol=24,nrow=sample_used);
	sampled_BVLS=matrix(NA,ncol=24,nrow=sample_used);
	dec_pow=rep(0,sample_used)
	
	for(t in 1:4){
		for(f in 1:6){
			matrix_index=(t-1)*6+f
			draw=all_samples[[(Gind-1)*24+(t-1)*6+f]]$draw
			sampled_MLS[,matrix_index]=draw$scale0[1,(Nsamples-sample_used+1):Nsamples,1]*gamma(1+1/draw$nu0[1,(Nsamples-sample_used+1):Nsamples,1])
			sampled_SDLS[,matrix_index]=draw$scale0[1,(Nsamples-sample_used+1):Nsamples,1]^2*(gamma(1+2/draw$nu0[1,(Nsamples-sample_used+1):Nsamples,1])-gamma(1+1/draw$nu0[1,(Nsamples-sample_used+1):Nsamples,1])^2)
			sampled_BVLS[,matrix_index]=1/draw$tnu0[1,(Nsamples-sample_used+1):Nsamples,1]
		}
	}

	for(k in 1:sample_used){
		confmat_tmp=matrix(0,24,24)
		for(worm in 1:nrow(subdata)){
			for(matrix_index in 1:24){
				compvec[matrix_index]=dnorm(subdata$Age[worm],sampled_MLS[k,matrix_index],sqrt(sampled_SDLS[k,matrix_index]+sampled_BVLS[k,matrix_index]),log=T)
			}
			dec_index=which.max(compvec)
			t_orig=which(subdata$T2[worm]==tlist)
			f_orig=which(subdata$Food[worm]==flist)
			index_orig=(t_orig-1)*6+f_orig
			confmat[dec_index,index_orig]=confmat[dec_index,index_orig]+1;
			confmat_tmp[dec_index,index_orig]=confmat[dec_index,index_orig]+1;
		}
		
		dec_pow[k]=sum(diag(confmat_tmp))/sum(confmat_tmp)
	}
	return(list(CM=confmat,DP=dec_pow))

}

decode_ALL <- function(){
	res_F <- vector("list",16) # 4 genes X 4 temperature
	res_T <- vector("list",24) # 4 genes X 6 food
	res_FT <- vector("list",4) # 4 genes

	cat("Decoding Food\n")
	for(g in 1:4){
		for(t in 1:4){
			cat(g,' ',t,"     \r")
			res_F[[(g-1)*4+t]]=Decoding_F(g,t)
		}
	}
	cat("\n")

	cat("Decoding Temperature\n")
	for(g in 1:4){
		for(f in 1:6){
			cat(g,' ',f,"     \r")
			res_T[[(g-1)*6+f]]=Decoding_T(g,f)
		}
	}

	cat("Decoding Food and Temperature\n")
	for(g in 1:4){
		cat(g,"     \r")
		res_FT[[g]]=Decoding_FT(g)
	}

	return(list(decF=res_F,decT=res_T,decFT=res_FT));
}

#all_decoding<-decode_ALL();
load("all_decoding_14052018.RData")

decodingF.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/LS_decodingF.pdf",sep=""),width=12,height=4)
	layout(matrix(1:4,nrow=1))
	par(mar=c(4,6,3,3))

	for(g in 1:4){
		tmp=cbind(all_decoding$decF[[(g-1)*4+1]]$DP,
				  all_decoding$decF[[(g-1)*4+2]]$DP,
 				  all_decoding$decF[[(g-1)*4+3]]$DP,
  				  all_decoding$decF[[(g-1)*4+4]]$DP)

		ylab="";if(g==1) ylab="Food Decoding power";
		boxplot(tmp,outline=F,cex.axis=1.5,ylim=c(0.2,.4),main=glist[g],
				cex.main=1.5,col=gcol[g],cex.axis=1.5,cex.lab=1.5,ylab=ylab,
				names=tlist)
	}
	if(print) dev.off();
}

decodingF.table <- function(print=F,folder=NA){
    tmp=matrix(NA,9,2)
	for(g in 1:3){
		tmp[(g-1)*3+1,1]=mean(all_decoding$decF[[(g-1)*4+1]]$DP);
		tmp[(g-1)*3+1,2]=sd(all_decoding$decF[[(g-1)*4+1]]$DP);

		tmp[(g-1)*3+2,1]=mean(all_decoding$decF[[(g-1)*4+3]]$DP);
		tmp[(g-1)*3+2,2]=sd(all_decoding$decF[[(g-1)*4+3]]$DP);
		
		tmp[(g-1)*3+3,1]=mean(all_decoding$decF[[(g-1)*4+4]]$DP);
		tmp[(g-1)*3+3,2]=sd(all_decoding$decF[[(g-1)*4+4]]$DP);

	}
	if(print) write.table(tmp,file=paste(folder,"/LS_decodingF.txt",sep=""),col.names=F,row.names=F)
		return(tmp)
}
decodingT.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/LS_decodingT.pdf",sep=""),width=12,height=4)
	layout(matrix(1:4,nrow=1))
	par(mar=c(10,6,3,3))

	for(g in 1:4){
		tmp=cbind(all_decoding$decT[[(g-1)*6+1]]$DP,
				  all_decoding$decT[[(g-1)*6+2]]$DP,
 				  all_decoding$decT[[(g-1)*6+3]]$DP,
 				  all_decoding$decT[[(g-1)*6+4]]$DP,
 				  all_decoding$decT[[(g-1)*6+5]]$DP,
  				  all_decoding$decT[[(g-1)*6+6]]$DP)

		ylab="";if(g==1) ylab="Temperature Decoding power";
		boxplot(tmp,outline=F,cex.axis=1.5,ylim=c(0.4,0.8),main=glist[g],
				cex.main=1.5,col=gcol[g],cex.axis=1.5,cex.lab=1.5,ylab=ylab,
				names=flist,las=2)
	}
	if(print) dev.off();
}

decodingFT.plot <- function(print=F,folder=NA){
	if(print) pdf(paste(folder,"/LS_decodingFT.pdf",sep=""),width=12,height=4)
	layout(matrix(1:4,nrow=1))
	par(mar=c(2,2,4,2));
	for(g in 1:4){
		image(x=1:24,y=1:24,z=all_decoding$decFT[[g]]$CM,col=colorRampPalette(colors=c(gcol[g],"white"))(100),main=glist[g])
		grid(4,4,lwd=2)
	}
	if(print) dev.off();
}

robTest <- function(){
	decs<-vector("list",4);
	nrows=length(all_decoding$decF[[1]]$DP)
	for(g in 1:4){
		decs[[g]] = matrix(NA,nrows,4)
		for(t in 1:4){
			decs[[g]][,t] = all_decoding$decF[[(g-1)*4+t]]$DP
		}
	}
	maxdiff <- matrix(NA,nrows,4)
	for(g in 1:4) maxdiff[,g]=apply(decs[[g]],1,function(x) diff(range(x))/mean(x))
	return(maxdiff)
}

robTest.plot <- function(print=FALSE,folder=NA){
	rt = robTest()
	if(print) pdf(paste(folder,"/LS_decodingF_robtest.pdf",sep=""),7,7)
	plot(ecdf(rt[,1]),main="", xlab="Delta/mean",ylab="")
	for(g in 2:4) lines(ecdf(rt[,g]),col=gcol[g])
	if(print) dev.off()
}

robTest2 <- function(){
	decs <- vector("list",3);
	nrows=length(all_decoding$decF[[1]]$DP)
	for(g in 1:3){
		decs[[g]] = matrix(NA,nrows,4)
		for(t in 1:4){
			decs[[g]][,t] = all_decoding$decF[[(g+1-1)*4+t]]$DP-all_decoding$decF[[t]]$DP
		}
	}

	compmat <- vector("list",3)

    for(g in 1:3){
		compmat[[g]]=matrix(NA,4,4)
		for(t1 in 1:4){
			for(t2 in 1:4){
				compmat[[g]][t1,t2] <- sum(decs[[g]][,t1]>decs[[g]][,t2])/nrows
			}
		}
	}
	return(compmat)
}

robTest3 <- function(){
	decs <- vector("list",3);
	nrows=length(all_decoding$decF[[1]]$DP)
	for(g in 1:3){
		decs[[g]] = matrix(NA,nrows,4)
		for(t in 1:4){
			decs[[g]][,t] = all_decoding$decF[[(g+1-1)*4+t]]$DP-all_decoding$decF[[t]]$DP
		}
	}

	compmat <- vector("list",3)

    for(g in 1:3){
		compmat[[g]] <- sum(decs[[g]][,1]>decs[[g]][,2] )
		                   # decs[[g]][,1]>decs[[g]][,3] |
							#decs[[g]][,1]>decs[[g]][,4] )/nrows
						#	decs[[g]][,2]>decs[[g]][,3] |
						#	decs[[g]][,2]>decs[[g]][,4] |
							#decs[[g]][,3]>decs[[g]][,4])

	}
	return(compmat)
}

cv.get <- function(){
	decs<-vector("list",4);
	nrows=length(all_decoding$decF[[1]]$DP)
	for(g in 1:4){
		decs[[g]] = matrix(NA,nrows,4)
		for(t in 1:4){
			decs[[g]][,t] = all_decoding$decF[[(g-1)*4+t]]$DP
		}
	}
	cv <- matrix(NA,nrows,4)
	for(g in 1:4) cv[,g] = apply(decs[[g]],1,function(x) (x)/mean(x))
	return(cv)
}

# Added on 04.01.2019: comparison of lifespan decoding power across genotypes.

LSDPcomp <- function(mode="food", g1=1,g2=2){
	s=length(all_decoding$decT[[1]]$DP)
   
	if(mode=="food"){
		pvec=rep(NA,4)
		for(i in 1:4) {
			pvec[i]=sum(all_decoding$decF[[(g1-1)*4+i]]$DP>all_decoding$decF[[(g2-1)*4+i]]$DP)/s
		}
		names(pvec)=tlist
	} else if(mode=="temp"){
		pvec=rep(NA,6)
		for(i in 1:6) {
			pvec[i]=sum(all_decoding$decT[[(g1-1)*4+i]]$DP>all_decoding$decT[[(g2-1)*4+i]]$DP)/s
		}
		names(pvec)=flist
	}

	return(pvec)
}

LSDPcomp.tables <- function(){
	# t1 WT>daf
	t=list();
	
	u1=LSDPcomp("food",1,2)
	u2=LSDPcomp("food",1,3)
	u3=LSDPcomp("food",1,4)
	u4=LSDPcomp("food",2,4)
	u5=LSDPcomp("food",3,4)

	t$foodcomp=rbind(u1,u2,u3,u4,u5);
	rownames(t$foodcomp)=c(paste(StrainList[1],StrainList[2:4],sep=">"),
	                       paste(StrainList[2:3],StrainList[4],sep=">"))

	u1=LSDPcomp("temp",1,2)
	u2=LSDPcomp("temp",1,3)
	u3=LSDPcomp("temp",1,4)
	u4=LSDPcomp("temp",2,4)
	u5=LSDPcomp("temp",3,4)


	t$tempcomp=rbind(u1,u2,u3,u4,u5);
	rownames(t$tempcomp)=c(paste(StrainList[1],StrainList[2:4],sep=">"),
	                       paste(StrainList[2:3],StrainList[4],sep=">"))

	return(t)
}


LS_EXPR_DP <- function(print=F,folder=NA){
	dec_LS=decodingF.table();
	dec_EX=read.table("~/workspace/mcmc/example/results_expression_04012019/Decoding_Power_fn.txt")
	if(print) pdf(paste(folder,"/LS_EXPR_DPcomp.pdf",sep=""),width=10,height=8)
	par(mar=c(5,5,2,2))
	plot(dec_LS[,1],dec_EX[,1],
	     col=rep(gcol[1:3],1,each=3),
		 pch=rep(0:2,3),cex=2,cex.axis=2,cex.lab=2,
		 xlim=c(0.2,0.5),ylim=c(0.2,0.5),
		 xlab="LS DP",
		 ylab="EXPR DP")
	arrows(x0=dec_LS[,1]-dec_LS[,2],x1=dec_LS[,1]+dec_LS[,2],
	       y0=dec_EX[,1],y1=dec_EX[,1],
		   col=rep(gcol[1:3],1,each=3),code=2,angle=90,length=0)
	arrows(x0=dec_LS[,1],x1=dec_LS[,1],
	       y0=dec_EX[,1]-dec_EX[,2],y1=dec_EX[,1]+dec_EX[,2],
		   col=rep(gcol[1:3],1,each=3),code=2,angle=90,length=0)
	legend("topright",legend=c(15,20,25),pch=0:2,cex=2)
	legend("bottomright",legend=StrainList[1:3],fill=gcol[1:3],cex=2)
    if(print) dev.off();
}
