library(rjags)

args <- commandArgs(trailingOnly=TRUE)
chainused=1000
data<-read.table("~/workspace/data/dataset-add-9/20170524_FEDt.dat");
names(data)<-c("Batch","date","GT","day","food","T","ADF","ASI","NSM","time");

Diff <- function(x, start) as.numeric(x - as.Date(cut(start, "year")))
flist <- c(1, 2e+07, 6.3e+07, 6.3e+08,2e+09, 1.1e+10);
fcontrol=2e+09;
glist <- c("QL196","QL404","QL402","QL435");
collist <- c("black","blue","red","purple");
tcol<-c("#FFFF00","#FFC800","#FF9600");
neuron_lab <- c("ADF","ASI","NSM")
tlist <-c(15,20,25);
lut<-read.table("~/workspace/Analysis/ExpressionAnalysis/R-scripts/combined_batch_list_v2.0.dat")
#lut$V1 <- sapply(as.character(lut$V1),function(x) strsplit(x,"-")[[1]][2])
fcol <- rgb(0,seq(0,1,length=6),seq(1,0,length=6))
beads_red <- read.table("/home/diana/workspace/NeuroShed/NeuroShedBeads_v3/beads_red_13022018.dat")
beads_red$DATE=as.Date(beads_red$DATE,"%Y-%m-%d")
start_date=beads_red$DATE[1]
beads_red$DATE=Diff(beads_red$DATE,beads_red$DATE[1])
beads_red$NVAL=beads_red$NVAL/1e3
#names(beads_red)<-c("index","NVAL","BG","NBEADS","FOLDER","DATE");
beads_green <- read.table("/home/diana/workspace/NeuroShed/NeuroShedBeads_v3/beads_green_13022018.dat")
beads_green$DATE=as.Date(beads_green$DATE,"%Y-%m-%d")
beads_green$DATE=Diff(beads_green$DATE,beads_green$DATE[1])
beads_green$NVAL=beads_green$NVAL/1e3
#names(beads_green)<-c("index","NVAL","BG","NBEADS","FOLDER","DATE");

getdata <- function(Genotype, Food, Temperature){
		tmpset<-subset(data,GT==Genotype & day==6 & (Batch %in% lut$V1) & (T == Temperature) & (food == Food), select=c("Batch","T","food","ADF","ASI","NSM"));
		g=match(Genotype,glist);		

		tabtemp<-table(tmpset$Batch);
        batches<-names(tabtemp)[tabtemp>0]

		for(i in 1:length(batches)){
			tmpset$bid[tmpset$Batch==batches[i]]=i	
		}

		tmpset[,c("ADF","ASI","NSM")]=tmpset[,c("ADF","ASI","NSM")]/1e6
		tmpset$DATE=Diff(as.Date(strtrim(tmpset$Batch,8),"%Y%m%d"),start_date)
		meanmat=apply(tmpset[,c("ADF","ASI","NSM")],2,mean)
		meanmat[1]=meanmat[1]/mean(beads_red$NVAL)
		meanmat[2]=meanmat[2]/mean(beads_green$NVAL)
		meanmat[3]=meanmat[3]/mean(beads_red$NVAL)
		cPostRed=getBeadsPosterior_red(1,1,1)
		cPostGreen=getBeadsPosterior_green(1,1,1)

		return(list(y=tmpset[,c("ADF","ASI","NSM")],
					ybid=tmpset$bid,
					ydate=tmpset$DATE,
					#meanmat=meanmat,
					cPostRed_a=cPostRed$a,
					cPostRed_b=cPostRed$b,
					cPostRed_mu=cPostRed$mu,
					cPostRed_lambda=cPostRed$lambda,
					cPostGreen_a=cPostGreen$a,
					cPostGreen_b=cPostGreen$b,
					cPostGreen_mu=cPostGreen$mu,
					cPostGreen_lambda=cPostGreen$lambda,
					#yC=data.frame(beads_red$NVAL,beads_green$NVAL,beads_red$DATE),
					NW=nrow(tmpset),
					NWC=nrow(beads_green),
					NB=length(batches)))
} 

getBeadsPosterior_red <- function(a0=1,b0=1,sigma0=1){
	n=nrow(beads_red)
	y=beads_red$NVAL
	X=cbind(1,beads_red$DATE)
	Lambda0 = diag(sigma0,2)
	Lambda  = t(X)%*%X+Lambda0
	beta_hat = solve(t(X)%*%X)%*%t(X)%*%y
	mu0 = c(1,0)
	mu  = solve(t(X)%*%X+Lambda0)%*%(t(X)%*%X%*%beta_hat+Lambda0%*%mu0)
	a=a0+n/2
	b=b0+0.5*(t(y)%*%y+t(mu0)%*%Lambda0%*%mu0-t(mu)%*%Lambda%*%mu)

	return(list(a=a,b=b,mu=mu,lambda=Lambda))
}

getBeadsPosterior_green <- function(a0=1,b0=1,sigma0=1){
	n=nrow(beads_green)
	y=beads_green$NVAL
	X=cbind(1,beads_green$DATE)
	Lambda0 = diag(sigma0,2)
	Lambda  = t(X)%*%X+Lambda0
	beta_hat = solve(t(X)%*%X)%*%t(X)%*%y
	mu0 = c(1,0)
	mu  = solve(t(X)%*%X+Lambda0)%*%(t(X)%*%X%*%beta_hat+Lambda0%*%mu0)
	a=a0+n/2
	b=b0+0.5*(t(y)%*%y+t(mu0)%*%Lambda0%*%mu0-t(mu)%*%Lambda%*%mu)

	return(list(a=a,b=b,mu=mu,lambda=Lambda))
}

model_beads = "model {
			for(i in 1:NWC){
				yC[i,1:2] ~ dmnorm( C[1,1:2] + C[2,1:2]*yC[i,3], Ty)
			}
		
			C[1,1] = mu[1]
			C[1,2] = mu[2]
			C[2,1] = mu[3]
			C[2,2] = mu[4]

			mu ~ dmnorm(c(1,1,0,0),TC)
			
			Ty ~ dwish(R2,4)
			TC ~ dwish(R4,6)
			
		    
			for(i in 1:2){
				for(j in 1:2){
				   	R2[i,j]=equals(i,j)*1e-3
				}
			}
			
			for(i in 1:4){
				for(j in 1:4){
				   	R4[i,j]=equals(i,j)*1e-3
				}
			}

	}"

model7="model {


	     # NORMALIZATION MODEL

		    rCt ~ dgamma(cPostRed_a,cPostRed_b)
		    gCt ~ dgamma(cPostGreen_a,cPostGreen_b)
			C[1:2,1] ~ dmnorm(cPostRed_mu,rCt*cPostRed_lambda)
			C[1:2,2] ~ dmnorm(cPostGreen_mu,rCt*cPostGreen_lambda)

         # EXPRESSION MODEL
		 ## For each neuron we use a different covariance (precision) matrix for the 18 conditions.
		 ## This reflect the underlying assumption that the mean response profiles of the neurons are not
         ## correlated.

            rat[1] ~ dmnorm(globmean[1],TglobADF)
            rat[2] ~ dmnorm(globmean[2],TglobASI)
            rat[3] ~ dmnorm(globmean[3],TglobNSM)

		 ## The mean normalized expression values for each trial are then sampled from a normal distribution
         ## with precision matrix which reflect the experimental variability from trial to trial. Here we assume that this experimental variability, characterized by T.rat depends on the environmental condition.

				T.rat[1:3,1:3] ~ dwish(R3.rat,5)
			    TE[1:3,1:3] ~ dwish(R3,5)

            for(i in 1:NB){
				rr[i,1:3] ~ dmnorm(rat[1:3],T.rat)
				mu.E[i,1] = rr[i,1]*(C[1,1] + C[2,1]*ydate[i])
				mu.E[i,2] = rr[i,2]*(C[1,2] + C[2,2]*ydate[i]) 
				mu.E[i,3] = rr[i,3]*(C[1,1] + C[2,1]*ydate[i])
			}

         ## Once assigned the average normalized expression, the observations are sampled from a normal distribution with covariance matrix which includes the correlation among the three neurons which also depends on the environmental condition.

            for(i in 1:NW){
				y[i,1:3] ~ dmnorm(mu.E[ybid[i],1:3],TE)
			}

				for(j in 1:3){
					globmean[j] ~ dunif(5,100)
				}
				TglobADF ~ dgamma(2,40)
				TglobASI ~ dgamma(2,40)
				TglobNSM ~ dgamma(2,40)

		    for(i in 1:3){
				for(j in 1:3){
				   	R3[i,j]=equals(i,j)*1e-0
					R3.rat[i,j]=equals(i,j)*1e-0
				}
			}
		    
		}"

gensamples_single <- function(g,f,t,iter=20000,ch=1,thin=1){
				datatest<-getdata(glist[g],flist[f],tlist[t]);
				varnames=c("rat","TglobADF","TglobASI","TglobNSM","T.rat","TE","mu.E","C");
				burn_in=1;
				steps=iter
				fileConn=file("model7.tmp")
				writeLines(model7,fileConn);
				close(fileConn)
				
				m=jags.model(file="model7.tmp",data=datatest,n.chains=ch,
#				             inits=list(.RNG.name="base::Wichmann-Hill",.RNG.seed=as.numeric(seed)));
				             inits=list(.RNG.name="base::Wichmann-Hill",.RNG.seed=as.numeric(args[1])));
				update(m,burn_in)
				draw=jags.samples(m,steps,thin=thin,variable.names=varnames);
	return(draw)
}

gensample_beads <- function(iter=10000,thin=1){
	varnames=c("C");
	burn_in=1000;
	steps=iter
	fileConn=file("model_beads.tmp")
	writeLines(model_beads,fileConn);
	close(fileConn)

	m=jags.model(file="model_beads.tmp",data=list(yC=data.frame(beads_red$NVAL,beads_green$NVAL,beads_red$DATE),NWC=nrow(beads_red)))
	update(m,burn_in)
	draw=jags.samples(m,steps,thin=thin,variable.names=varnames);
	return(draw)
}


gensamples <- function(iter=40000,ch=2,thin=1){
	res_list=vector("list",length=4)
	
	for(g in 1:4){
		res_list[[g]]$rat=array(dim=c(18,3,iter/thin,1))
		res_list[[g]]$TglobADF=array(dim=c(18,iter/thin,1))
		res_list[[g]]$TglobASI=array(dim=c(18,iter/thin,1))
		res_list[[g]]$TglobNSM=array(dim=c(18,iter/thin,1))
		res_list[[g]]$T.rat=array(dim=c(3,3,18,iter/thin,1))
		res_list[[g]]$TE=array(dim=c(3,3,18,iter/thin,1))
		res_list[[g]]$C=array(dim=c(2,2,18,iter/thin,1))

		for(t in 1:3){
			for(f in 1:6){
				ind=(t-1)*6+f
				cat(glist[g],' ',tlist[t],' ',flist[f],"\n")
				res_tmp<-gensamples_single(g,f,t,iter,ch,thin)
				res_list[[g]]$rat[ind,,,1]=res_tmp$rat
				res_list[[g]]$TglobADF[ind,,1]=res_tmp$TglobADF
				res_list[[g]]$TglobASI[ind,,1]=res_tmp$TglobASI
				res_list[[g]]$TglobNSM[ind,,1]=res_tmp$TglobNSM
				res_list[[g]]$T.rat[,,ind,,1]=res_tmp$T.rat
				res_list[[g]]$TE[,,ind,,1]=res_tmp$TE
				res_list[[g]]$C[,,ind,,1]=res_tmp$C
			}
		}			
				
	}
	return(res_list)
}

all_samples=gensamples(50000,1,10)
save(all_samples,file=paste("all_samples_28082018_ch",args[1],".RData",sep=""))

