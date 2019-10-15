require(Hmisc)
NPAR=9
NMOD=512

tlist<-c("15","20","25");
flist<-c("sb","27","67","68","29","110");
mat<-matrix(0,nrow=18,ncol=512);
parmat<-matrix(NA,ncol=11,nrow=18);
parmat2<-matrix(NA,ncol=11,nrow=18); # most frequent model
for(t in 1:3){
    for(f in 1:6){
        tab<-read.table(paste("results/T",tlist[t],"_",flist[f],".dat",sep=""))
        names(tab)<-c("t1","t2","t3","t4","t5","t6","t7","t8","t9","model","chi2");
        minchi2=10;
        for(k in 1:nrow(tab)){
            if(mat[f+(t-1)*6,tab[k,10]+1]<(9-tab[k,11])/9.){
                mat[f+(t-1)*6,tab[k,10]+1]=(9-tab[k,11])/9.;
            }
        }
        parmat[f+(t-1)*6,]=as.numeric(tab[which.min(tab$chi2),]);
        mode=as.numeric(names(sort(table(tab$model),decreasing=T)[1]))
        parmat2[f+(t-1)*6,]=as.numeric(tab[mode+1,])
    }
}
2
best3<-(apply(mat,1,function(x) order(x,decreasing=T)[1:3]))
parmatrast=matrix(as.numeric(cut(parmat[,1:9],breaks=c(0,.01,1,2,600))),ncol=9,nrow=18)
parmatrast2=matrix(as.numeric(cut(parmat2[,1:9],breaks=c(0,.01,1,2,1600))),ncol=9,nrow=18)

image(parmatrast,col=rev(c(rgb(0,0,1),rgb(.1,.5,1),"orange","red")),xaxt='n',yaxt='n')

prob <- function(W,gt=1){
  wxx=W[1];
  wxy=W[2];
  wxz=W[3];
  wyx=W[4];
  wyy=W[5];
  wyz=W[6];
  wzx=W[7];
  wzy=W[8];
  wzz=W[9];

  if(gt==2) {
      wxx=1; wxy=1; wxz=1;
      wzx=1; wzy=1; wzz=1;
  }
  if(gt==3){
      wyx=1; wyy=1; wyz=1;
  }
  if(gt==4){
      wxx=1; wxy=1; wxz=1;
      wyx=1; wyy=1; wyz=1;
      wzx=1; wzy=1; wzz=1;
  }

  res=matrix(
         c(-3  , 1      , 1       , 0                  , 1       , 0                  , 0                  , 0,
           wxx , -2-wxx , 0       , 1                  , 0       , 1                  , 0                  , 0,
           wyy , 0      , -wyy-2  , 1                  , 0       , 0                  , 1                  , 0,
           0   , wxy*wyy, wxx*wyx , -1-wxy*wyy-wxx*wyx , 0       , 0                  , 0                  , 1,
           wzz , 0      , 0       , 0                  , -2-wzz  , 1                  , 1                  , 0,
           0   , wzz*wxz, 0       , 0                  , wzx*wxx , -1-wzx*wxx-wzz*wxz , 0                  , 1, 
           0   , 0      , wzz*wyz , 0                  , wyy*wzy , 0                  , -1-wyy*wzy-wzz*wyz , 1,
           0   , 0      , 0       , wxz*wyz*wzz        , 0       , wxy*wyy*wzy        , wzx*wyx*wxx        , -wxy*wyy*wzy-wzx*wyx*wxx-wxz*wyz*wzz),
         ncol=8,nrow=8);
  res[8,]=rep(1,8);

  pall <- solve(res,c(0,0,0,0,0,0,0,1))  
  
  return(list(c(pall[2]+pall[4]+pall[6]+pall[8], pall[3]+pall[4]+pall[7]+pall[8], pall[5]+pall[6]+pall[7]+pall[8]), pall,res))
}


