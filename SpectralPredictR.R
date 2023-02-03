#Kyle R Kovach
#SpectroPredictR
#December 1, 2020

#----Load Libraries----
# install.packages("doParallel")
# install.packages("spectrolab")
# install.packages("dplyr")
# install.packages("")

require(doParallel)
require(spectrolab)
require(dplyr)
require(tidyr)

#----Initialize----
set.seed(112345)
cores=detectCores()-1

specfolder=choose.dir(default="",caption="Please select main folder containing subdirectories of spectral files.")

split_path=function(path) {
  if (dirname(path) %in% c("_", path)) return(basename(path))
  return(c(basename(path),split_path(dirname(path))))
}

#----Import Spectra----
for (j in c("asd","sed","sig"))
{
  if (exists('finaloutput')){rm(finaloutput)}
  folderlist=unique(sub("/[^/]+$","",list.files(path=specfolder,
                                                pattern=paste0(".",j,"$"),
                                                recursive=TRUE,
                                                full.names = TRUE)))
  if (length(folderlist)!=0){
    cl=makeCluster(cores)
    registerDoParallel(cl)
    
    finaloutput=foreach(i=folderlist,
                        .combine=rbind,
                        .packages=c('spectrolab'))%dopar%
      {
        spectra=read_spectra(path=i,format=j)
        data.frame(spectra,folder=i,format=j,check.names=FALSE)->sedadjust
        if (j=="sed"){
          sedadjustT=sedadjust[,3:ncol(sedadjust)]
          sedadjustH=sedadjust[,1:2]
          seqcol=seq(from=1,to=2150,by=2)
          sedadjustTe=sedadjustT[,seqcol]/100
          sedadjustTo=sedadjustT[,-seqcol]
          sedadjustF=cbind(sedadjustH,sedadjustTe,sedadjustTo)
          TestOutput=sedadjustF[,names(sedadjust)]
        }
      };FST=as_spectra(finaloutput[,1:2152],name_idx=1)
    
    stopCluster(cl)
    
    if (j=="asd"||j=="sig"){
      splice_val=ifelse(j=="asd",c(1000, 1800),c(990, 1900))
      FST=match_sensors(
        FST,
        splice_at=splice_val
      )
    };FSTdf=data.frame(FST,check.names=FALSE)
    if (!exists('FFO')){FFO=FSTdf[0,]}
    if (!exists('RT')){RT=finaloutput[0,c(1,2153:2154)]}
    FFO=rbind(FFO,FSTdf);RT=rbind(RT,finaloutput[,c(1,2153:2154)])
  }};FSTfinal=as_spectra(FFO,name_idx=1);FSTfinal_norm=normalize(FSTfinal)

FST_f=data.frame(FSTfinal);FST_f=cbind(RT,FSTfinal)
FST_n=data.frame(FSTfinal_norm);FST_n=cbind(RT,FSTfinal_norm)

#----Process Spectra and Vector Normalize----
headcount=ncol(FST_f)-2151
specstart=headcount+1

cl=makeCluster(cores);registerDoParallel(cl);sampledata_VN=foreach(i=1:nrow(FST_f),
                                                                   .combine=rbind)%dopar%
  {
    a=FST_f[i,]
    b=sqrt(rowSums(a[,specstart:ncol(a)]^2))
    d=a[,-c(1:headcount)]
    e=a[,c(1:headcount)]
    f=d/b
    cbind(e,f)
  }
stopCluster(cl)

sampledata_wav=sampledata_VN[,-c(1:headcount)]
sampledata_head=sampledata_VN[,c(1:headcount)]
sampledata_5nm_wav=sampledata_wav[,seq_len(ncol(sampledata_wav)) %% 5 == 1]
sampledata_5nm=cbind(sampledata_head,sampledata_5nm_wav)
sampledata_VN=sampledata_5nm

#----Predict Traits----

modeldirectory=list.files(path=choose.dir(),
                          pattern=".csv$",
                          recursive=FALSE,
                          full.names = TRUE,
                          caption="Please select main folder containing either fresh or dry models (based on spectra being processed).")

cl=makeCluster(cores);registerDoParallel(cl);finished_output=foreach(model=modeldirectory,
                                                                     .combine=rbind)%dopar%
  {
    modelname=split_path(model)[1]
    modelcoeffs_mat=data.matrix(read.csv(model, header = TRUE, row.names = 1,sep = ","))
    sampledata_VN_tail=sampledata_VN[,-c(1:headcount)]
    sampledata_VN_head=sampledata_VN[,c(1:headcount)]
    sampledata_VN_tail$constant=1
    sampledata_VN_tail=cbind(intercept=sampledata_VN_tail[,ncol(sampledata_VN_tail)],sampledata_VN_tail[,c(1:ncol(sampledata_VN_tail)-1)])
    sampledata_VN_tail_mat=as.matrix(sampledata_VN_tail)
    modeloutput=sampledata_VN_tail_mat %*% t(modelcoeffs_mat)
    modeloutput_bind=cbind(sampledata_VN_head,as.data.frame(modeloutput))
    data.frame(modelname=modelname,
               sampledata_VN_head,
               t_mean=rowMeans(modeloutput),
               t_std=apply(modeloutput,1,sd),
               as.data.frame(modeloutput),
               check.names=FALSE)
  };stopCluster(cl)

predsub=finished_output[,c(1:2,6:7)]

predspread=pivot_wider(data=predsub,id_cols=sample_name,names_from=modelname,values_from=c("t_mean","t_std"))
predunlist=as.data.frame(unnest(predspread))
colnames(predunlist)=sub(".csv", "", colnames(predunlist))
write.csv(predunlist,"SpectroPredictR_finished_model_output.csv",row.names=FALSE)