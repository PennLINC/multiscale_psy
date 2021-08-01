# matlab spits out FC features to 14 decimal places... there is enough to load without that boggis
library(vroom)
fc<-vroom('/cbica/projects/pinesParcels/results_psy_master_fcfeats_Sub.csv')

# set colnames from first row
colnames(fc)<-unlist(fc[1,])

# remove row of colnames
fc<-fc[-c(1),]

# set to match demographics
colnames(fc)[6889]<-'scanid'

# re-tuck 707 into spot
scanid<-seq(1,790)
scanid[1:706]<-fc$scanid[1:706]
# null value
scanid[707]<-99999999
scanid[707:789]<-fc$scanid[708:790]
fc$scanid<-scanid

# round ridiculous number of decimal points
fc[] <- lapply(fc, function(x) {
  if(is.character(x)) round(as.numeric(as.character(x)),digits=3) else x
})

# re-write the rounded version
write.csv(fc,'/cbica/projects/pinesParcels/results_psy/aggregated_data/fc/master_fcfeats_Sub_rounded.csv')
