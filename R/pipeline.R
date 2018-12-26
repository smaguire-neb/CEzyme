
process.data<- function(data.dir,channel){
  data<-read.data(data.dir)
  stds<-get.std.peak(data)
  models<-get.models(stds)
  data.sized<-apply.size(data,models,stds)
  peaks<-discover.peaks(data.sized,channel = "ch1")
  areas<-get.peak.area(peaks,data.sized,channel = "ch1")
  return(list(data.sized,peaks,areas))
}
