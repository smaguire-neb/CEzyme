#' Find the peaks in the standards channel
#' @description This function locates candidate standard peaks using the findpeaks
#'    function in the pracma package. The inital peak thershold is set at
#'    1/6 of the maximum value in the standards channel. If 9 peaks are found at this
#'    thershold then the standards are returned. If less than 9 peaks are found then
#'    the thershold is iterativly lowered until at least 9 peaks are found. If more than
#'    9 peaks are found then models are fit for all possible sets of 9 peaks. The set of
#'    peaks with the highest R squared value is chosen as the true peaks. This function
#'    is called by get.std.peak.
#' @param ch5 Column that contains the data from channel 5 (Std channel)
#' @return Dplyr data tables containing the location and size of the
#'    standard peaks for each file.
#' @seealso \code{\link{get.std.peak,find.best.peak,read.data}}
#' @export
get.peaks<-function(ch5){
  std.peaks<-findpeaks(ch5,minpeakheight = max(ch5)/6,nups=2)
  if(nrow(std.peaks) == 9) {
    return(std.peaks)
  } else if(nrow(std.peaks<9)) {
    while(nrow(std.peaks) < 9 ){
      min<- (max(ch5)/4) - 250
      std.peaks<-as_tibble(findpeaks(ch5,minpeakheight = min,nups=2))
    }
  }
  if(nrow(std.peaks) > 9){
    std.peaks <- find.best.peaks(std.peaks) %>%
      rename(V1=height,V2=position,V3=start,V4=end)
  }
  return(std.peaks)
}

#' Find the peaks in the standards channel
#' @description This function is used when more than 9 standard peaks are found.
#'    It fits all possible sets of linear models with 9 datapoints and choses the one that
#'    has the lowest deviance.
#' @param cand.peaks table containing the candidate peaks where the number of peaks is more than 9.
#' @return Dplyr table containing the best set of 9 peaks.
#' @seealso \code{\link{get.peaks,find.best.peak,read.data}}
#' @export
find.best.peaks<-function(cand.peaks){
  cand.peaks <-
    cand.peaks %>%
    as_tibble() %>%
    rename(height=V1,position=V2,start=V3,end=V4)
  combos<-combn(1:nrow(cand.peaks),9,simplify = F)
  r.squared<-
    map(combos, function(x) filter(cand.peaks, row_number() %in% x)) %>%
    map(function(x) mutate(x, std = c(15,20,25,35,50,62,80,110,120))) %>%
    map(function(x) lm(std~position,data=x)) %>%
    map(summary) %>%
    map_dbl("r.squared")
  best<-which(r.squared == max(r.squared))
  best.combo<-combos[[best]]
  return(cand.peaks[best.combo,])
}

#' Find the peaks in the standards channel
#' @description This function identifies the standard peaks in the
#'    Liz channel (Channel 5). It uses the get peaks function to find the
#'    correct peaks.
#' @param data.list List of data frames generated by the read.data function
#' @return List of dplyr data tables containing the location and size of the
#'    standard peaks for each file.
#' @seealso \code{\link{get.peaks,find.best.peak,read.data}}
#' @export

get.std.peak<-function(data.list){
  std.peaks<-
    map(data.list, function(x) select(x,ch5)) %>%
    map(unlist) %>%
    map(get.peaks) %>%
    map(as_tibble) %>%
    map(function(x) rename(x,height=V1,position=V2,start=V3,end=V4)) %>%
    map(function(x) mutate(x,std = c(15,20,25,35,50,62,80,110,120)))
  std.peaks
}
