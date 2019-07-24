#Discover Life Site Climate (Temperature) Data Acquisition
#For Pheno-mismatch moth & caterpillar phenology
#Elise Larsen, Ries Lab, 2019

#Script to Access Daymet Temperature Data for Discover Life Mothing sites listed in moth-sites.csv file
#(Moth data from https://www.discoverlife.org/)
#climate data from Daymet https://daymet.ornl.gov/

#libraries needed: None
library(daymetr)
library(dplyr)

#data inputs
#site data frame created from https://www.discoverlife.org/moth/report.html
dlsites<-read.csv("data/moth-sites.csv", header=T)
#function to calculate degree days (single-sine approximation)
source("scripts/degday1.R")



#create data frame for data in Discover Life Moth Research - specific to "table 5"
#html example of data table: https://www.discoverlife.org/moth/data/table5_34.0_-83.4.html
dl.climate<-data.frame(year=numeric(0),julian.day=numeric(0),tmaxC=numeric(0),tminC=numeric(0),site=character(0))


#Loop through sites and collect temperature data 
#Temporal rang begins jan 1 of the first year of data collection at that site to the year prior to today's year.
#Daymet years always have 365 days; daymet doesn't account for leap years

for(i in 1:nrow(dlsites)) {
  daymetX<-download_daymet(site = dlsites$site[i], lat = dlsites$lat[i], lon = dlsites$long[i],start = dlsites$startyear[i], end = as.numeric(format(Sys.Date(), "%Y"))-1,path = tempdir(), internal = TRUE, silent = FALSE, force = FALSE)
  daymetX$data$site<-as.character(dlsites$site[i])
  dl.climate<-rbind(dl.climate,daymetX$data[,c(1,2,7,8,10)])
  remove(daymetX)
}
names(dl.climate)<-c("year","julian.day","tmaxC","tminC","site")
#For all days, calculate growing degree day using accumulation thresholds 10C (minimum) and 30C (maximum)
#No growing degree days are accumulated outside this 10C - 30C range
dl.climate$dayGDD<-0
for(i in 1:nrow(dl.climate)) {
  dl.climate$dayGDD[i]<-round(degreedays(dl.climate$tmin[i],dl.climate$tmax[i],10,30),1)
}
#accumulate growing degree days across days within each year at each site
dl.climate<-dl.climate %>%
  group_by(site,year) %>%
  mutate(cumGDD=cumsum(dayGDD))

#save(dl.climate,file="data/moth_temp+gdd.RData")
#write.csv(dl.climate,file="data/moth_temp+gdd.csv")

