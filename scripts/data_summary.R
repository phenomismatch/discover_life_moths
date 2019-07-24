#Discover Life Moth Data Summary
#For Pheno-mismatch moth & caterpillar phenology
#Elise Larsen, Ries Lab, 2019
##Currently a collection of data summary tables and visualizations
##Currently, site is a factor & order is alphabetic - need to switch order to match latitude


#Script to Access Discover Life Mothing Data using sites listed in moth-sites.csv file
#Data from https://www.discoverlife.org/

#libraries needed: 
library(dplyr)
library(ggplot2)
library(usmap)

#data inputs:
dl.sites<-read.csv("data/moth-sites.csv", header=T)
load("data/moth-abundance.RData")
# or read.csv("moth-abundance.csv", header=T)
load("data/moth-temp+gdd.RData")
# or read.csv("moth-temp+gdd.csv", header=T)


#NOTE: The moth abundance data table conflates days with no moths with days when data were not collected (0 photos)
#Filter to days with moth photos 
moth.abundance<-subset(moth.abundance, photos>0)


#Map sites
usa <- map_data("usa") 
ggplot() + geom_polygon(data = usa, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3) + geom_point(data = dl.sites, aes(x = long, y = lat), color = "yellow", size = 3)

#Tally # days of observations by site, month, & year
temporalxsite<-moth.abundance %>%
  group_by(site,year,month) %>%
  tally()

ggplot(aes(year), data=temporalxsite) + geom_bar(aes(fill = site))
ggplot(aes(month), data=temporalxsite) + geom_bar(aes(fill = site))
#to do: add labels, switch site level order to match latitude

#Tally # days of observations by site, month, & year filtered to April-July
temporaltarget<-moth.abundance %>%
  group_by(site,year,month) %>%
  filter(month %in% c(4:8)) %>%  #include april through august
  tally()

ggplot(aes(year), data=temporaltarget) + geom_bar(aes(fill = site))
ggplot(aes(month), data=temporaltarget) + geom_bar(aes(fill = site))
#to do: add labels, switch site level order to match latitude


#ID year range for each site
siteyears<-moth.abundance %>%
  group_by(site) %>%
  summarize(minyear=min(year), maxyear=max(year))

#Accumulate # moths across days for each site*year
moth.data<-moth.abundance %>%
  group_by(year, site) %>%
  mutate(cumAbund=cumsum(photos))

#Calculate proportion of moths relative to total for site*year
moth.data<-moth.data %>%
  group_by(year, site) %>%
  mutate(propMoth= cumAbund/max(cumAbund))

head(moth.data[,c(1:5,7,15:18)])


#add temperature anomaly (preliminary)
dl.climate<-dl.climate %>%
  group_by(site,julian.day) %>%
  mutate(typGDD=mean(cumGDD),anomaly=(cumGDD-typGDD)/typGDD)

#Merge moth & climate data
moths<-merge(moth.data, dl.climate, c("site","year","julian.day"))

#VIS: Plot of when data exist for each site
##### IN SPRING MEETING PRESENTATION
ggplot(data=climate.data, aes(julian.day,site)) +
  geom_point(aes(julian.day,site, color=cumGDD))

#Moth accumulation plots
ggplot(data=moth.data, aes(julian.day,propMoth, color=year)) +
  geom_point(aes(julian.day,propMoth, color=year))

ggplot(data=moths, aes(julian.day,propMoth, color=cumGDD)) +
  geom_point(aes(julian.day,propMoth, color=cumGDD))

ggplot(data=moths, aes(julian.day,propMoth, color=anomaly)) +
  geom_point(aes(julian.day,propMoth, color=anomaly))


## Preliminary comparison of 50th percentile of moth abundance & climate proxy
## In current formuation, Spatial variation much greater than interannual variation
# Calculate when half of moth photos are accumulated
midcurve<-data.frame(site=numeric(0),year=numeric(0),DOY10=numeric(0),DOY50=numeric(0))
for(j in unique(moth.data$site)) {
  temp2<-moth.data[moth.data$site==j,]
  for(i in unique(temp2$year)) {
    temp1<-temp2[temp2$year==i,]
    doy1<-temp1$julian.day[which(temp1$propMoth>0.499)[1]-1]
    doy2<-temp1$julian.day[which(temp1$propMoth>0.499)[1]]
    pr1<-temp1$propMoth[which(temp1$propMoth>0.499)[1]-1]
    pr2<-temp1$propMoth[which(temp1$propMoth>0.499)[1]]
    if(which(temp1$propMoth>0.099)[1]==1) {
      doy10a<-90
      doy10b<-temp1$julian.day[which(temp1$propMoth>0.099)[1]]
      pr10a<-0
      pr10b<-temp1$propMoth[which(temp1$propMoth>0.099)[1]]
    } else {
      doy10a<-temp1$julian.day[which(temp1$propMoth>0.099)[1]-1]
      doy10b<-temp1$julian.day[which(temp1$propMoth>0.099)[1]]
      pr10a<-temp1$propMoth[which(temp1$propMoth>0.099)[1]-1]
      pr10b<-temp1$propMoth[which(temp1$propMoth>0.099)[1]]
    }
    midcurve[nrow(midcurve)+1,]<-c(j,i,round((.1-pr10a)/((pr10b-pr10a)/(doy10b-doy10a)))+doy10a,round((.5-pr1)/((pr2-pr1)/(doy2-doy1)))+doy1 )
  }
}
midcurve$sitename<-factor(midcurve$site)
midcurve$year<-as.numeric(midcurve$year)
midcurve$DOY10<-as.numeric(midcurve$DOY10)
midcurve$DOY50<-as.numeric(midcurve$DOY50)

## Start with simple temperature proxy - accumulated GDD by June 1
gdd152<-dl.climate[dl.climate$julian.day==152,]

test.data1<-merge(midcurve, gdd152)

#Simple exploratory linear model
lm.moth1<-lm(DOY50~cumGDD, data=test.data1)
regression.moth = summary(lm.moth1) #save regression summary as variable
names(regression.moth) #get names so we can index this data
a= regression.moth$coefficients["(Intercept)","Estimate"] #grab values
b= regression.moth$coefficients["cumGDD","Estimate"]
#abline(a,b) #add the regression line
summary(lm.moth1)

###VIS FROM SPRING PHENO MEETING
ggplot(test.data1, aes(x=cumGDD, y=DOY50, color=sitename)) + 
  geom_point(aes(cumGDD,DOY50,color=sitename)) + 
  geom_smooth(method=lm, color='#2C3E50', level=0.95) +
  ggtitle("50th Percentile Moth Abundance") +
  labs(x="GDD accumulated by June 1", y="Julian Day of Half Moth Accumulation (April-July)", color="Site")



