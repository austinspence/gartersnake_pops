rm(list=ls())
setwd("D:\\Project_CacheValley_Thamnophis")
library(reshape2)
library(ggplot2)
library(dplyr)
library(car)
library(e1071)
library(corrplot)
library(lubridate)
library(gtable)
library(zoo)
library(stringr)

d <- read.csv("cache-thamnophis-data.csv")

names(d) <- c("id","date","year","season","subsite","site","juv.sex",
              "sex","sex.juv.only","cort_bl","cort_ps","mass","svl","tl",
              "bka_bl","bka_ps","test","species","glu_bl","glu_ps","notes")

d$id <- as.factor(as.character(d$id))
d$date <- as.Date(as.character(d$date),format="%m/%d/%Y")

# deal with sexes & ages
d$sexage <- paste(d$juv.sex,d$sex,d$sex.juv.only)
d$age[d$sex == "NEO"|d$sex == "JUV"] <- 'juv'
d$age <- 'adult'; d$age[grepl("JUV",d$sexage) == T] <- 'juv'
d <- dplyr::select(d, -c(juv.sex,sex,sex.juv.only))
d$sex <- NA
d$sex[grepl("M",d$sexage) == T] <- 'male'
d$sex[grepl("F",d$sexage) == T] <- 'female'
d <- dplyr::select(d, -c(sexage))

d$svl[d$svl == 44] <- 440
d$mass[d$mass == '77(65)'] <- 65 ## 77 is mass with food (10 minnows), 65 is mass without
d$mass <- as.numeric(d$mass)
d$glu_bl[d$glu_bl == 'LO'] <- NA ## TODO: check this with Lori
d$glu_bl <- as.numeric(d$glu_bl)
d$glu_ps[d$glu_ps == 'LO'] <- NA ## TODO: check this with Lori
d$glu_ps <- as.numeric(d$glu_ps)
d$bka_bl[d$bka_bl > 100] <- 100 ## TODO: check this with Lori
d$bka_ps[d$bka_ps > 100] <- 100 ## TODO: check this with Lori
d$test[d$test > 300] <- NA ## TODO: check if taking out these two outliers is OK with all
d$cort_bl[d$cort_bl > 300] <- NA ## TODO: check if taking out this outliers is OK with all

d$cort_react <- d$cort_ps-d$cort_bl
d$bka_react <- d$bka_ps-d$bka_bl
Z<-lm(svl~mass,data=d)
d[names(Z$residuals),"bodycon"]<-Z$residuals

hist(d$svl[d$age == 'juv']); min(d$svl[d$age == 'juv'],na.rm=T); max(d$svl[d$age == 'juv'],na.rm=T)
hist(d$svl[d$age == 'adult']); min(d$svl[d$age == 'adult'],na.rm=T); max(d$svl[d$age == 'adult'],na.rm=T)
# TODO: check that most JUV are smaller than most adults

ggplot(d,aes(x=svl,fill=age))+geom_histogram()

## meteorological data

precip <- read.csv('Logan_MonthlyTotalPrecip_2011_2016-in.csv')
temp <- read.csv('avg-temp.csv')
clim <- read.csv('logan-climate-data-combined.csv')

m.precip <- melt(precip,id.vars = 'Year')
m.precip.ann <- m.precip[m.precip$variable=='Annual',]
d$annprecip <- m.precip.ann[match(d$year,m.precip.ann$Year),"value"]
d$month <-  month.abb[month(as.POSIXlt(d$date, format="%Y-%m-%d"))]
d$moyr <- paste(d$month,d$year)
d$seasyr <- paste(d$season,d$year)
d$seasyr <- factor(d$seasyr, levels = c('Spring 2011','Fall 2011','Spring 2012',
                                        'Fall 2012','Spring 2013','Fall 2013',
                                        'Spring 2014','Fall 2014','Spring 2015',
                                        'Fall 2015','Spring 2016'))
m.precip$moyr <- paste(m.precip$variable,m.precip$Year)
d$moprecip <- m.precip[match(d$moyr,m.precip$moyr),"value"]
d$moprecip <- d$moprecip*2.54


clim$day <- str_pad(as.character(clim$day),2,pad="0")
#sapply(clim$month,function(x) grep(paste("?i",x,sep=""),month.abb))
clim$date <- paste0(clim$year,"-",clim$month,"-",clim$day)
clim$dateP <- as.Date(strptime(clim$date,format="%Y-%b-%d"))
#clim$dateP$zone <- NULL
table(complete.cases(clim))

clim$tempC <- (clim$temp-32)*(5/9)
clim$precip <- clim$precip.in*2.54
d$temp <- clim[match(d$date,clim$dateP),"tempC"]
d$precip <- clim[match(d$date,clim$dateP),"precip"]
d$snow_depth <- clim[match(d$date,clim$dateP),"snow_depth"]

clim$temp.1 <- c(NA,head(clim['tempC'],dim(clim)[1]-1)[[1]]) # shift down by one
clim$temp.3 <- c(rep(NA,2),rollmean(clim$tempC,3,align='right')) # 3-day moving average
clim$temp.5 <- c(rep(NA,4),rollmean(clim$tempC,5,align='right')) # 5-day moving average
clim$temp.10 <- c(rep(NA,9),rollmean(clim$tempC,10,align='right')) # 10-day moving average

clim$precip.1 <- c(NA,head(clim['precip'],dim(clim)[1]-1)[[1]]) # shift down by one
clim$precip.3 <- c(rep(NA,2),rollmean(clim$precip,3,align='right')) # 3-day moving average
clim$precip.5 <- c(rep(NA,4),rollmean(clim$precip,5,align='right')) # 5-day moving average
clim$precip.10 <- c(rep(NA,9),rollmean(clim$precip,10,align='right')) # 10-day moving average

clim$precip.sum.2 <- rollsumr(clim$precip,k=2,fill=NA)
clim$precip.sum.3 <- rollsumr(clim$precip,k=3,fill=NA)
clim$precip.sum.5 <- rollsumr(clim$precip,k=5,fill=NA)
clim$precip.sum.10 <- rollsumr(clim$precip,k=10,fill=NA)

d$temp.1 <- clim[match(d$date,clim$dateP),"temp.1"]
d$temp.3 <- clim[match(d$date,clim$dateP),"temp.3"]
d$temp.5 <- clim[match(d$date,clim$dateP),"temp.5"]
d$temp.10 <- clim[match(d$date,clim$dateP),"temp.10"]

d$precip.1 <- clim[match(d$date,clim$dateP),"precip.1"]
d$precip.3 <- clim[match(d$date,clim$dateP),"precip.3"]
d$precip.5 <- clim[match(d$date,clim$dateP),"precip.5"]
d$precip.10 <- clim[match(d$date,clim$dateP),"precip.10"]

d$precip.sum.2 <- clim[match(d$date,clim$dateP),"precip.sum.2"]
d$precip.sum.3 <- clim[match(d$date,clim$dateP),"precip.sum.3"]
d$precip.sum.5 <- clim[match(d$date,clim$dateP),"precip.sum.5"]
d$precip.sum.10 <- clim[match(d$date,clim$dateP),"precip.sum.10"]

###################################

table(d$species,d$sex)
table(d$species,d$sex,d$age)
table(d$species,d$sex,d$site)

d2 <- na.omit(d, cols=c("sex", "species"))
d3 <- droplevels(d[d$site != 'SH',])
d4 <- droplevels(d3[d3$species != 'UNK',])
d5 <- droplevels(d4[d4$sex == 'male',])

table <- as.data.frame(table(d5$seasyr,d5$site,d5$species))

summ <- d5 %>%
  group_by(seasyr,site,species) %>%
  summarise(mean.svl=round(mean(svl,na.rm=T),0),
            sd.svl = round(sd(svl,na.rm=T),0),
            se.svl = round(sd.svl/sqrt(length(svl)),0),
            ci95.svl = round(se.svl*1.96,0),
            N.svl = length(svl),
            mean.mass=round(mean(mass,na.rm=T),0),
            sd.mass = round(sd(mass,na.rm=T),0),
            se.mass = round(sd.mass/sqrt(length(mass)),0),
            ci95.mass = round(se.mass*1.96,0),
            N.mass = length(mass))
write.csv(summ,"summary-svl-mass.csv")

sir <- d3[d3$species == "SIR",]
ele <- d3[d3$species == "ELE",]

sirF <- droplevels(sir[sir$sex == 'female',])
sirM <- droplevels(sir[sir$sex == 'male',])
eleF <- droplevels(ele[ele$sex == 'female',])
eleM <- droplevels(ele[ele$sex == 'male',])

d4M <- droplevels(d4[d4$sex == 'male',])
d4F <- droplevels(d4[d4$sex == 'female',])

#### checks for normality (nothing is normal)

normal <- apply(d[,c(7:14,16:17)], 2, shapiro.test)
capture.output(normal,file="normal.txt")

normal2 <- apply(ele[,c(7:14,16:17)], 2, shapiro.test)
capture.output(normal2,file="normal_elegans.txt")

normal3 <- apply(sir[,c(7:14,16:17)], 2, shapiro.test)
capture.output(normal3,file="normal_sirtalis.txt")

# skewness should be close to 0 (-2 to 2), kurtosis close to 3

skewness(d$svl,na.rm = T)
kurtosis(d$svl,na.rm = T)
hist(d$svl)
scatterplot(d$svl,d$mass)

#### plots

ggplot(d,aes(x=date,y=bodycon,color=species))+geom_point()+facet_grid(species~sex)

ggplot(d,aes(x=date,y=cort_bl,color=species))+geom_point()+facet_grid(species~sex)

ggplot(d,aes(x=date,y=cort_bl,color=species))+geom_point()+facet_grid(year~species)

ggplot(d,aes(x=date,y=bka_bl))+geom_point()+facet_grid(species~sex)

ggplot(d,aes(x=date,y=test,color=species))+geom_point()+facet_grid(year~species)

ggplot(d,aes(x=date,y=cort_react,color=species))+geom_point()+facet_grid(species~sex)

#### stats

# summary(lm(cort_bl~date*sex,ele))
# summary(lm(cort_bl~date*sex,sir))

summary(lm(svl~site,eleM))
summary(lm(svl~site,sirM))

summary(lm(mass~site,eleM))
summary(lm(mass~site,sirM))

summary(lm(bodycon~site,eleM))
summary(lm(bodycon~site,sirM))

summary(lm(cort_bl~site,eleM))
summary(lm(cort_bl~site,sirM))

summary(lm(bka_bl~site,eleM))
summary(lm(bka_bl~site,sirM))

summary(lm(glu_bl~site,eleM))
summary(lm(glu_bl~site,sirM))

summary(lm(test~site,eleM))
summary(lm(test~site,sirM))

summary(lm(cort_ps~site,eleM))
summary(lm(cort_ps~site,sirM))

summary(lm(bka_ps~site,eleM))
summary(lm(bka_ps~site,sirM))

summary(lm(glu_ps~site,eleM))
summary(lm(glu_ps~site,sirM))

summary(lm(cort_react~site,eleM))
summary(lm(cort_react~site,sirM))

summary(lm(bka_react~site,eleM))
summary(lm(bka_react~site,sirM))

summary(lm(cort_bl~species,d4M))
summary(lm(bka_bl~species,d4M))
summary(lm(glu_bl~species,d4M))
summary(lm(test~species,d4M))
summary(lm(cort_ps~species,d4M))
summary(lm(bka_ps~species,d4M))
summary(lm(glu_ps~species,d4M))
summary(lm(cort_react~species,d4M))
summary(lm(bka_react~species,d4M))
summary(lm(svl~species,d4M))
summary(lm(mass~species,d4M))
summary(lm(bodycon~species,d4M))

# d4Mm <- melt(d4M,idvars=c('id','season','subsite','site','species','notes','age','sex',
#                           'month','moyr','date','year'))
# d4Mm <- d4Mm[d4Mm$variable %in% vars,]
#d4Mm$vars <- as.factor(d4Mm$vars,) # TODO: reorder factor levels as in vars

vars <- c('cort_bl','bka_bl','glu_bl','test','cort_ps','bka_ps','glu_ps','cort_react','bka_react','mass','svl','bodycon')

e <- reshape2::melt(d4M, id.vars = c('id','season','subsite','site','species','notes','age','sex',
                                     'month','moyr','date','year'))
e <- e[e$variable %in% vars,]

e$variable <- factor(e$variable,
                     levels = c('cort_bl','bka_bl','glu_bl','test','cort_ps','bka_ps','glu_ps','cort_react','bka_react','mass','svl','bodycon'),
                     labels = c('CORT (baseline)','BKA (baseline)','Glucose (baseline)','T','CORT (post-stress)','BKA (post-stress)','Glucose (post-stress)','CORT (??)','BKA (??)','Mass (g)','Snout-vent Length (mm)','Body Condition'))

e$season <- factor(e$season, levels = c('Spring','Fall'))

e$species <- factor(e$species, levels = c('ELE','SIR'), labels = c('T. elegans','T. sirtalis'))

e$value <- as.numeric(e$value)
  
ggplot(e,aes(x=value,fill=species))+
  facet_wrap(variable~.,scales="free")+
  geom_histogram()

ggplot(e,aes(y=value,group=species,fill=species))+
  facet_wrap(variable~.,scales="free")+
  geom_boxplot()+
  scale_fill_manual(values = c('#EBC08B','#BC3C1F'))+
  theme(axis.text.x =  element_blank(), axis.ticks.x =  element_blank(), 
        axis.title.y = element_blank(), legend.text = element_text(face='italic'))

eS <- droplevels(e[e$species=="SIR",])
eE <- droplevels(e[e$species=="ELE",])

ggplot(eS,aes(y=value,group=season,fill=season))+
  facet_wrap(variable~.,scales="free")+
  geom_boxplot()+
  scale_fill_manual(values = c('#a1d76a','#f1a340'))+
  theme(axis.text.x =  element_blank(), axis.ticks.x =  element_blank(), axis.title.y = element_blank())
  
ggplot(eE,aes(y=value,group=season,fill=season))+
  facet_wrap(variable~.,scales="free")+
  geom_boxplot()+
  scale_fill_manual(values = c('#a1d76a','#f1a340'))+
  theme(axis.text.x =  element_blank(), axis.ticks.x =  element_blank(), axis.title.y = element_blank())

# trade-offs?
summary(lm(cort_bl~bka_bl*sex,ele))
summary(lm(cort_bl~bka_ps*sex,ele))
summary(lm(cort_ps~bka_bl*sex,ele))
summary(lm(cort_ps~bka_ps*sex,ele))
summary(lm(cort_react~bka_bl*sex,ele))
summary(lm(cort_react~bka_ps*sex,ele))
summary(lm(test~bka_bl,ele))
summary(lm(test~cort_bl,ele))
summary(lm(test~cort_react,ele))
summary(lm(glu_bl~bka_bl,ele))
summary(lm(glu_bl~cort_bl,ele))
summary(lm(glu_bl~cort_react,ele))
summary(lm(glu_ps~bka_ps,ele))
summary(lm(glu_ps~cort_ps,ele))
summary(lm(glu_ps~cort_react,ele))

# trade-offs?
summary(lm(cort_bl~bka_bl*sex,sir))
summary(lm(cort_bl~bka_ps*sex,sir))
summary(lm(cort_ps~bka_bl*sex,sir))
summary(lm(cort_ps~bka_ps*sex,sir))
summary(lm(cort_react~bka_bl*sex,sir))
summary(lm(cort_react~bka_ps*sex,sir))        
summary(lm(test~bka_bl,sir))
summary(lm(test~cort_bl,sir))
summary(lm(test~cort_react,sir))
summary(lm(glu_bl~bka_bl,sir))
summary(lm(glu_bl~cort_bl,sir))
summary(lm(glu_bl~cort_react,sir))
summary(lm(glu_ps~bka_ps,sir))
summary(lm(glu_ps~cort_ps,sir))
summary(lm(glu_ps~cort_react,sir))

# correlations
contvars <- c('cort_bl','bka_bl','glu_bl','test','cort_ps','bka_ps','glu_ps','cort_react','bka_react','mass','svl','bodycon')
ele.c <- ele[contvars]
colnames(ele.c) <- c('CORT (baseline)','Bactericial ability (baseline)','Glucose (baseline)','T','CORT (post-stress)','Bactericial ability (post-stress)','Glucose (post-stress)','CORT (??)','Bactericial ability (??)','Mass (g)','Snout-vent Length (mm)','Body Condition')
res <- cor(ele.c, method = "pearson", use = "complete.obs")
round(res, 2)
corrplot(res, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)

sir.c <- sir[contvars]
colnames(sir.c) <- c('CORT (baseline)','Bactericial ability (baseline)','Glucose (baseline)','T','CORT (post-stress)','Bactericial ability (post-stress)','Glucose (post-stress)','CORT (??)','Bactericial ability (??)','Mass (g)','Snout-vent Length (mm)','Body Condition')
res2 <- cor(sir.c, method = "pearson", use = "complete.obs")
round(res2, 2)
corrplot(res2, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)

# weather data plots

summary(lm(bka_bl~moprecip*sex,ele))
summary(lm(bka_bl~moprecip*sex,sir))
ggplot(d,aes(x=moprecip,y=bka_bl,color=species))+geom_point()+facet_grid(species~sex)

summary(lm(cort_bl~moprecip*sex,ele))
summary(lm(cort_bl~moprecip*sex,sir))
ggplot(d,aes(x=moprecip,y=cort_bl,color=species))+geom_point()+facet_grid(species~sex)

# summary(lm(cort_react~moprecip*sex,ele))
# summary(lm(cort_react~moprecip*sex,sir))
# ggplot(d,aes(x=moprecip,y=cort_react,color=species))+geom_point()

summary(lm(test~moprecip,ele))
summary(lm(test~moprecip,sir))
ggplot(d,aes(x=moprecip,y=test))+geom_point()

summary(lm(glu_bl~moprecip*sex,ele))
summary(lm(glu_bl~moprecip*sex,sir))
ggplot(d,aes(x=moprecip,y=glu_bl,color=species))+geom_point()+facet_grid(species~sex)


my_y_title <- expression(paste(italic("T. elegans")," bactericidal ability (%; baseline)"))
summary(lm(bka_bl~moprecip,eleM))
ggplot(eleM,aes(x=moprecip,y=bka_bl))+geom_point()+geom_smooth(method="lm")+
  labs(y=my_y_title,x="Monthly precipitation (cm)")
summary(lm(bka_bl~moprecip,sirM))
ggplot(sirM,aes(x=moprecip,y=bka_bl))+geom_point()

## weather data moving averages

ggplot(eleM,aes(x=temp.10,y=bka_bl))+geom_point()
ggplot(sirM,aes(x=precip.10,y=test))+geom_point()
ggplot(eleM,aes(x=temp.10,y=test))+geom_point()
ggplot(eleM,aes(x=precip.5,y=test))+geom_point()
ggplot(eleM,aes(x=precip.sum.5,y=test))+geom_point()
ggplot(sirM,aes(x=precip.sum.10,y=test))+geom_point()
ggplot(sirM,aes(x=precip.10,y=glu_bl))+geom_point()
ggplot(sirM,aes(x=precip.10,y=glu_ps))+geom_point()

my_y_title_2 <- expression(paste(italic("T. sirtalis")," CORT (ng/mL; baseline)"))
summary(lm(cort_bl~temp.10,sirM))
summary(lm(cort_bl~temp.10*site,sirM))
summary(lm(cort_bl~temp.10*season,sirM))
summary(lm(cort_bl~temp.10*season*site,sirM))
ggplot(sirM,aes(x=temp.10,y=cort_bl))+geom_point()+geom_smooth(method="lm")+
  labs(y=my_y_title_2,x="Mean Temperature (°C; past 10 days)")+xlim(2.5,15)

my_y_title_3 <- expression(paste(italic("T. sirtalis")," testosterone (ng/mL; baseline)"))
summary(lm(test~precip.sum.10,sirM))
summary(lm(test~precip.sum.10*site,sirM))
summary(lm(test~precip.sum.10*season,sirM))
summary(lm(test~precip.sum.10*season*site,sirM))
ggplot(sirM,aes(x=precip.sum.10,y=test,color=season))+geom_point()+geom_smooth(method="lm")+
  xlim(0,2.6)+
  scale_color_manual(breaks = c("Spring","Fall"),values = c('#f1a340','#a1d76a'))+
  labs(y=my_y_title_3,x="Cumulative precipitation (cm; past 10 days)")
ggplot(sirM,aes(x=precip.sum.10,y=test,color=season))+geom_point()+geom_smooth(method="lm")+
  xlim(0,3)+facet_grid(~season)+
  scale_color_manual(values = c('#a1d76a','#f1a340'))+
  labs(y=my_y_title_3,x="Cumulative precipitation (cm; past 10 days)")

summary(lm(glu_ps~temp*season,sirM))
summary(lm(glu_ps~temp.1*season,sirM))
summary(lm(glu_ps~temp.3*season,sirM))
summary(lm(glu_ps~temp.5*season,sirM))
summary(lm(glu_ps~temp.10*season,sirM))

summary(lm(glu_ps~precip*season,sirM))
summary(lm(glu_ps~precip.1*season,sirM))
summary(lm(glu_ps~precip.3*season,sirM))
summary(lm(glu_ps~precip.5*season,sirM))
summary(lm(glu_ps~precip.10*season,sirM))

summary(lm(glu_ps~precip.sum.2*season,sirM))
summary(lm(glu_ps~precip.sum.3*season,sirM))
summary(lm(glu_ps~precip.sum.5*season,sirM))
summary(lm(glu_ps~precip.sum.10*season,sirM))

summary(lm(glu_ps~temp*season,eleM))
summary(lm(glu_ps~temp.1*season,eleM))
summary(lm(glu_ps~temp.3*season,eleM))
summary(lm(glu_ps~temp.5*season,eleM))
summary(lm(glu_ps~temp.10*season,eleM))

summary(lm(glu_ps~precip*season,eleM))
summary(lm(glu_ps~precip.1*season,eleM))
summary(lm(glu_ps~precip.3*season,eleM))
summary(lm(glu_ps~precip.5*season,eleM))
summary(lm(glu_ps~precip.10*season,eleM))

summary(lm(glu_ps~precip.sum.2*season,eleM))
summary(lm(glu_ps~precip.sum.3*season,eleM))
summary(lm(glu_ps~precip.sum.5*season,eleM))
summary(lm(glu_ps~precip.sum.10*season,eleM))

## seasonality

summary(lm(svl~site*season,eleM))
summary(lm(svl~site*season,sirM))

summary(lm(mass~site*season,eleM))
summary(lm(mass~site*season,sirM))

summary(lm(bodycon~site*season,eleM))
summary(lm(bodycon~site*season,sirM))

summary(lm(cort_bl~site*season,eleM))
summary(lm(cort_bl~site*season,sirM))

summary(lm(bka_bl~site*season,eleM))
summary(lm(bka_bl~site*season,sirM))

summary(lm(glu_bl~site*season,eleM))
summary(lm(glu_bl~site*season,sirM))

summary(lm(test~site*season,eleM))
summary(lm(test~site*season,sirM))

summary(lm(cort_ps~site*season,eleM))
summary(lm(cort_ps~site*season,sirM))

summary(lm(bka_ps~site*season,eleM))
summary(lm(bka_ps~site*season,sirM))

summary(lm(glu_ps~site*season,eleM))
summary(lm(glu_ps~site*season,sirM))

summary(lm(cort_react~site*season,eleM))
summary(lm(cort_react~site*season,sirM))

summary(lm(bka_react~site*season,eleM))
summary(lm(bka_react~site*season,sirM))

## season & year

summary(lm(svl~site*season,eleM))
summary(lm(mass~site*season,eleM))
summary(lm(bodycon~site*season,eleM))
summary(lm(cort_bl~site*season,eleM))
summary(lm(bka_bl~site*season,eleM))
summary(lm(glu_bl~site*season,eleM))
summary(lm(test~site*season,eleM))
summary(lm(cort_ps~site*season,eleM))
summary(lm(bka_ps~site*season,eleM))
summary(lm(glu_ps~site*season,eleM))
summary(lm(cort_react~site*season,eleM))
summary(lm(bka_react~site*season,eleM))

## manova

op <- options(contrasts = c("contr.helmert", "contr.poly"))

npk2.aov <- manova(cbind(cort_bl,bka_bl,glu_bl,test,cort_ps,bka_ps,glu_ps,cort_react,bka_react,mass,svl,bodycon) ~ species * site * season, d4M)
summary(npk2.aov)



