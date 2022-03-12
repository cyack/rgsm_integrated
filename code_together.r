####################
#### STEP 1: Read in all data
#####################



t2asir<-read.csv("ASIRQuery4Combined.csv")
tdry<-read.csv("rivereyes_v2.csv")
angQ<-read.csv("abq_gage_08330000.csv")
sanaQ<-read.csv("SanAcacia_gage_08354900.csv")
t2vie<-read.csv("releases 8_26_2019.csv")
converter<-read.csv("fws2br_rmconverter.csv")
tmeso<-read.csv("mesohab_sum.csv")
trescue<-read.csv("Fish Rescue 2009-2020.csv")
phiR<-read.csv("rescue_surv.csv")
ee1<-read.csv("sum_ee1.csv")
ee2<-read.csv("sum_ee2.csv")
ee3<-read.csv("sum_ee3.csv")
ee4<-read.csv("sum_ee4.csv")
ee5<-read.csv("sum_ee5.csv")
aQ<-read.csv("aQ.csv",header=T)[,c(3:19)]
sQ<-read.csv("sQ.csv",header=T)[,c(3:19)]
##out of sample data
oosD<-read.csv("oos_rivereyes.csv")
oosR<-read.csv("oos_releases.csv")
oosc<-read.csv("oos_catch.csv") 



######################
#### STEP 2: Define variables and reformat data
######################



# define two grids
grid<-c(55.4,116,170,210.1) # start and end points for 3 main river segments
vgrid<-seq(55.4,210.2,0.2) # represents all 200 meter segments in study area
crossgrid<-findInterval(vgrid[-length(vgrid)]+.05,grid)
gridst_end<-rbind(range(which(crossgrid==1)),range(which(crossgrid==2)),range(which(crossgrid==3)))
years<-c(2002:2018) ##years used to fit model
Nyears<-length(years)
# reformat drying data
tdry$ustrata<-findInterval(tdry$URM,vgrid)
tdry$dstrata<-findInterval(tdry$LRM,vgrid)
tdry$jul<-NA
for (i in 1:length(tdry$jul)){
	tdry$jul[i]<-julian(as.Date(paste(tdry$Year[i],tdry$month[i],tdry$dom[i],sep="-")),as.Date(paste(tdry$Year[i],"03","31",sep="-")))[[1]]}
fdryday<-matrix(NA,nrow=(Nyears),ncol=(length(vgrid)-1))
ldryday<-matrix(NA,nrow=(Nyears),ncol=(length(vgrid)-1))
for (k in 1:(length(vgrid)-1)){
	for (j in 1:Nyears){
		temp2<-subset(tdry$jul,k<=tdry$ustrata&k>=tdry$dstrata&tdry$Year==(2001+j))
		fdryday[j,k]<-ifelse(length(temp2)==0,NA,min(temp2))
		ldryday[j,k]<-ifelse(length(temp2)==0,NA,max(temp2))
	}}
# reformat augmented fish release data
vkeep<-match(c("Date","Month","Year","C","L","S","Number","RM_up","RM_down"),names(t2vie))
tvie<-t2vie[,vkeep]
tvie$strata_up<-findInterval(converter[match(tvie$RM_up,converter[,1]),2],vgrid)
tvie$strata_down<-findInterval(converter[match(tvie$RM_down,converter[,1]),2],vgrid)
tlen<-length(tvie$strata_down)
vie<-tvie
for (i in 1:tlen){
	temp<-tvie$strata_up[i]-tvie$strata_down[i]
	vie$Number[i]<-tvie$Number[i]/(1+temp)
	if (temp>0){
		for (j in 1:temp){
			temp2<-tvie[i,]
			temp2$strata_down<-temp2$strata_down+j
			temp2$Number<-tvie$Number[i]/(1+temp)
			vie<-rbind(vie,temp2)
			}}}
tvie<-vie
tvie$C<-ifelse(tvie$Year>2007|(tvie$Year==2007&tvie$Month>9),paste(tvie$C),"")
tvie$L<-ifelse(tvie$Year>2007|(tvie$Year==2007&tvie$Month>9),paste(tvie$L),"")
tvie$S<-ifelse(tvie$Year>2007|(tvie$Year==2007&tvie$Month>9),paste(tvie$S),"")
tvie$VIE<-paste(tvie$C,tvie$L,sep="")
Vcodes<-sort(unique(tvie$VIE))
tvie$Vc<-match(tvie$VIE,Vcodes)
tvie$period<-ifelse(tvie$Month<4,8*(tvie$Year-2002)+1,ifelse(tvie$Month>9,8*(tvie$Year-2002)+9,8*(tvie$Year-2002)+tvie$Month-3))
tvie$winmon<-ifelse(tvie$Month<4,4-tvie$Month,ifelse(tvie$Month>9,16-tvie$Month,0))
vkeep<-match(c("Number","strata_down","Vc","period","winmon"),names(tvie))
vie<-tvie[,vkeep]
vie<-subset(vie,vie$period!=130&vie$period!=132) ###complicate things too much and not that many fish
rm(tvie,vkeep)
####
Nstrata<-3
Nvstrata<-length(vgrid)-1
nw_um<-matrix(0,nrow=30,ncol=Nstrata) #summarized at coarse grid
w_um<-matrix(0,nrow=6,ncol=Nstrata)
nmons<-matrix(0,nrow=6,ncol=Nstrata)
t1<-subset(vie,vie$Vc==1&8*(vie$period/8-floor(vie$period/8))==1)
t2<-subset(vie,vie$Vc==1&8*(vie$period/8-floor(vie$period/8))!=1)
t3<-subset(vie,vie$Vc>1)
tw_m<-array(0,dim=c((Nyears-5),Nvstrata,2)) #summarized at fine grid
w_m<-array(0,dim=c((Nyears-5),Nstrata,2)) #summarized at coarse grid
for (i in 1:6){
	for (j in 1:Nstrata){
		temp<-subset(t1,((t1$period-1)/8+1)==i&crossgrid[t1$strata_down]==j)
		if (dim(temp)[1]>0){
			w_um[i,j]<-sum(temp[,1])
			nmons[i,j]<-round(sum(temp[,1]*temp[,5])/sum(temp[,1]))
			}}}
for (i in 1:6){
	for (t in 1:5){
		for (j in 1:Nstrata){
			temp<-subset(t2,(floor(t2$period/8)+1)==i&(8*(t2$period/8-floor(t2$period/8))-1)==t&crossgrid[t2$strata_down]==j)
			if (dim(temp)[1]>0){
			nw_um[((i-1)*5+t),j]<-sum(temp[,1])
			}}}}
for (i in 1:(Nyears-5)){
	for (j in 1:Nvstrata){
		temp<-subset(t3,((t3$period-1)/8-5)==i&t3$strata_down==j)
		if (dim(temp)[1]>0){
			tw_m[i,j,1]<-sum(subset(temp[,1],temp$winmon<3))
			tw_m[i,j,2]<-sum(subset(temp[,1],temp$winmon>3))
			}}
	for (j in 1:Nstrata){
		w_m[i,j,1]<-sum(tw_m[i,gridst_end[j,1]:gridst_end[j,2],1])
		w_m[i,j,2]<-sum(tw_m[i,gridst_end[j,1]:gridst_end[j,2],2])
		}}
w_um[6,3]<-sum(nw_um[25:26,3])
##### done reformating vie release data	

#simplify and reformat monitoring database
keep<-match(c("ProjectName","DateSampled","SiteID","year", "month","RMStart","HabitatNumber","Habitat","SamplingEffort","RepeatedSamplingNumber",
"Species","DepletionNumber","Gear","NumberCaptured","AgeClass","LengthSL","LengthMinSL","LengthMaxSL","VIEColor","VIELocation"),names(t2asir))
tasir<-t2asir[,keep]
tasir<-subset(tasir,tasir$Habitat!=""&tasir$Gear!="larval"&tasir$year>2001)
ghab<-c("","BW","MCPLPO","MCPO","MCSHPLPO","MCSHPO","PO","SCPLPO","SCPO","SCSHPLPO","SCSHPO","SHPO","MCED","MCSHED","SCED","SCSHED") # sampled pool habitats
ghabAv<-c("BW NS","PO NS","SHPO NS") #not sampled pool habitats - only measured during population estimation
bhabAv<-c("RU NS","SHRU NS") #not sampled run/riffle habitat.
tasir$cHab<-ifelse(is.na(match(tasir$Habitat,ghab))==FALSE,1,ifelse(is.na(match(tasir$Habitat,ghabAv))==FALSE,2,
	ifelse(is.na(match(tasir$Habitat,bhabAv))==FALSE,4,3)))
tasir$strata<-findInterval(tasir$RMStart,vgrid)
tasir$period<-8*(tasir$year-2002)+tasir$month-3
tasir<-subset(tasir,tasir$period>0)
extractdom<-function(x){
	lx<-nchar(x)
	tx<-numeric()
	for (i in 1:lx){
		tx[i]<-substr(x,i,i)}
	wx<-which(tx=="/")
	as.numeric(substr(x,wx[1]+1,wx[2]-1))}
tasir$dom<-NA
tasir$jul<-NA
tasir$dry<-NA
for (j in 1:length(tasir[,1])){
	tasir$dom[j]<-extractdom(paste(tasir$DateSampled[j]))
	tasir$jul[j]<-julian(as.Date(paste(tasir$year[j],tasir$month[j],tasir$dom[j],sep="-")),as.Date(paste(tasir$year[j],"03","31",sep="-")))[[1]]
	tf1<-fdryday[(tasir$year[j]-2001),tasir$strata[j]]
	tl1<-ldryday[(tasir$year[j]-2001),tasir$strata[j]]
	tasir$dry[j]<-ifelse(is.na(tf1)==TRUE|(tf1>tasir$jul[j])|(tl1<tasir$jul[j]),0,1)}
# subset different asir data for more reformating
###start with regular monitoring data
t2mon<-subset(tasir,tasir$ProjectName=="Hybognathus Amarus Population Monitoring"&tasir$month>3&tasir$month<11&tasir$cHab!=4&tasir$SamplingEffort!=""&tasir$SamplingEffort!="#N/A"&tasir$dry==0)
t2mon$uni_id<-paste(t2mon$year,t2mon$jul,t2mon$RMStart,t2mon$HabitatNumber)
t2mon$C<-ifelse(t2mon$year>2007,substr(paste(t2mon$VIEColor),1,1),"")
t2mon$L<-ifelse(t2mon$year>2007,substr(paste(t2mon$VIELocation),1,1),"")
t2mon$VC<-match(paste(t2mon$C,t2mon$L,sep=""),Vcodes)
t2mon$VC[which(t2mon$AgeClass==0&t2mon$VC>1)]<-1
# summarize data to seine haul
tmon<-data.frame(uni=sort(unique(t2mon$uni_id)))
tmon$RMStart<-t2mon$RMStart[match(tmon[,1],t2mon$uni_id)]
tmon$period<-t2mon$period[match(tmon[,1],t2mon$uni_id)]
tmon$jul<-t2mon$jul[match(tmon[,1],t2mon$uni_id)]
tmon$year<-t2mon$year[match(tmon[,1],t2mon$uni_id)]
tmon$month<-t2mon$month[match(tmon[,1],t2mon$uni_id)]
tmon$date<-t2mon$DateSampled[match(tmon[,1],t2mon$uni_id)]
tmon$angQ<-angQ$cfs[match(tmon$date,angQ$Date)]
tmon$sanaQ<-sanaQ$cfs[match(tmon$date,sanaQ$Date)]
tmon$HaulNo<-t2mon$HabitatNumber[match(tmon[,1],t2mon$uni_id)]
tmon$cHab<-t2mon$cHab[match(tmon[,1],t2mon$uni_id)]
tmon$effort<-as.numeric(paste(t2mon$SamplingEffort[match(tmon[,1],t2mon$uni_id)]))
tmon$strata<-t2mon$strata[match(tmon[,1],t2mon$uni_id)]
tmon$cQ<-ifelse(tmon$strata>573,tmon$angQ,tmon$sanaQ)
tmon$hybama0<-0
tmon$hybama1<-0
tmon$hybama2<-0
tmon$hybama3<-0
tmon$hybama4<-0
tmon$hybama5<-0
tmon$hybama6<-0
tmon$hybama7<-0
tmon$hybama8<-0
tmon$hybama9<-0
tmon$hybama10<-0
tmon$hybama11<-0
tmon$hybama12<-0
tmon$hybama13<-0
tmon$hybama14<-0
tmon$hybama01<-0
tstart<-which(names(tmon)=="hybama0")
t3mon<-subset(t2mon,t2mon$Species=="HYBAMA")
for (i in 1:length(t3mon[,1])){
	tx<-match(t3mon$uni_id[i],tmon$uni)
	ty<-ifelse(is.na(t3mon$AgeClass[i])==TRUE,1+tstart+length(Vcodes),ifelse(
		t3mon$AgeClass[i]==0,tstart,t3mon$VC[i]+tstart))
	tz<-t3mon$NumberCaptured[i]
	tmon[tx,ty]<-tmon[tx,ty]+tz
	}
tmon$type<-ifelse(tmon$cHab==1,1,2)
tmon$SPt_id<-paste(tmon$year,tmon$jul,tmon$strata,tmon$type)
tmon<-subset(tmon,tmon$cQ<1000) #remove seine hauls were discharge was greater than 1000 cfs
####subset data to summarize all hauls on the same day, in the same habitat and same river segment
mon<-data.frame(spt=sort(unique(tmon$SPt_id)))
mon$period<-tmon$period[match(mon[,1],tmon$SPt_id)]
mon$jul<-tmon$jul[match(mon[,1],tmon$SPt_id)]
mon$year<-tmon$year[match(mon[,1],tmon$SPt_id)]
mon$month<-tmon$month[match(mon[,1],tmon$SPt_id)]
mon$strata<-tmon$strata[match(mon[,1],tmon$SPt_id)]
mon$type<-tmon$type[match(mon[,1],tmon$SPt_id)]
mon$effort<-0
mon$cQ<-0
mon$hybama0<-0
mon$hybama1<-0
mon$hybama2<-0
mon$hybama3<-0
mon$hybama4<-0
mon$hybama5<-0
mon$hybama6<-0
mon$hybama7<-0
mon$hybama8<-0
mon$hybama9<-0
mon$hybama10<-0
mon$hybama11<-0
mon$hybama12<-0
mon$hybama13<-0
mon$hybama14<-0
mon$hybama01<-0
mean_names<-c("cQ")
sum_names<-c("effort","hybama0","hybama1","hybama2","hybama3","hybama4","hybama5",
	"hybama6","hybama7","hybama8","hybama9","hybama10","hybama11","hybama12",
	"hybama13","hybama14","hybama01")
tsm<-match(sum_names,names(tmon))
sm<-match(sum_names,names(mon))
tmn<-match(mean_names,names(tmon))
smn<-match(mean_names,names(mon))
for (i in 1:length(mon[,1])){
	temp<-subset(tmon[,tsm],tmon$SPt_id==mon$spt[i])
	mon[i,sm]<-colSums(temp)
	temp<-subset(tmon[,tmn],tmon$SPt_id==mon$spt[i])
	mon[i,smn]<-mean(temp)
	}
vie_names<-c("hybama2","hybama3","hybama4","hybama5","hybama6","hybama7","hybama8",
	"hybama9","hybama10","hybama11","hybama12","hybama13","hybama14")
mon$hybama_vie<-rowSums(mon[,vie_names])
mon$c_id<-paste(mon$year,mon$jul,crossgrid[mon$strata],mon$type)

## fit vie dispersal analysis described in appendix S1 to determine weights 
wm<-tw_m[,,1]+tw_m[,,2]
wm2<-tapply(wm[1,],crossgrid,sum)
for (j in 2:12){
	wm2<-rbind(wm2,tapply(wm[j,],crossgrid,sum))}
tt<-subset(mon,mon$year>2007)
tt$temp<-NA
for (i in 1:length(tt$temp)){
	tt$temp[i]<-wm2[(tt$year[i]-2007),crossgrid[tt$strata[i]]]}
tt2<-subset(tt,tt$temp!=0&tt$jul<60)
fitcauchy<-function(par){ # in units of 200 m sites
	a<-exp(par[2:3])
	predV<-numeric()	
	for (i in 1:401){
		t2<-ifelse(crossgrid==crossgrid[tt2$strata[i]],1,0)
		t3<-wm[(tt2$year[i]-2007),]*t2
		t4<-which(t3!=0)
		t4b<-which(t2!=0)
		t5<-numeric()
		for (j in 1:length(t4)){t5[j]<-sum(dcauchy(t4b,t4[j]+par[5],exp(par[1])))}
		t6<-sum(t3[t4]*dcauchy(tt2$strata[i],t4+par[5],exp(par[1]))/t5)
		predV[i]<-t6*tt2$effort[i]*a[tt2$type[i]]/400000}
	-1*sum(dnbinom(tt2$hybama_vie,mu=predV,size=exp(par[4]),log=TRUE))}
m<-optim(c(3,2,0,-2,0),fitcauchy,method="BFGS",hessian=TRUE)
mweights<-matrix(NA,nrow=12,ncol=Nvstrata)
for (i in 1:12){
	for (s in 1:3){
		t3<-wm[i,]*ifelse(crossgrid==s,1,0)
		t4<-which(t3!=0)
		if (length(t4)>0){
		t5<-numeric()
		for (j in 1:length(t4)){t5[j]<-sum(dcauchy(c(gridst_end[s,1]:gridst_end[s,2]),t4[j]+m$par[5],exp(m$par[1])))}
	t6<-numeric()
	for (k in gridst_end[s,1]:gridst_end[s,2]){t6[(k+1-gridst_end[s,1])]<-sum(t3[t4]*dcauchy(k,t4+m$par[5],exp(m$par[1]))/t5)}
	mweights[i,gridst_end[s,1]:gridst_end[s,2]]<-t6/sum(t6)
	}}}
mon$vieW<-0
for (i in 1:length(mon[,1])){
	mon$vieW[i]<-ifelse(mon$year[i]<2008,NA,mweights[(mon$year[i]-2007),mon$strata[i]])}
mon$Cstrata<-crossgrid[mon$strata]
monC<-data.frame(cid=sort(unique(mon$c_id)))
monC$jul<-mon$jul[match(monC[,1],mon$c_id)]
monC$year<-mon$year[match(monC[,1],mon$c_id)]
monC$month<-mon$month[match(monC[,1],mon$c_id)]
monC$Cstrata<-crossgrid[mon$strata[match(monC[,1],mon$c_id)]]
monC$type<-mon$type[match(monC[,1],mon$c_id)]
monC$effort<-0
monC$cQ<-0
monC$hybama0<-0
monC$hybama1<-0
monC$hybama01<-0
monC$hybama_vie<-0
monC$vieW<-0
monC$period<-mon$period[match(monC[,1],mon$c_id)]

mean_names<-c("cQ")
sum_names<-c("effort","hybama0","hybama1","hybama_vie","hybama01","vieW")
tsm<-match(sum_names,names(mon))
sm<-match(sum_names,names(monC))
tmn<-match(mean_names,names(mon))
smn<-match(mean_names,names(monC))
for (i in 1:length(monC[,1])){
	temp<-subset(mon[,tsm],mon$c_id==monC$cid[i])
	monC[i,sm]<-colSums(temp)
	temp<-subset(mon[,tmn],mon$c_id==monC$cid[i])
	monC[i,smn]<-mean(temp)
	}
#### Done reformating monthly April to October monitoring data
#### reformat November data
t2monre<-subset(tasir,tasir$ProjectName=="Hybognathus Amarus Population Monitoring Repeated"&tasir$VIEColor==""&tasir$VIELocation==""&tasir$dry==0)
t2monre$uni_id<-paste(t2monre$year,t2monre$jul,t2monre$RMStart,t2monre$HabitatNumber)
tmonre<-data.frame(uni_YSH=sort(unique(t2monre$uni_id)))
tmonre$Date<-t2monre$Date[match(tmonre[,1],t2monre$uni_id)]
tmonre$jul<-t2monre$jul[match(tmonre[,1],t2monre$uni_id)]
tmonre$angQ<-angQ$cfs[match(tmonre$Date,angQ$Date)]
tmonre$sanaQ<-sanaQ$cfs[match(tmonre$Date,sanaQ$Date)]
tmonre$RMStart<-t2monre$RMStart[match(tmonre[,1],t2monre$uni_id)]
tmonre$year<-t2monre$year[match(tmonre[,1],t2monre$uni_id)]
tmonre$HabNo<-t2monre$HabitatNumber[match(tmonre[,1],t2monre$uni_id)]
tmonre$cHab<-t2monre$cHab[match(tmonre[,1],t2monre$uni_id)]
tmonre$effort<-as.numeric(paste(t2monre$SamplingEffort[match(tmonre[,1],t2monre$uni_id)]))
tmonre$Cstrata<-crossgrid[t2monre$strata[match(tmonre[,1],t2monre$uni_id)]]
tmonre$cQ<-ifelse(tmonre$Cstrata==3,tmonre$angQ,tmonre$sanaQ)
tmonre$hybama_1<-0
tmonre$hybama_2<-0
tmonre$hybama_3<-0
tmonre$hybama_4<-0

tstart<-which(names(tmonre)=="hybama_1")
t3monre<-subset(t2monre,t2monre$Species=="HYBAMA")
for (i in 1:length(t3monre[,1])){
	tx<-match(t3monre$uni_id[i],tmonre$uni_YSH)
	ty<-t3monre$RepeatedSamplingNumber[i]+tstart-1
	tz<-t3monre$NumberCaptured[i]
	tmonre[tx,ty]<-tmonre[tx,ty]+tz
	}
tmonre$type<-ifelse(tmonre$cHab==1,1,2)
tmonre$period<-t2monre$period[match(tmonre[,1],t2monre$uni_id)]
tmonre$SPt_id<-paste(tmonre$year,tmonre$jul,tmonre$Cstrata,tmonre$type)
## summarize by day, river segement and habitat type
monre<-data.frame(spt=sort(unique(tmonre$SPt_id)))
monre$period<-tmonre$period[match(monre[,1],tmonre$SPt_id)]
monre$year<-tmonre$year[match(monre[,1],tmonre$SPt_id)]
monre$month<-tmonre$month[match(monre[,1],tmonre$SPt_id)]
monre$Cstrata<-tmonre$Cstrata[match(monre[,1],tmonre$SPt_id)]
monre$type<-tmonre$type[match(monre[,1],tmonre$SPt_id)]
monre$jul<-tmonre$jul[match(monre[,1],tmonre$SPt_id)]
monre$effort<-0
monre$cQ<-0
monre$hybama_1<-0
mean_names<-c("cQ")
sum_names<-c("effort","hybama_1")
tsm<-match(sum_names,names(tmonre))
sm<-match(sum_names,names(monre))
tmn<-match(mean_names,names(tmonre))
smn<-match(mean_names,names(monre))
for (i in 1:length(monre[,1])){
	temp<-subset(tmonre[,tsm],tmonre$SPt_id==monre$spt[i])
	monre[i,sm]<-colSums(temp)
	temp<-subset(tmonre[,tmn],tmonre$SPt_id==monre$spt[i])
	monre[i,smn]<-mean(temp)
	}
monre<-subset(monre,monre$cQ<1000)
#### only looking at these data for habitat availabilityestimates
t2est<-subset(tasir,tasir$ProjectName=="Hybognathus Amarus Population Estimation")
t2est$uni_YSH<-paste(t2est$year,t2est$RMStart,t2est$HabitatNumber)
test<-data.frame(uni_YSH=unique(t2est$uni_YSH))
test$RMStart<-t2est$RMStart[match(test[,1],t2est$uni_YSH)]
test$year<-t2est$year[match(test[,1],t2est$uni_YSH)]
test$Date<-t2est$Date[match(test[,1],t2est$uni_YSH)]
test$angQ<-angQ$cfs[match(test$Date,angQ$Date)]
test$sanaQ<-sanaQ$cfs[match(test$Date,sanaQ$Date)]
test$HabNo<-t2est$HabitatNumber[match(test[,1],t2est$uni_YSH)]
test$cHab<-t2est$cHab[match(test[,1],t2est$uni_YSH)]
test$effort<-as.numeric(paste(t2est$SamplingEffort[match(test[,1],t2est$uni_YSH)]))
test$Cstrata<-crossgrid[t2est$strata[match(test[,1],t2est$uni_YSH)]]
test$cQ<-ifelse(test$Cstrata==3,test$angQ,test$sanaQ)
#
test$uni_YS<-paste(test$year,test$RMStart)
av<-data.frame(uni_YS=unique(test$uni_YS))
av$RMStart<-test$RMStart[match(av[,1],test$uni_YS)]
av$Date<-test$Date[match(av[,1],test$uni_YS)]
av$angQ<-angQ$cfs[match(av$Date,angQ$Date)]
av$sanaQ<-sanaQ$cfs[match(av$Date,sanaQ$Date)]
av$year<-test$year[match(av[,1],test$uni_YS)]
av$Cstrata<-test$Cstrata[match(av[,1],test$uni_YS)]
av$cQ<-ifelse(av$Cstrata==3,av$angQ,av$sanaQ)
av$gHab<-0
av$bHab<-0
for (i in 1:length(av[,1])){
	temp<-subset(test,test$uni_YS==av[i,1])
	av$gHab[i]<-sum(subset(temp$effort,temp$cHab<3))
	av$bHab[i]<-sum(subset(temp$effort,temp$cHab>2))
	}
av$TotHab<-av$gHab+av$bHab
av$pro<-av$gHab/av$TotHab
av<-subset(av,av$year>2008) ### only use 2009-2011 because 2008 data is incomplete
#### reformat data from Braun et al. 2015
tmeso$Chab<-ifelse(tmeso$mesohab_cl==1|tmeso$mesohab_cl==7|tmeso$mesohab_cl==8|tmeso$mesohab_cl==14|tmeso$mesohab_cl==15,1,3)
tmeso$Cstrata<-findInterval(tmeso$RM,grid)
tmeso$angQ<-angQ$cfs[match(tmeso$Date,angQ$Date)]
tmeso$sanaQ<-sanaQ$cfs[match(tmeso$Date,sanaQ$Date)]
tmeso$cQ<-ifelse(tmeso$Cstrata==3,tmeso$angQ,tmeso$sanaQ)
tmeso$uni_RD<-paste(tmeso$RM,tmeso$Date)
meso<-data.frame(uni_RD=unique(tmeso$uni_RD))
meso$RM<-tmeso$RM[match(meso[,1],tmeso$uni_RD)]
meso$site_abv<-tmeso$site_abv[match(meso[,1],tmeso$uni_RD)]
meso$Date<-tmeso$Date[match(meso[,1],tmeso$uni_RD)]
meso$Cstrata<-tmeso$Cstrata[match(meso[,1],tmeso$uni_RD)]
meso$cQ<-tmeso$cQ[match(meso[,1],tmeso$uni_RD)]
meso$gHab<-0
meso$bHab<-0
for (j in 1:length(meso[,1])){
	temp<-subset(tmeso,tmeso$uni_RD==meso$uni_RD[j])
	meso$gHab[j]<-sum(subset(temp$Shape_Area,temp$Chab==1))
	meso$bHab[j]<-sum(subset(temp$Shape_Area,temp$Chab==3))
	}
meso$TotHab<-meso$gHab+meso$bHab
meso$pro<-meso$gHab/meso$TotHab
meso<-subset(meso,meso$Cstrata!=4)
mAV<-rbind(cbind(meso$Cstrata,meso$cQ/1000,meso$gHab,meso$gHab+meso$bHab),
	cbind(av$Cstrata,av$cQ/1000,av$gHab,av$gHab+av$bHab)) 
NAVsamps<-dim(mAV)[1]
#### done reformating and combining mesohabitat availability data

#### further summarize catch data 
# isolate data for age - 1+ fish only
tt<-subset(monC,(monC$hybama01==0&monC$month<10)|monC$month<7)
mon1<-cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,ifelse(tt$month<6,tt$hybama1+tt$hybama01,tt$hybama1))
mon1_effQ<-cbind(tt$effort,tt$cQ/1000)
colnames(mon1)<-c("year","julian","Cstrata","habitat","catch")
Nobs_mon1<-dim(mon1)[1]

# isolate data for augmented fish only
###lump different vie marks and start at period 49 for coarse strata 1 and 2 and 97 for coarse strata 3...
tt<-subset(mon,mon$Cstrata<3&mon$period>48&mon$jul<60)
monV<-cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama_vie)
monV_effQ<-cbind(tt$effort,tt$cQ/1000,tt$vieW)
tt<-subset(mon,mon$Cstrata>2&mon$period>96&mon$jul<60)
monV<-rbind(monV,cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama_vie))
monV_effQ<-rbind(monV_effQ,cbind(tt$effort,tt$cQ/1000,tt$vieW))
monV<-subset(monV,is.na(monV_effQ[,3])==FALSE)
monV_effQ<-subset(monV_effQ,is.na(monV_effQ[,3])==FALSE)
Nobs_monV<-dim(monV)[1]
colnames(monV)<-c("year","julian","Cstrata","habitat","catch")

# isolate data for age 0 fish only
tt<-subset(monC,monC$hybama01==0&monC$month>6&monC$month<10)
mon0<-cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama0)
colnames(mon0)<-c("year","julian","Cstrata","habitat","catch")
mon0_effQ<-cbind(tt$effort,tt$cQ/1000)
Nobs_mon0<-dim(mon0)[1]

#  isolate data where age -0  and  age - 1+ fish were not separated
tt<-subset(monC,monC$hybama01>0&monC$month>6&monC$month<10)
mon01<-cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama0+tt$hybama1+tt$hybama01)
mon01_effQ<-cbind(tt$effort,tt$cQ/1000)
#now do oct and november data -
tt<-subset(monC,monC$month==10)
mon01<-rbind(mon01,cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama0+tt$hybama1+tt$hybama01))
mon01_effQ<-rbind(mon01_effQ,cbind(tt$effort,tt$cQ/1000))
# add november data
mon01<-rbind(mon01,cbind((monre$year-2001),215,monre$Cstrata,monre$type,monre$hybama_1))
mon01_effQ<-rbind(mon01_effQ,cbind(monre$effort,monre$cQ/1000))
colnames(mon01)<-c("year","julian","Cstrata","habitat","catch")
Nobs_mon01<-dim(mon01)[1]

##read in results of population estimation drawn from Dudley et al., 2012
Nz<-c(1108430,1387948,267272,122381)
lNz<-log(Nz)
cvz<-c(0.3,0.259,0.372,0.376)
Cz<-exp(1.96*sqrt(log(1+cvz^2)))
lCz<-log(Cz)/1.96

###read in, subset and reformat fish rescue data using drying data
trescue<-subset(trescue,trescue[,1]>2008&trescue[,1]<2019)
trescue$jul<-NA
for (i in 1:length(trescue$jul)){
	trescue$jul[i]<-julian(as.Date(paste(trescue$year[i],trescue$month[i],trescue$dom[i],sep="-")),as.Date(paste(trescue$year[i],"03","31",sep="-")))[[1]]}
trescue$strata<-findInterval(trescue$rm,vgrid)
tR0<-matrix(NA,nrow=(Nyears-7),ncol=Nvstrata)
tR1<-matrix(NA,nrow=(Nyears-7),ncol=Nvstrata)
for (t in 1:(Nyears-7)){
	for (j in 1:Nvstrata){
		if(length(subset(trescue$jul,trescue$year==(t+2008)&trescue$strata==j))>0){
			tday1<-max(c(min(subset(trescue$jul,trescue$year==(t+2008)&trescue$strata==j))+14,fdryday[(t+7),j]+4),na.rm=T)
			if (tday1>0){
				temp<-subset(trescue,trescue$year==(t+2008)&trescue$strata==j&trescue$jul<tday1)
				tR1[t,j]<-sum(as.numeric(paste(temp$adult.alive)),na.rm=T)
				tR0[t,j]<-sum(temp$yoy.alive,na.rm=T)
	}}}}
#
Ntotjul<-215
StrataLen_rm<-array(NA,dim=c(Nyears,Ntotjul,3))
prop_nd<-array(NA,dim=c(Nyears,Ntotjul,3))
cum_nd<-array(NA,dim=c(Nyears,Ntotjul,3))
R0<-matrix(NA,ncol=4,nrow=(Nyears-7)*Ntotjul*Nstrata)
R1<-matrix(NA,ncol=4,nrow=(Nyears-7)*Ntotjul*Nstrata)
R0[,1]<-rep(c(8:Nyears),each=Ntotjul*Nstrata)
R1[,1]<-rep(c(8:Nyears),each=Ntotjul*Nstrata)
R0[,2]<-rep(c(1:Ntotjul),(Nyears-7)*Nstrata)
R1[,2]<-rep(c(1:Ntotjul),(Nyears-7)*Nstrata)
R0[,3]<-rep(rep(c(1:Nstrata),each=Ntotjul),(Nyears-7))
R1[,3]<-rep(rep(c(1:Nstrata),each=Ntotjul),(Nyears-7))
maxStrataLen<-table(crossgrid)/5
cum_phiR<-array(NA,dim=c(Nyears,Ntotjul,3))
# calculate weighted average survival of rescued fish
phiR$jul<-NA	
for (j in 1:12){phiR$jul[j]<-julian(as.Date(paste(2000,phiR$month[j],phiR$day[j],sep="-")),as.Date("2000-03-31"))[[1]]}
tm<-glm(cbind(phiR[,4],phiR[,3]-phiR[,4])~phiR[,5],family="binomial")
pred_phiR<-plogis(coef(tm)[1]+coef(tm)[2]*c(1:Ntotjul))
for (t in 1:Nyears){
	for (d in 1:Ntotjul){
		StrataLen_rm[t,d,1]<-length(which(crossgrid==1&(is.na(fdryday[t,])==T|fdryday[t,]>d|ldryday[t,]<d)))/5
		StrataLen_rm[t,d,2]<-length(which(crossgrid==2&(is.na(fdryday[t,])==T|fdryday[t,]>d|ldryday[t,]<d)))/5
		StrataLen_rm[t,d,3]<-length(which(crossgrid==3&(is.na(fdryday[t,])==T|fdryday[t,]>d|ldryday[t,]<d)))/5
		t1<-which(crossgrid==1&fdryday[t,]==d)
		t2<-which(crossgrid==2&fdryday[t,]==d)
		t3<-which(crossgrid==3&fdryday[t,]==d)
		prop_nd[t,d,1]<-length(t1)/5/maxStrataLen[1]
		prop_nd[t,d,2]<-length(t2)/5/maxStrataLen[2]
		prop_nd[t,d,3]<-length(t3)/5/maxStrataLen[3]
		if (t>7&length(t1)>0){
			R0[which(R0[,1]==t&R0[,2]==d&R0[,3]==1),4]<-sum(tR0[(t-7),t1])
			R1[which(R1[,1]==t&R1[,2]==d&R1[,3]==1),4]<-sum(tR1[(t-7),t1])
			}
		if (t>7&length(t2)>0){
			R0[which(R0[,1]==t&R0[,2]==d&R0[,3]==2),4]<-sum(tR0[(t-7),t2])
			R1[which(R1[,1]==t&R1[,2]==d&R1[,3]==2),4]<-sum(tR1[(t-7),t2])
			}
		if (t>7&length(t3)>0){
			R0[which(R0[,1]==t&R0[,2]==d&R0[,3]==3),4]<-sum(tR0[(t-7),t3])
			R1[which(R1[,1]==t&R1[,2]==d&R1[,3]==3),4]<-sum(tR1[(t-7),t3])
			}
		}
	cum_nd[t,,1]<-cumsum(prop_nd[t,,1])
	cum_nd[t,,2]<-cumsum(prop_nd[t,,2])
	cum_nd[t,,3]<-cumsum(prop_nd[t,,3])
	cum_phiR[t,,1]<-cumsum(prop_nd[t,,1]*pred_phiR)
	cum_phiR[t,,2]<-cumsum(prop_nd[t,,2]*pred_phiR)
	cum_phiR[t,,3]<-cumsum(prop_nd[t,,3]*pred_phiR)
	}
R0<-subset(R0,is.na(R0[,4])==FALSE&R0[,2]>91)
R1<-subset(R1,is.na(R1[,4])==FALSE)
NR0<-length(R0[,1])
NR1<-length(R1[,1])
# convert to rkm
StrataLen<-array(NA,dim=c(Nyears,Ntotjul,3))
StrataLen[,,1]<-92.1*StrataLen_rm[,,1]/maxStrataLen[1]
StrataLen[,,2]<-85.5*StrataLen_rm[,,2]/maxStrataLen[2]
StrataLen[,,3]<-65*StrataLen_rm[,,3]/maxStrataLen[3]

## reformat summarize expert ellicitation data
ee<-list(ee1=ee1,ee2=ee2,ee3=ee3,ee4=ee4,ee5=ee5)
# function to calculate standard error from quantile, upper, lower and mean info
calcsig<-function(up,down,mean,q,INT=c(0,1000)){
	sig<-function(x){abs(pnorm(up,mean,x)-pnorm(down,mean,x)-q)}
	optimize(sig,interval=INT)$minimum}
## 5a - create prior on movement out of reach
#extract and convert info to sd
emove<-matrix(NA,nrow=5,ncol=2)
for (i in 1:5){
	temp<-ee[[i]][ee[[i]][,1]=="5a",]
	emove[i,1]<-temp[1,6]
	emove[i,2]<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
	}
## 6a,6b,6c - create priors on river width
#extract and convert info to sd
Nexperts<-5
ewidths<-array(NA,dim=c(5,3,2,8))
for (i in 1:5){
	temp<-ee[[i]][ee[[i]][,1]=="6a - river width",]
	temp2<-ee[[i]][ee[[i]][,1]=="6b - river width",]
	temp3<-ee[[i]][ee[[i]][,1]=="6c - river width",]
	for (j in 1:8){
		ewidths[i,1,1,j]<-temp[j,6]
		ewidths[i,2,1,j]<-temp2[j,6]
		ewidths[i,3,1,j]<-temp3[j,6]
		ewidths[i,1,2,j]<-calcsig(temp[j,4],temp[j,3],temp[j,6],temp[j,5]/100)
		ewidths[i,2,2,j]<-calcsig(temp2[j,4],temp2[j,3],temp2[j,6],temp2[j,5]/100)
		ewidths[i,3,2,j]<-calcsig(temp3[j,4],temp3[j,3],temp3[j,6],temp3[j,5]/100)
		}}
refQ<-c(5,50,100,150,200,250,500,1000)/1000

## fucntions to calculate the expert informed flow covariate
lin_int<-function(x,xlow,xhi,ylow,yhi){(yhi-ylow)*(x-xlow)/(xhi-xlow)+ylow}
q2hab<-function(q,pars){
	p<-c(0,pars,pars[length(pars)])
	qs<-c(0,5,50,100,150,200,250,1000,1500,2000,2500,3000,4000,5000,6000,7000,2*10^4)
	t1<-findInterval(q,qs)
	lin_int(q,qs[t1],qs[(t1+1)],p[t1],p[(t1+1)])
	}
prop<-function(pars){
	t<-seq(2,132,10)
	tout<-pars[1]
	for (i in 2:131){
		t1<-findInterval(i,t)
		tout[i]<-lin_int(i,t[t1],t[(t1+1)],pars[t1],pars[(t1+1)])}
	tout[132]<-pars[14]
	tout/sum(tout)}
calc_cov<-function(Q,q2fPARS,kappa,D,prPARS){
	thab<-numeric()
	for (i in 1:214){thab[i]<-q2hab(Q[i],pars=q2fPARS)}
	tprop<-prop(prPARS)
	t2hab<-numeric()
	for (i in 1:132){t2hab[i]<-min(thab[i:(i+D)])}
	tstart<-min(c(which(Q>kappa),which(Q[-1]-Q[-length(Q)]>100),132))
	out<-sum(tprop[1:tstart])*t2hab[tstart]+sum(tprop[tstart:132]*t2hab[tstart:132])
	return(out)}

# calculate predictions by different experts under different designs	
preds<-array(NA,dim=c(17,3,9,5))
riversegment<-c("3e","3c","3a")
design<-matrix(0,nrow=9,ncol=4)
design[1,1]<-1
design[2,1]<-(-1)
design[3,2]<-1
design[4,2]<-(-1)
design[5,3]<-1
design[6,3]<-(-1)
design[7,4]<-1
design[8,4]<-(-1)

for (i in 1:5){
	for (d in 1:9){
		temp<-subset(ee[[i]],ee[[i]][,1]=="1"&is.na(ee[[i]][,5])==FALSE)
		t1<-numeric()
		for (t in 1:6){
			se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
			t1[t]<-temp[t,6]+se*design[d,1]
			se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
			t1[(t+7)]<-temp[(t+6),6]-se*design[d,1]
			}
		t1[7]<-1
		t1[14]<-0
		temp<-subset(ee[[i]],ee[[i]][,1]=="2")
		se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
		t2<-temp[1,6]+se*design[d,2]
		temp<-subset(ee[[i]],ee[[i]][,1]=="4")
		se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
		t4<-temp[1,6]+se*design[d,4]
		for (r in 1:3){
			temp<-subset(ee[[i]],ee[[i]][,1]==riversegment[r]&is.na(ee[[i]][,5])==FALSE)
			t3<-numeric()
			for (t in 1:6){
				se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
				t3[t]<-temp[t,6]+se*design[d,3]
				}	
			for (t in 8:15){
				se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
				t3[t]<-temp[t,6]-se*design[d,3]
				}
			t3[7]<-temp[7,6]
			for (j in 1:17){		
				if (r==3) {q<-aQ[,j]} else {q<-sQ[,j]}
				preds[j,r,d,i]<-calc_cov(q,t3,t2,t4,t1)
				}}}}


######################
#### STEP 3: Fit first round of models - running models takes many hours
######################	

			
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
ct<-function(x){x-mean(x)}
# first expert
# design 1
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,1,1]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_1_1 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 2
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,2,1]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_1_2 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 3
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,3,1]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_1_3 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 4
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,4,1]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_1_4 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 5
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,5,1]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_1_5 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 6
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,6,1]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_1_6 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 7
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,7,1]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_1_7 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 8
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,8,1]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_1_8 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 9
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,9,1]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_1_9 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# second expert
# design 1
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,1,2]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_2_1 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 2
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,2,2]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_2_2 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 3
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,3,2]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_2_3 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 4
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,4,2]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_2_4 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 5
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,5,2]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_2_5 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 6
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,6,2]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_2_6 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 7
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,7,2]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_2_7 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 8
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,8,2]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_2_8 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 9
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,9,2]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_2_9 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# third expert
# design 1
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,1,3]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_3_1 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 2
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,2,3]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_3_2 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 3
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,3,3]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_3_3 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 4
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,4,3]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_3_4 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 5
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,5,3]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_3_5 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 6
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,6,3]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_3_6 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 7
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,7,3]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_3_7 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 8
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,8,3]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_3_8 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 9
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,9,3]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_3_9 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# fourth expert
# design 1
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,1,4]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_4_1 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 2
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,2,4]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_4_2 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 3
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,3,4]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_4_3 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 4
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,4,4]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_4_4 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 5
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,5,4]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_4_5 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 6
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,6,4]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_4_6 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 7
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,7,4]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_4_7 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 8
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,8,4]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_4_8 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 9
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,9,4]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_4_9 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# fifth expert
# design 1
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,1,5]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_5_1 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 2
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,2,5]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_5_2 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 3
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,3,5]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_5_3 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 4
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,4,5]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_5_4 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 5
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,5,5]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_5_5 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 6
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,6,5]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_5_6 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 7
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,7,5]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_5_7 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 8
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,8,5]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_5_8 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
# design 9
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds[,,9,5]),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M_5_9 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 1000) 
#check convergence and rerun any models were max rhat is greater than 1.1
c(max(summary(M_1_1)$summary[,10],na.rm=T),max(summary(M_1_2)$summary[,10],na.rm=T),max(summary(M_1_3)$summary[,10],na.rm=T),max(summary(M_1_4)$summary[,10],na.rm=T),max(summary(M_1_5)$summary[,10],na.rm=T),max(summary(M_1_6)$summary[,10],na.rm=T),max(summary(M_1_7)$summary[,10],na.rm=T),max(summary(M_1_8)$summary[,10],na.rm=T),max(summary(M_1_9)$summary[,10],na.rm=T),
max(summary(M_2_1)$summary[,10],na.rm=T),max(summary(M_2_2)$summary[,10],na.rm=T),max(summary(M_2_3)$summary[,10],na.rm=T),max(summary(M_2_4)$summary[,10],na.rm=T),max(summary(M_2_5)$summary[,10],na.rm=T),max(summary(M_2_6)$summary[,10],na.rm=T),max(summary(M_2_7)$summary[,10],na.rm=T),max(summary(M_2_8)$summary[,10],na.rm=T),max(summary(M_2_9)$summary[,10],na.rm=T),
max(summary(M_3_1)$summary[,10],na.rm=T),max(summary(M_3_2)$summary[,10],na.rm=T),max(summary(M_3_3)$summary[,10],na.rm=T),max(summary(M_3_4)$summary[,10],na.rm=T),max(summary(M_3_5)$summary[,10],na.rm=T),max(summary(M_3_6)$summary[,10],na.rm=T),max(summary(M_3_7)$summary[,10],na.rm=T),max(summary(M_3_8)$summary[,10],na.rm=T),max(summary(M_3_9)$summary[,10],na.rm=T),
max(summary(M_4_1)$summary[,10],na.rm=T),max(summary(M_4_2)$summary[,10],na.rm=T),max(summary(M_4_3)$summary[,10],na.rm=T),max(summary(M_4_4)$summary[,10],na.rm=T),max(summary(M_4_5)$summary[,10],na.rm=T),max(summary(M_4_6)$summary[,10],na.rm=T),max(summary(M_4_7)$summary[,10],na.rm=T),max(summary(M_4_8)$summary[,10],na.rm=T),max(summary(M_4_9)$summary[,10],na.rm=T),
max(summary(M_5_1)$summary[,10],na.rm=T),max(summary(M_5_2)$summary[,10],na.rm=T),max(summary(M_5_3)$summary[,10],na.rm=T),max(summary(M_5_4)$summary[,10],na.rm=T),max(summary(M_5_5)$summary[,10],na.rm=T),max(summary(M_5_6)$summary[,10],na.rm=T),max(summary(M_5_7)$summary[,10],na.rm=T),max(summary(M_5_8)$summary[,10],na.rm=T),max(summary(M_5_9)$summary[,10],na.rm=T))
# then calculate R2's to be reported in appendix S3
calc_R2<-function(mod){
	eps1<-apply(extract(mod,"lbeta_eps")[[1]][,1,],1,var)
	lpred1<-apply(log(extract(mod,"Rf")[[1]][,,1]),1,var)
	eps2<-apply(extract(mod,"lbeta_eps")[[1]][,2,],1,var)
	lpred2<-apply(log(extract(mod,"Rf")[[1]][,,2]),1,var)
	eps3<-apply(extract(mod,"lbeta_eps")[[1]][,3,],1,var)
	lpred3<-apply(log(extract(mod,"Rf")[[1]][,,3]),1,var)
	r2_1<-1-eps1/lpred1
	r2_2<-1-eps2/lpred2
	r2_3<-1-eps3/lpred3
	return(c(mean(r2_1),mean(r2_2),mean(r2_3)))}
#
r2s<-round(rbind(calc_R2(M_1_1),calc_R2(M_1_2),calc_R2(M_1_3),calc_R2(M_1_4),calc_R2(M_1_5),calc_R2(M_1_6),calc_R2(M_1_7),calc_R2(M_1_8),calc_R2(M_1_9),
calc_R2(M_2_1),calc_R2(M_2_2),calc_R2(M_2_3),calc_R2(M_2_4),calc_R2(M_2_5),calc_R2(M_2_6),calc_R2(M_2_7),calc_R2(M_2_8),calc_R2(M_2_9),
calc_R2(M_3_1),calc_R2(M_3_2),calc_R2(M_3_3),calc_R2(M_3_4),calc_R2(M_3_5),calc_R2(M_3_6),calc_R2(M_3_7),calc_R2(M_3_8),calc_R2(M_3_9),
calc_R2(M_4_1),calc_R2(M_4_2),calc_R2(M_4_3),calc_R2(M_4_4),calc_R2(M_4_5),calc_R2(M_4_6),calc_R2(M_4_7),calc_R2(M_4_8),calc_R2(M_4_9),
calc_R2(M_5_1),calc_R2(M_5_2),calc_R2(M_5_3),calc_R2(M_5_4),calc_R2(M_5_5),calc_R2(M_5_6),calc_R2(M_5_7),calc_R2(M_5_8),calc_R2(M_5_9)),3)


######################
#### STEP 4: Fit second round of models based on results of first round- running models takes many hours
######################	

# calculate new set of larval carry capacity indices based on expert 3 (best performing in round 1)	
preds2<-array(NA,dim=c(17,3))
temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
t1<-numeric()
for (t in 1:6){
	se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
	t1[t]<-temp[t,6]-se
	se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
	t1[(t+7)]<-temp[(t+6),6]+se
	}
t1[7]<-1
t1[14]<-0
temp<-subset(ee[[3]],ee[[3]][,1]=="2")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t2<-temp[1,6]
temp<-subset(ee[[3]],ee[[3]][,1]=="4")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t4<-temp[1,6]-se
#
for (r in 1:3){
	temp<-subset(ee[[3]],ee[[3]][,1]==riversegment[r]&is.na(ee[[3]][,5])==FALSE)
	t3<-numeric()
	for (t in 1:6){
		se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
		t3[t]<-temp[t,6]+se
		}	
	for (t in 8:15){
		se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
		t3[t]<-temp[t,6]-se
		}
	t3[7]<-temp[7,6]
	for (j in 1:17){
		if (r==3) {q<-aQ[,j]} else {q<-sQ[,j]}
		preds2[j,r]<-calc_cov(q,t3,t2,t4,t1)
		}}
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(preds2),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M2_1 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 2000) 
#allow for interannual variability in survival
M2_1re <- stan("int3re.stan",data = rgsm.data,chains = 3, iter = 10000) 

## compare to model based on average flow in may-june
mjflow_sana<-subset(sanaQ,(sanaQ$month==5|sanaQ$month==6)&sanaQ$year>2001&sanaQ$year<2019)
mjflow_ang<-subset(angQ,(angQ$month==5|angQ$month==6)&angQ$year>2001&angQ$year<2019)
simpflow<-cbind(tapply(mjflow_sana$cfs,mjflow_sana$year,mean),tapply(mjflow_sana$cfs,mjflow_sana$year,mean),tapply(mjflow_ang$cfs,mjflow_ang$year,mean))/1000
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
   ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(simpflow),
   mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
   w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)
M2_2 <- stan("int3.stan",data = rgsm.data,chains = 3, iter = 2000) 
M2_2re <- stan("int3re.stan",data = rgsm.data,chains = 3, iter = 10000) 
###


######################
#### STEP 5: Do model checks and out of sample comparisons and make figures 2 - 6 and S1 - S4 from best model
######################	


# calculations for fig S1
lp0<-extract(M2_1re,"lpmon0")[[1]]
lp1<-extract(M2_1re,"lpmon1")[[1]]
lp01<-extract(M2_1re,"lpmon01")[[1]]
lpV<-extract(M2_1re,"lpmonV")[[1]]
lR0<-extract(M2_1re,"lpR0")[[1]]
lR1<-extract(M2_1re,"lpR1")[[1]]
sz<-extract(M2_1re,"sz")[[1]]
rsz<-extract(M2_1re,"rsz")[[1]]
#
pc0<-lp0
pc1<-lp1
pc01<-lp01
pcV<-lpV
pr0<-lR0
pr1<-lR1
#
q_pc0<-numeric()
q_pc1<-numeric()
q_pc01<-numeric()
q_pcV<-numeric()
q_pr0<-numeric()
q_pr1<-numeric()
qf<-function(obs,pred){
	temp<-which(obs==sort(pred))
	if (length(temp)>0) {sample(temp,1)/length(pred)} else {
		temp2<-findInterval(obs,sort(pred))
			if (temp2==0) {1/length(pred)} else { temp2/length(pred)}}}
#
for (i in 1:(dim(lp0)[2])){pc0[,i]<-rnbinom(15000,mu=exp(lp0[,i]),size=sz)
	q_pc0[i]<-qf(mon0[i,5],pc0[,i])}
for (i in 1:(dim(lp1)[2])){pc1[,i]<-rnbinom(15000,mu=exp(lp1[,i]),size=sz)
	q_pc1[i]<-qf(mon1[i,5],pc1[,i])}
for (i in 1:(dim(lp01)[2])){pc01[,i]<-rnbinom(15000,mu=exp(lp01[,i]),size=sz)
	q_pc01[i]<-qf(mon01[i,5],pc01[,i])}
for (i in 1:(dim(lpV)[2])){pcV[,i]<-rnbinom(15000,mu=exp(lpV[,i]),size=sz)
	q_pcV[i]<-qf(monV[i,5],pcV[,i])}
for (i in 1:(dim(lR0)[2])){pr0[,i]<-rnbinom(15000,mu=exp(lR0[,i]),size=rsz)
	q_pr0[i]<-qf(R0[i,4],pr0[,i])}
for (i in 1:(dim(lR1)[2])){pr1[,i]<-rnbinom(15000,mu=exp(lR1[,i]),size=rsz)
	q_pr1[i]<-qf(R1[i,4],pr1[,i])}
#combine monitoring data for plotting
q_mon<-c(q_pc0,q_pc1,q_pc01,q_pcV)
pc_mon<-cbind(pc0,pc1,pc01,pcV)
obsc_mon<-c(mon0[,5],mon1[,5],mon01[,5],monV[,5])
Q_mon<-c(mon0_effQ[,2],mon1_effQ[,2],mon01_effQ[,2],monV_effQ[,2])
#ditto for rescue
q_res<-c(q_pr0,q_pr1)
pc_res<-cbind(pr0,pr1)
obsc_res<-c(R0[,4],R1[,4])
jul_res<-c(R0[,2],R1[,2])
# code to make actual Fig S1 plot
par(mfrow=c(4,2))
par(mar=c(4,4,1,1))
xymax<-max(c(obsc_mon,apply(pc_mon,2,mean)))
plot(1+obsc_mon,1+apply(pc_mon,2,mean),ylim=c(1,1+xymax),xlim=c(1,1+xymax),axes=FALSE,xlab="",ylab="",pch=19,col=rgb(0,0,0,.1),main="Monitoring data",log="xy")
axis(1,at=c(1,2,11,101,1001,10001),labels=c(0,1,10,100,1000,10000),pos=1)
axis(2,at=c(1,2,11,101,1001,10001),labels=c(0,1,10,100,1000,10000),pos=1,las=T)
mtext("Observed catch",1,2,cex=0.7)
mtext("Predicted catch",2,2,cex=0.7)
text(1000,3,"R2 = 0.76")
text(.3,3000,"A)",xpd=T)
#
xymax<-max(c(obsc_res,apply(pc_res,2,mean)))
plot(1+obsc_res,1+apply(pc_res,2,mean),ylim=c(1,1+xymax),xlim=c(1,1+xymax),axes=FALSE,xlab="",ylab="",pch=19,col=rgb(0,0,0,.1),main="Rescue data",log="xy")
axis(1,at=c(1,2,11,101,1001,10001),labels=c(0,1,10,100,1000,10000),pos=1)
axis(2,at=c(1,2,11,101,1001,10001),labels=c(0,1,10,100,1000,10000),pos=1,las=T)
mtext("Observed catch",1,2,cex=0.7)
mtext("Predicted catch",2,2,cex=0.7)
text(8000,3,"R2 = 0.77")
text(.15,30000,"B)",xpd=T)
# % zeros
iszero<-function(x){ifelse(x==0,1,0)}
hist(apply(iszero(pc_mon),1,mean),axes=FALSE,xlab="",ylab="",pch=19,main="")
abline(v=mean(iszero(obsc_mon)),col="red",lty=2,lwd=2)
axis(1,pos=0)
axis(2,pos=min(apply(iszero(pc_mon),1,mean)),at=c(0,750,1500,2250),labels=seq(0,0.15,.05),las=T)
mtext("Proportion of catch equal to zero",1,2,cex=0.7)
mtext("% of simulations",2,2,cex=0.7)
text(mean(iszero(obsc_mon))+.01,1500,"Observed",col="red")
text(0.52,2600,"C)",xpd=T)
# % zeros
hist(apply(iszero(pc_res),1,mean),axes=FALSE,xlab="",ylab="",pch=19,main="")
abline(v=mean(iszero(obsc_res)),col="red",lty=2,lwd=2)
axis(1,pos=0)
axis(2,pos=min(apply(iszero(pc_res),1,mean)),at=c(0,1500,3000,4500),labels=seq(0,0.3,.1),las=T)
mtext("Proportion of catch equal to zero",1,2,cex=0.7)
mtext("% of simulations",2,2,cex=0.7)
text(mean(iszero(obsc_res))-0.03,4000,"Observed",col="red")
text(0.14,4500,"D)",xpd=T)
#qq
plot(sort(q_mon),c(1:length(q_mon))/length(q_mon),pch=19,col=rgb(0,0,0,.01),main="",xlab="",ylab="",axes=FALSE)
curve(1*x,add=T,lty=2,lwd=2)
axis(1,pos=0)
axis(2,pos=0,las=T)
mtext("Observed quantile",1,2,cex=0.7)
mtext("Expected quantile",2,2,cex=0.7)
text(-0.16,1.1,"E)",xpd=T)
#qq
plot(sort(q_res),c(1:length(q_res))/length(q_res),pch=19,col=rgb(0,0,0,.1),main="",xlab="",ylab="",axes=FALSE)
curve(1*x,add=T,lty=2,lwd=2)
axis(1,pos=0)
axis(2,pos=0,las=T)
mtext("Observed quantile",1,2,cex=0.7)
mtext("Expected quantile",2,2,cex=0.7)
text(-0.16,1.1,"F)",xpd=T)
#
plot(Q_mon,q_mon,pch=19,col=rgb(0,0,0,.1),main="",xlab="",ylab="",axes=FALSE)
axis(1,at=c(0,0.5,1),labels=c(0,500,1000),pos=0)
axis(2,pos=0,las=T)
mtext("Discharge (cfs)",1,2,cex=0.7)
mtext("Observed quantile",2,2,cex=0.7)
text(-0.16,1.1,"G)",xpd=T)
#
plot(jul_res,q_res,pch=19,col=rgb(0,0,0,.1),main="",xlab="",ylab="",axes=FALSE)
axis(1,pos=0)
axis(2,pos=0,las=T)
mtext("Days after April 1st",1,2,cex=0.7)
mtext("Observed quantile",2,2,cex=0.7)
text(-27,1.1,"H)",xpd=T)

# checking priors vs. posteriors
# code for figs S2 and S3 - can do for any fitted model with time variation in mortality rates - set mod = M2_1re for actual plots in appendix
checkpp_rep<-function(mod){
  par(mfrow=c(1,3))
  plot(density(extract(mod,"a")[[1]]),xlim=c(300,1500),axes=FALSE,main="",xlab="a",col="red")
  curve(dnorm(x,675,135),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=300,las=T)
  plot(density(extract(mod,"beta_stk")[[1]]),xlim=c(0.5,3),axes=FALSE,main="",xlab="beta_stk",col="red")
  curve(dnorm(x,1,0.1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0.5,las=T)
  plot(density(extract(mod,"beta_2")[[1]]),xlim=c(0.5,3),axes=FALSE,main="",xlab="beta_2",col="red")
  curve(dnorm(x,2,0.1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0.5,las=T)
 } 

checkpp_sd<-function(mod){
  par(mfrow=c(4,4))
  for (k in 1:3){
     plot(density(extract(mod,"mu_lM0")[[1]][,k]),xlim=c(-8,-3),axes=FALSE,main="",xlab=paste("mu_lM0",k),col="red")
	curve(dunif(x,-8,-3),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=-8,las=T)
	}
  for (k in 1:3){
     plot(density(extract(mod,"mu_lM1")[[1]][,k]),xlim=c(-8,-3),axes=FALSE,main="",xlab=paste("mu_lM1",k),col="red")
	curve(dunif(x,-8,-3),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=-8,las=T)
	}
  for (k in 1:3){
     plot(density(extract(mod,"mu_lMw")[[1]][,k]),xlim=c(-8,-3),axes=FALSE,main="",xlab=paste("mu_lMw",k),col="red")
	curve(dunif(x,-8,-3),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=-8,las=T)
	}
    plot(density(extract(mod,"irphi")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="irphi",col="red")
	curve(dunif(x,0,1),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=0,las=T)
	plot(density(extract(mod,"p0")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="p0",col="red")
	curve(dunif(x,0,1),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=0,las=T)
	   plot(density(extract(mod,"p1")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="p1",col="red")
	curve(dunif(x,0,1),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=0,las=T)
   plot(density(extract(mod,"rp0")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="rp0",col="red")
	curve(dunif(x,0,1),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=0,las=T)
   plot(density(extract(mod,"rp1")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="rp1",col="red")
	curve(dunif(x,0,1),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=0,las=T)
  plot(density(extract(mod,"sz")[[1]]),main="",xlab="sz",col="red",axes=FALSE,xlim=c(0.25,0.6))
  	curve(dunif(x,0,1),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=0.25,las=T)
   plot(density(extract(mod,"rsz")[[1]]),main="",xlab="rsz",col="red",axes=FALSE,xlim=c(0.25,0.6))
   	curve(dunif(x,0,1),col="blue",add=T,lty=2)
	axis(1,pos=0)
	axis(2,pos=0.25,las=T)
	}

###### calculations for figure 2
t1<-matrix(NA,ncol=14,nrow=6)
for (i in 1:5){
	temp<-subset(ee[[i]],ee[[i]][,1]=="1"&is.na(ee[[i]][,5])==FALSE)
		for (t in 1:6){
			t1[i,t]<-temp[t,6]
			t1[i,(t+7)]<-temp[(t+6),6]
			}
		t1[i,7]<-1
		t1[i,14]<-0
	}
temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
for (t in 1:6){
	se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
	t1[6,t]<-temp[t,6]-se
	se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
	t1[6,(t+7)]<-temp[(t+6),6]+se
	}
t1[6,7]<-1
t1[6,14]<-0
#
t3<-matrix(NA,ncol=15,nrow=6)
for (i in 1:5){
	temp<-subset(ee[[i]],ee[[i]][,1]==riversegment[3]&is.na(ee[[i]][,5])==FALSE)
	for (t in 1:15){t3[i,t]<-temp[t,6]}
	}
temp<-subset(ee[[3]],ee[[3]][,1]==riversegment[3]&is.na(ee[[3]][,5])==FALSE)
for (t in 1:6){
	se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
	t3[6,t]<-temp[t,6]+se
	}	
for (t in 8:15){
	se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
	t3[6,t]<-temp[t,6]-se
	}
t3[6,7]<-temp[7,6]
#
t2<-numeric()
t2se<-numeric()
t4<-numeric()
t4se<-numeric()
for (i in 1:5){
	temp<-subset(ee[[i]],ee[[i]][,1]=="2")
	t2[i]<-temp[1,6]
	t2se[i]<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
	temp<-subset(ee[[i]],ee[[i]][,1]=="4")
	t4[i]<-temp[1,6]
	t4se[i]<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
	}
# code to actually plot Fig 2
layout(matrix(c(1,1,1,2,3,3,3,4,4,5,5,6,6,6),byrow=T,ncol=7))
par(mar=c(4,5,1,1))
plot(prop(t1[6,]),axes=FALSE,ylab="",xlab="",type="l",lwd=3,ylim=c(0,.025),col=rgb(0,0,0,.6))
points(prop(t1[1,]),type="l",col=rgb(1,0,0,.4),lwd=2)
points(prop(t1[2,]),type="l",col=rgb(0,1,0,.4),lwd=2)
points(prop(t1[3,]),type="l",col=rgb(0,0,1,.4),lwd=2)
points(prop(t1[4,]),type="l",col=rgb(0,1,1,.4),lwd=2)
points(prop(t1[5,]),type="l",col=rgb(1,0,1,.4),lwd=2)
axis(2,las=T,pos=1)
mtext("Proportion of eggs laid",2,3,cex=0.7)
axis(1,c(1,32,62,93,123),c("Mar-1","Apr-1","May-1","Jun-1","Jul-1"),pos=0)
axis(1,c(1,132),c("",""),tck=FALSE,pos=0)
mtext("Month-Day",1,2,cex=0.7)
text(-33,.02625,"A)",xpd=T)
#legend
par(mar=c(1,1,1,1))
plot(1,axes=FALSE,ylab="",xlab="",ylim=c(0,1),xlim=c(0,1),col="white")
legend(0,1,title="Expert",c("1 - mean","2 - mean","3 - mean","4 - mean","5 - mean", "3 - best","no expert"),pch=c(rep(19,6),8),
	col=c(rgb(1,0,0,.4),rgb(0,1,0,.4),rgb(0,0,1,.4),rgb(0,1,1,.4),rgb(1,0,1,.4),rgb(0,0,0,.6),rgb(0,0,0,.6)),bty="n",border=NA,y.intersp=1.5)
#
par(mar=c(4,5,1,1))
plot(c(5:7000),q2hab(c(5:7000),t3[6,])/min(q2hab(c(5:7000),t3[6,])),axes=FALSE,lwd=3,ylab="",xlab="",type="l",col=rgb(0,0,0,.6),log="xy",ylim=c(1,1000))
points(c(5:7000),q2hab(c(5:7000),t3[1,])/min(q2hab(c(5:7000),t3[1,])),type="l",col=rgb(1,0,0,.4),lwd=2)
points(c(5:7000),q2hab(c(5:7000),t3[2,])/min(q2hab(c(5:7000),t3[2,])),type="l",col=rgb(0,1,0,.4),lwd=2)
points(c(5:7000),q2hab(c(5:7000),t3[3,])/min(q2hab(c(5:7000),t3[3,])),type="l",col=rgb(0,0,1,.4),lwd=2)
points(c(5:7000),q2hab(c(5:7000),t3[4,])/min(q2hab(c(5:7000),t3[4,])),type="l",col=rgb(0,1,1,.4),lwd=2)	
points(c(5:7000),q2hab(c(5:7000),t3[5,])/min(q2hab(c(5:7000),t3[5,])),type="l",col=rgb(1,0,1,.4),lwd=2)
mtext("Relative Amount of Larval habitat",2,2.5,cex=0.7)
mtext("Discharge",1,2,cex=0.7)
axis(2,at=c(1,10,100,1000),labels=c("1","10","10^2","10^3"),las=T,pos=5)
axis(1,c(5,50,500,5000),pos=1)
axis(1,c(0,7000),c("",""),tck=FALSE,pos=1)
text(.95,1400,"B)",xpd=T)
#
plot(t2,pch=19,col=c(rgb(1,0,0,.4),rgb(0,1,0,.4),rgb(0,0,1,.4),rgb(0,1,1,.4),rgb(1,0,1,.4)),xlab="",ylab="",ylim=c(0,1000),xlim=c(0.5,5.5),axes=FALSE)
segments(c(1:6),t2+t2se,c(1:6),t2-t2se,col=c(rgb(1,0,0,.4),rgb(0,1,0,.4),rgb(0,0,1,.4),rgb(0,1,1,.4),rgb(1,0,1,.4)))
axis(2,pos=0.5,las=T)
axis(1,at=c(0.5,1:5,7),labels=c("",1:5,""),pos=0)
points(3,t2[3],col=rgb(0,0,0,.6),pch=19)
mtext("Discharge that cues spawning (cfs)",2,3,cex=0.7)
mtext("Expert",1,2,cex=0.7)
text(-1.6,1075,"C)",xpd=T)
#
plot(t4,pch=19,col=c(rgb(1,0,0,.4),rgb(0,1,0,.4),rgb(0,0,1,.4),rgb(0,1,1,.4),rgb(1,0,1,.4)),xlab="",ylab="",ylim=c(0,35),xlim=c(0.5,5.5),axes=FALSE)
segments(c(1:6),t4+t4se,c(1:6),t4-t4se,col=c(rgb(1,0,0,.4),rgb(0,1,0,.4),rgb(0,0,1,.4),rgb(0,1,1,.4),rgb(1,0,1,.4)))
axis(2,pos=0.5,las=T)
axis(1,at=c(0.5,1:5,7),labels=c("",1:5,""),pos=0)
points(3,t4[3]-t4se[3],col=rgb(0,0,0,.6),pch=19)
mtext("Duration (days)",2,2.5,cex=0.7)
mtext("Expert",1,2,cex=0.7)
text(-1.2,37,"D)",xpd=T)
#
plot(c(seq(0.7,1.3,0.1),seq(1.7,2.3,0.1),seq(2.7,3.3,0.1)),c(r2s[seq(9,45,9),1],calc_R2(M2_1re)[1],calc_R2(M2_2re)[1],r2s[seq(9,45,9),2],calc_R2(M2_1re)[2],
	calc_R2(M2_2re)[2],r2s[seq(9,45,9),3],calc_R2(M2_1re)[3],calc_R2(M2_2re)[3]),
	col=rep(c(rgb(1,0,0,.4),rgb(0,1,0,.4),rgb(0,0,1,.4),rgb(0,1,1,.4),rgb(1,0,1,.4),col=rgb(0,0,0,.6),col=rgb(0,0,0,.6)),3),
	pch=c(rep(19,6),8),axes=FALSE,ylab="",xlab="",ylim=c(0,1),xlim=c(0.6,3.5))
axis(2,pos=0.6,las=T)
axis(1,at=c(0,1:3,4),labels=c("","San Acacia","Isleta","Angostura",""),pos=0)
mtext("Variance explained (R^2)",2,2.25,cex=0.7)
mtext("River Segment",1,2,cex=0.7)
text(0,1.06,"E)",xpd=T)

#### Out of sample predictions
### format out of sample data
oosD$ustrata<-findInterval(oosD$URM,vgrid)
oosD$dstrata<-findInterval(oosD$LRM,vgrid)
oosD$jul<-NA
for (i in 1:length(oosD$jul)){
	oosD$jul[i]<-julian(as.Date(paste(oosD$Year[i],oosD$month[i],oosD$dom[i],sep="-")),as.Date(paste(oosD$Year[i],"03","31",sep="-")))[[1]]}
fdryday_oos<-matrix(NA,nrow=2,ncol=(length(vgrid)-1))
ldryday_oos<-matrix(NA,nrow=2,ncol=(length(vgrid)-1))
for (k in 1:(length(vgrid)-1)){
	for (j in 1:2){
		temp2<-subset(oosD$jul,k<=oosD$ustrata&k>=oosD$dstrata&oosD$Year==(2018+j))
		fdryday_oos[j,k]<-ifelse(length(temp2)==0,NA,min(temp2))
		ldryday_oos[j,k]<-ifelse(length(temp2)==0,NA,max(temp2))
	}}
StrataLen_rm_oos<-array(NA,dim=c(2,Ntotjul,3))
prop_nd_oos<-array(NA,dim=c(2,Ntotjul,3))
cum_nd_oos<-array(NA,dim=c(2,Ntotjul,3))
cum_phiR_oos<-array(NA,dim=c(2,Ntotjul,3))
for (t in 1:2){
	for (d in 1:Ntotjul){
		StrataLen_rm_oos[t,d,1]<-length(which(crossgrid==1&(is.na(fdryday_oos[t,])==T|fdryday_oos[t,]>d|ldryday_oos[t,]<d)))/5
		StrataLen_rm_oos[t,d,2]<-length(which(crossgrid==2&(is.na(fdryday_oos[t,])==T|fdryday_oos[t,]>d|ldryday_oos[t,]<d)))/5
		StrataLen_rm_oos[t,d,3]<-length(which(crossgrid==3&(is.na(fdryday_oos[t,])==T|fdryday_oos[t,]>d|ldryday_oos[t,]<d)))/5
		t1<-which(crossgrid==1&fdryday_oos[t,]==d)
		t2<-which(crossgrid==2&fdryday_oos[t,]==d)
		t3<-which(crossgrid==3&fdryday_oos[t,]==d)
		prop_nd_oos[t,d,1]<-length(t1)/5/maxStrataLen[1]
		prop_nd_oos[t,d,2]<-length(t2)/5/maxStrataLen[2]
		prop_nd_oos[t,d,3]<-length(t3)/5/maxStrataLen[3]
		}
	cum_nd_oos[t,,1]<-cumsum(prop_nd_oos[t,,1])
	cum_nd_oos[t,,2]<-cumsum(prop_nd_oos[t,,2])
	cum_nd_oos[t,,3]<-cumsum(prop_nd_oos[t,,3])
	cum_phiR_oos[t,,1]<-cumsum(prop_nd_oos[t,,1]*pred_phiR)
	cum_phiR_oos[t,,2]<-cumsum(prop_nd_oos[t,,2]*pred_phiR)
	cum_phiR_oos[t,,3]<-cumsum(prop_nd_oos[t,,3]*pred_phiR)
	}
StrataLen_oos<-array(NA,dim=c(2,Ntotjul,3))
StrataLen_oos[,,1]<-92.1*StrataLen_rm_oos[,,1]/maxStrataLen[1]
StrataLen_oos[,,2]<-85.5*StrataLen_rm_oos[,,2]/maxStrataLen[2]
StrataLen_oos[,,3]<-65*StrataLen_rm_oos[,,3]/maxStrataLen[3]
#
vkeep2<-match(c("Date","Month","Year","C","L","S","Number","RM"),names(oosR))
oosR<-oosR[,vkeep2]
oosR$strata<-findInterval(converter[match(oosR$RM,converter[,1]),2],vgrid)
oosR$period<-ifelse(oosR$Month<4,8*(oosR$Year-2002)+1,ifelse(oosR$Month>9,8*(oosR$Year-2002)+9,8*(oosR$Year-2002)+oosR$Month-3))
oosR$winmon<-ifelse(oosR$Month<4,4-oosR$Month,ifelse(oosR$Month>9,16-oosR$Month,0))
#
tw_m_oos<-array(0,dim=c(2,Nvstrata,2)) #summarized at fine grid
w_m_oos<-array(0,dim=c(2,Nstrata,2)) #summarized at coarse grid
#
for (i in 1:2){
	for (j in 1:Nvstrata){
		temp<-subset(oosR,((oosR$period-1)/8-16)==i&oosR$strata==j)
		if (dim(temp)[1]>0){
			tw_m_oos[i,j,1]<-sum(subset(temp$Number,temp$winmon<3))
			tw_m_oos[i,j,2]<-sum(subset(temp$Number,temp$winmon>3))
			}}
	for (j in 1:Nstrata){
		w_m_oos[i,j,1]<-sum(tw_m_oos[i,gridst_end[j,1]:gridst_end[j,2],1])
		w_m_oos[i,j,2]<-sum(tw_m_oos[i,gridst_end[j,1]:gridst_end[j,2],2])
		}}
##
keep<-match(c("Date.Collected","RM_Start","Haul","Habitat","Effort_m.2","year","month",
"Species_Codes","Gear","Length..SL.mm.","AgeClass","SumOfSPEC"),names(oosc))
oosc<-oosc[,keep]
oosc<-subset(oosc,oosc$Habitat!=""&oosc$Gear!="larval")
ghab<-c("","BW","MCPLPO","MCPO","MCSHPLPO","MCSHPO","PO","SCPLPO","SCPO","SCSHPLPO","SCSHPO","SHPO","MCED","MCSHED","SCED","SCSHED")
oosc$cHab<-ifelse(is.na(match(oosc$Habitat,ghab))==FALSE,1,3)
oosc$strata<-findInterval(oosc$RM_Start,vgrid)
oosc$period<-8*(oosc$year-2002)+oosc$month-3
extractdom<-function(x){
	lx<-nchar(x)
	tx<-numeric()
	for (i in 1:lx){
		tx[i]<-substr(x,i,i)}
	wx<-which(tx=="/")
	as.numeric(substr(x,wx[1]+1,wx[2]-1))}
oosc$dom<-NA
oosc$jul<-NA
oosc$dry<-NA
for (j in 1:length(oosc[,1])){
	oosc$dom[j]<-extractdom(paste(oosc$Date.Collected[j]))
	oosc$jul[j]<-julian(as.Date(paste(oosc$year[j],oosc$month[j],oosc$dom[j],sep="-")),as.Date(paste(oosc$year[j],"03","31",sep="-")))[[1]]
	tf1<-fdryday_oos[(oosc$year[j]-2018),oosc$strata[j]]
	tl1<-ldryday_oos[(oosc$year[j]-2018),oosc$strata[j]]
	oosc$dry[j]<-ifelse(is.na(tf1)==TRUE|(tf1>oosc$jul[j])|(tl1<oosc$jul[j]),0,1)}
oosc$uni_id<-paste(oosc$year,oosc$jul,oosc$RM_Start,oosc$Haul)
#
toosc<-data.frame(uni=sort(unique(oosc$uni_id)))
toosc$RMStart<-oosc$RM_Start[match(toosc[,1],oosc$uni_id)]
toosc$period<-oosc$period[match(toosc[,1],oosc$uni_id)]
toosc$jul<-oosc$jul[match(toosc[,1],oosc$uni_id)]
toosc$year<-oosc$year[match(toosc[,1],oosc$uni_id)]
toosc$month<-oosc$month[match(toosc[,1],oosc$uni_id)]
toosc$date<-oosc$Date.Collected[match(toosc[,1],oosc$uni_id)]
toosc$angQ<-angQ$cfs[match(toosc$date,angQ$Date)]
toosc$sanaQ<-sanaQ$cfs[match(toosc$date,sanaQ$Date)]
toosc$HaulNo<-oosc$Haul[match(toosc[,1],oosc$uni_id)]
toosc$type<-ifelse(oosc$cHab[match(toosc[,1],oosc$uni_id)]==1,1,2)
toosc$effort<-as.numeric(paste(oosc$Effort_m.2[match(toosc[,1],oosc$uni_id)]))
toosc$strata<-oosc$strata[match(toosc[,1],oosc$uni_id)]
toosc$cQ<-ifelse(toosc$strata>573,toosc$angQ,toosc$sanaQ)
toosc$hybama<-0
thybama<-subset(oosc,oosc$Species_Codes=="HYBAMA")
for (i in 1:length(toosc[,1])){
	temp<-subset(thybama,thybama$uni_id==toosc$uni[i])
	toosc$hybama[i]<-sum(temp$SumOfSPEC)}
toosc$Cstrata<-crossgrid[toosc$strata]
toosc$c_id<-paste(toosc$year,toosc$jul,toosc$Cstrata,toosc$type)
##
oos_C<-data.frame(cid=sort(unique(toosc$c_id)))
oos_C$jul<-toosc$jul[match(oos_C[,1],toosc$c_id)]
oos_C$year<-toosc$year[match(oos_C[,1],toosc$c_id)]
oos_C$month<-toosc$month[match(oos_C[,1],toosc$c_id)]
oos_C$Cstrata<-toosc$Cstrata[match(oos_C[,1],toosc$c_id)]
oos_C$type<-toosc$type[match(oos_C[,1],toosc$c_id)]
oos_C$effort<-0
oos_C$cQ<-0
oos_C$hybama<-0
oos_C$period<-toosc$period[match(oos_C[,1],toosc$c_id)]
mean_names<-c("cQ")
sum_names<-c("effort","hybama")
tsm<-match(sum_names,names(toosc))
sm<-match(sum_names,names(oos_C))
tmn<-match(mean_names,names(toosc))
smn<-match(mean_names,names(oos_C))
for (i in 1:length(oos_C[,1])){
	temp<-subset(toosc[,tsm],toosc$c_id==oos_C$cid[i])
	oos_C[i,sm]<-colSums(temp)
	temp<-subset(toosc[,tmn],toosc$c_id==oos_C$cid[i])
	oos_C[i,smn]<-mean(temp)
	}
Nobs_oosC<-dim(oos_C)[1]
# function to make out of sample forecasts
forecast_oos_re<-function(mod,Xout){
	Sp_N0<-extract(mod,"Sp_N")[[1]][,18,,]#last dimension is age / aug
	irphi<-extract(mod,"irphi")[[1]]
	beta_2<-extract(mod,"beta_2")[[1]]
	beta_stk<-extract(mod,"beta_stk")[[1]]
	a<-extract(mod,"a")[[1]]
	mu_lM0<-extract(mod,"mu_lM0")[[1]]
	mu_lM1<-extract(mod,"mu_lM1")[[1]]
	mu_lMw<-extract(mod,"mu_lMw")[[1]]
	sd_lM<-extract(mod,"sd_lM")[[1]]
	move<-extract(mod,"move")[[1]]
	rp0<-extract(mod,"rp0")[[1]]
	rp1<-extract(mod,"rp1")[[1]]
	mu_lbeta<-extract(mod,"mu_lbeta")[[1]]
	sd_lbeta<-extract(mod,"sd_lbeta")[[1]]
	B_lbeta<-extract(mod,"B_lbeta")[[1]]
	A0_perpool<-extract(mod,"A0_perpool")[[1]]
	AtQ_perpool<-extract(mod,"AtQ_perpool")[[1]]
	sl_width<-extract(mod,"sl_width")[[1]]
	bankfull<-extract(mod,"bankfull")[[1]]
	alpha0_max<-extract(mod,"alpha0_max")[[1]]
	alpha1_max<-extract(mod,"alpha1_max")[[1]]
	alpha0_int<-extract(mod,"alpha0_int")[[1]]
	alpha1_int<-extract(mod,"alpha1_int")[[1]]
	p0<-extract(mod,"p0")[[1]]
	p1<-extract(mod,"p1")[[1]]
	sz<-extract(mod,"sz")[[1]]
	sd_lbeta<-extract(mod,"sd_lbeta")[[1]]
	iter<-length(p0)
	#
	Mw<-array(NA,dim=c(iter,3,2))
	M0<-array(NA,dim=c(iter,3,2))
	M1<-array(NA,dim=c(iter,3,2))
	for (k in 1:Nstrata){
		Mw[,k,1]<-exp(rnorm(iter,mu_lMw[,k],sd_lM))
		Mw[,k,2]<-exp(rnorm(iter,mu_lMw[,k],sd_lM))
		M0[,k,1]<-exp(rnorm(iter,mu_lM0[,k],sd_lM))
		M0[,k,2]<-exp(rnorm(iter,mu_lM0[,k],sd_lM))
		M1[,k,1]<-exp(rnorm(iter,mu_lM1[,k],sd_lM))
		M1[,k,2]<-exp(rnorm(iter,mu_lM1[,k],sd_lM))
		}
	#
	FN<-array(NA,dim=c(iter,2,2,3))# iter, years,size classes, river segments		
	pC<-matrix(NA,ncol=Nobs_oosC,nrow=iter)
	C<-matrix(NA,ncol=Nobs_oosC,nrow=iter)
	SpN<-array(NA,dim=c(iter,2,3,3))# iter, years,size classes, river segments
	effS<-array(NA,dim=c(iter,2,3))# iter, years, river segments
	for (k in 1:Nstrata){
		SpN[,1,1,k]<-Sp_N0[,k,1]
		SpN[,1,2,k]<-Sp_N0[,k,2]
		SpN[,1,3,k]<-(w_m_oos[1,k,1]*exp(-60*Mw[,k,1])+w_m_oos[1,k,2]*exp(-151*Mw[,k,1]))*irphi
		effS[,1,k]=SpN[,1,1,k]+SpN[,1,2,k]*beta_2+SpN[,1,3,k]*beta_stk
		tN1<-rowSums(SpN[,1,1:2,k])
		tRf<-exp(rnorm(iter,mu_lbeta[,k],sd_lbeta)+B_lbeta*Xout[1,k])
		tN0=a*(effS[,1,k])/(1+a*(effS[,1,k])/tRf)		
		FN[,1,1,k]= tN0*exp(-M0[,k,1]*124)*((1-cum_nd_oos[1,91,k])+(move-1)*(cum_nd_oos[1,215,k]-cum_nd_oos[1,91,k])+(1-move)*rp0*(cum_phiR_oos[1,215,k]-cum_phiR_oos[1,91,k]))
		FN[,1,2,k]=tN1*exp(-M1[,k,1]*215)*(1+(move-1)*cum_nd_oos[1,215,k]+(1-move)*rp1*cum_phiR_oos[1,215,k])			
		#
		SpN[,2,1:3,k]=cbind(FN[,1,1,k]*exp(-Mw[,k,2]*150),FN[,1,2,k]*exp(-Mw[,k,2]*150),(w_m_oos[2,k,1]*exp(-60*Mw[,k,2])+w_m_oos[2,k,2]*exp(-151*Mw[,k,2]))*irphi)
		effS[,2,k]=SpN[,2,1,k]+SpN[,2,2,k]*beta_2+SpN[,2,3,k]*beta_stk
		tN1<-rowSums(SpN[,2,1:2,k])
		tRf<-exp(rnorm(iter,mu_lbeta[,k],sd_lbeta)+B_lbeta*Xout[2,k])
		tN0=a*(effS[,2,k])/(1+a*(effS[,2,k])/tRf)		
		FN[,2,1,k]= tN0*exp(-M0[,k,2]*124)*((1-cum_nd_oos[2,91,k])+(move-1)*(cum_nd_oos[2,215,k]-cum_nd_oos[2,91,k])+(1-move)*rp0*(cum_phiR_oos[2,215,k]-cum_phiR_oos[2,91,k]))
		FN[,2,2,k]=tN1*exp(-M1[,k,2]*215)*(1+(move-1)*cum_nd_oos[2,215,k]+(1-move)*rp1*cum_phiR_oos[2,215,k])			
		}			
	for (i in 1:Nobs_oosC){
		q<-oos_C$cQ[i]/1000
		totpool=exp(A0_perpool+AtQ_perpool*q)*200*sl_width[oos_C[i,5]]*q/(1+sl_width[oos_C[i,5]]*q/bankfull[oos_C[i,5]])
		totrun=(1-exp(A0_perpool+AtQ_perpool*q))*200*sl_width[oos_C[i,5]]*q/(1+sl_width[oos_C[i,5]]*q/bankfull[oos_C[i,5]])
		talpha0=c(alpha0_int*q/(1+alpha0_int*q/alpha0_max),1)
		talpha1=c(alpha1_int*q/(1+alpha1_int*q/alpha1_max),1)
		y<-ifelse(oos_C[i,3]==2019,1,2)
		pC[,i]=(p0*oos_C$effort[i]*talpha0[oos_C$type[i]]*.2*FN[,y,1,oos_C[i,5]]*exp(-M0[,oos_C[i,5],(oos_C[i,3]-2018)]*(oos_C[i,2]-91))*
			((1-cum_nd_oos[y,91,oos_C[i,5]])+(move-1)*(cum_nd_oos[y,oos_C[i,2],oos_C[i,5]]-cum_nd_oos[y,91,oos_C[i,5]])+
			(1-move)*rp0*(cum_phiR_oos[y,oos_C[i,2],oos_C[i,5]]-cum_phiR_oos[y,91,oos_C[i,5]])))/
			(StrataLen_oos[y,oos_C[i,2],oos_C[i,5]]*(talpha0[1]*totpool+totrun))+
			(p1*oos_C$effort[i]*talpha1[oos_C$type[i]]*.2*FN[,y,2,oos_C[i,5]]*
			exp(-M1[,oos_C[i,5],(oos_C[i,3]-2018)]*oos_C[i,2])*(1+(move-1)*cum_nd_oos[y,oos_C[i,2],oos_C[i,5]]+(1-move)*rp1*cum_phiR_oos[y,oos_C[i,2],oos_C[i,5]]))/
			(StrataLen_oos[y,oos_C[i,2],oos_C[i,5]]*(talpha1[1]*totpool+totrun))
		C[,i]<-rnbinom(iter,mu=pC[,i],size=sz)
		}
	w19<-which(oos_C$year==2019)
	pc19<-rowSums(pC[,w19])
	ef19<-sum(oos_C$effort[w19])
	pc20<-rowSums(pC[,-w19])
	ef20<-sum(oos_C$effort[-w19])
	c19<-rowSums(C[,w19])
	c20<-rowSums(C[,-w19])
	return(list(ef19=ef19,ef20=ef20,pc19=pc19,pc20=pc20,c19=c19,c20=c20,SpN=SpN,FN=FN,effS=effS))}
#
aQ_out<-cbind(subset(angQ[,5],angQ[,2]==2019&angQ[,3]>3&angQ[,3]<10),subset(angQ[,5],angQ[,2]==2020&angQ[,3]>3&angQ[,3]<10))
sQ_out<-cbind(subset(sanaQ[,5],sanaQ[,2]==2019&sanaQ[,3]>3&sanaQ[,3]<10),subset(sanaQ[,5],sanaQ[,2]==2020&sanaQ[,3]>3&sanaQ[,3]<10))
preds_out<-array(NA,dim=c(2,3))
temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
t1<-numeric()
for (t in 1:6){
	se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
	t1[t]<-temp[t,6]-se
	se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
	t1[(t+7)]<-temp[(t+6),6]+se
	}
t1[7]<-1
t1[14]<-0
temp<-subset(ee[[3]],ee[[3]][,1]=="2")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t2<-temp[1,6]
temp<-subset(ee[[3]],ee[[3]][,1]=="4")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t4<-temp[1,6]-se
for (r in 1:3){
	temp<-subset(ee[[3]],ee[[3]][,1]==riversegment[r]&is.na(ee[[3]][,5])==FALSE)
	t3<-numeric()
	for (t in 1:6){
		se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
		t3[t]<-temp[t,6]+se
		}	
	for (t in 8:15){
		se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
		t3[t]<-temp[t,6]-se
		}
	t3[7]<-temp[7,6]
	for (j in 1:2){
		if (r==3) {q<-aQ_out[,j]} else {q<-sQ_out[,j]}
		preds_out[j,r]<-calc_cov(q,t3,t2,t4,t1) # only uses present habitat
		}}
#### out of sample forecast from integrated model with larval carrying capacity
f1<-forecast_oos_re(M2_1re,preds_out-mean(preds2))
# out of sample forecasts from models based on may-june average flow
mjflow_sana_out<-subset(sanaQ,(sanaQ$month==5|sanaQ$month==6)&sanaQ$year>2018)
mjflow_ang_out<-subset(angQ,(angQ$month==5|angQ$month==6)&angQ$year>2018)
simpflow_out<-cbind(tapply(mjflow_sana_out$cfs,mjflow_sana_out$year,mean),tapply(mjflow_sana_out$cfs,mjflow_sana_out$year,mean),tapply(mjflow_ang_out$cfs,mjflow_ang_out$year,mean))/1000
f2<-forecast_oos_re(M2_2re,simpflow_out-mean(simpflow))
# calculate observation
obs_catch<-tapply(oos_C$hybama,oos_C$year,sum)
# calculate prediction using biop
tbQ19<-mean(subset(angQ[,5],(angQ[,3]==5|angQ[,3]==6)&angQ[,2]==2019))
tbQ20<-mean(subset(angQ[,5],(angQ[,3]==5|angQ[,3]==6)&angQ[,2]==2020))
biop19<-(10^(rnorm(1000000,-0.1477 + (0.0004*tbQ19)-(tbQ19^2) * 0.000000014284,.1)) - 1)*f1$ef19/100 #approximate sd of .1
biop20<-(10^(rnorm(1000000,-0.1477 + (0.0004*tbQ20)-(tbQ20^2) * 0.000000014284,.1)) - 1)*f1$ef20/100 #approximate sd of .1
# out of sample figure
q025<-function(x){quantile(x,.025)}
q975<-function(x){quantile(x,.975)}
q10<-function(x){quantile(x,.1)}
q90<-function(x){quantile(x,.9)}

# code to make Fig 3
par(mfrow=c(2,1))
par(mar=c(4,4,1,1))
plot(c(mean(f1$c19),mean(f2$c19),mean(biop19)),pch=19,xlab="",ylab="",ylim=c(0,7500),axes=FALSE,xlim=c(0.5,3.5))
segments(c(1:3),c(q025(f1$c19),q025(f2$c19),q025(biop19)),c(1:3),c(q975(f1$c19),q975(f2$c19),q975(biop19)),lwd=2)
segments(c(1:3),c(q10(f1$c19),q10(f2$c19),q10(biop19)),c(1:3),c(q90(f1$c19),q90(f2$c19),q90(biop19)),lwd=4)
segments(0.5,obs_catch[1],3.75,obs_catch[1],lwd=2,lty=2,"red")
text(1.5,obs_catch[1]+500,"Observed",col="red")
mtext("2019 catch",2,2.5)
axis(2,c(0,2500,5000,7500),pos=0.5,las=T)
axis(1,at=c(0.5,1,2,3,4),c("","Integrated model","Integrated model","Regression model",""),pos=0)
axis(1,at=c(0.5,1,2,3,4),c("","with expert","with may through","with may through",""),pos=0,padj=1.5)
axis(1,at=c(0.5,1,2,3,4),c("","informed covariate","june discharge","june discharge",""),pos=0,padj=3)
text(0.05,8000,"A)",xpd=T)
#second panel
plot(c(mean(f1$c20),mean(f2$c20),mean(biop20)),pch=19,xlab="",ylab="",ylim=c(0,300),axes=FALSE,xlim=c(0.5,3.5))
segments(c(1:3),c(q025(f1$c20),q025(f2$c20),q025(biop20)),c(1:3),c(q975(f1$c20),q975(f2$c20),q975(biop20)),lwd=2)
segments(c(1:3),c(q10(f1$c20),q10(f2$c20),q10(biop20)),c(1:3),c(q90(f1$c20),q90(f2$c20),q90(biop20)),lwd=4)
segments(0.5,obs_catch[2],3.75,obs_catch[2],lwd=2,lty=2,"red")
text(1.5,obs_catch[2]+20,"Observed",col="red")
mtext("2020 catch",2,2.25)
axis(2,c(0,100,200,300),pos=0.5,las=T)
axis(1,at=c(0.5,1,2,3,4),c("","Integrated model","Integrated model","Regression model",""),pos=0)
axis(1,at=c(0.5,1,2,3,4),c("","with expert","with may through","with may through",""),pos=0,padj=1.5)
axis(1,at=c(0.5,1,2,3,4),c("","informed covariate","june discharge","june discharge",""),pos=0,padj=3)
text(0.05,320,"B)",xpd=T)

# calculations for, and code to plot, figure 4
aug<-c(rowSums(w_um),rowSums(w_m))[1:17]
sN<-extract(M2_1re,"Sp_N")[[1]]
N1<-extract(M2_1re,"N1")[[1]]
N0<-extract(M2_1re,"N0")[[1]]
fN<-extract(M2_1re,"F_N")[[1]]
fn<-fN[,,1,1]+fN[,,2,1]+fN[,,3,1]+fN[,,1,2]+fN[,,2,2]+fN[,,3,2]
n1<-N1[,-18,1]+N1[,-18,2]+N1[,-18,3]
sn<-sN[,,1,1]+sN[,,2,1]+sN[,,3,1]+sN[,,1,2]+sN[,,2,2]+sN[,,3,2]+sN[,,1,3]+sN[,,2,3]+sN[,,3,3]
par(mfrow=c(3,1))
par(mar=c(3,3,1,1))
plot(c(2002:2018),apply(fn,2,mean),pch=19,col="white",axes=FALSE,ylab="",xlab="",log="y",ylim=c(100,100000000),xlim=c(2001.5,2018.5))
segments(c(2002:2018)+.15,apply(fn,2,q975),c(2002:2018)+.15,apply(fn,2,q025),col="orange",lwd=6)
segments(c(2002:2018)-.15,apply(sn[,-18]-n1,2,q975),c(2002:2018)-.15,apply(sn[,-18]-n1,2,q025),col=rgb(0,.5,1,.4),lwd=6)
segments(c(2002:2018)-.15,apply(n1,2,q975),c(2002:2018)-.15,apply(n1,2,q025),col="green",lwd=6)
points(c(2008:2011)+.15,exp(lNz),pch=19,col=rgb(0,0,0,.3))
segments(c(2008:2011)+.15,exp(lNz+lCz*1.96),c(2008:2011)+.15,exp(lNz-lCz*1.96),col=rgb(0,0,0,.4),lwd=3)
axis(2,c(100,1000,10000,100000,1000000,10000000),c("10^2","10^3","10^4","10^5","10^6","10^7"),las=T,pos=2001.5)
axis(1,seq(2002,2018,4),c(seq(2002,2018,4)),pos=100)
axis(1,c(2001.5,2018.5),c("",""),tck=FALSE,pos=100)
mtext("Total Abundance",2,1.5,cex=0.7)
mtext("Year",1,2,cex=0.7)
legend(2001.75,100000000,"Spring - wild born" ,fill="green",bty="n",border=NA)
legend(2006.25,100000000,"Spring - augmented" ,fill=rgb(0,.5,1,.4),bty="n",border=NA)
legend(2010.75,100000000,"Fall - wild born" ,fill="orange",bty="n",border=NA)
legend(2014.25,100000000,"Population estimate" ,col=rgb(0,0,0,.3),pch=19,bty="n",border=NA)
text(2000,100000000,"A)",xpd=T)
#2nd panel
plot(c(2002:2018),rep(0.5,17),pch=19,col="white",axes=FALSE,ylab="",xlab="",ylim=c(0,1.3),xlim=c(2001.5,2018.5))
segments(c(2002:2018)+.15,0,c(2002:2018)+.15,apply((fN[,,1,1]+fN[,,2,1]+fN[,,3,1])/fn,2,mean),col="orange",lwd=3)
segments(c(2002:2018)+.15,apply((fN[,,1,1]+fN[,,2,1]+fN[,,3,1])/fn,2,mean),c(2002:2018)+.15,1,col="red",lwd=3)
segments(c(2002:2018)-.15,0,c(2002:2018)-.15,apply((sN[,-18,1,1]+sN[,-18,2,1]+sN[,-18,3,1])/sn[,-18],2,mean),col=rgb(0,1,0,.7),lwd=3)
segments(c(2002:2018)-.15,apply((sN[,-18,1,1]+sN[,-18,2,1]+sN[,-18,3,1])/sn[,-18],2,mean),c(2002:2018)-.15,apply((n1)/sn[,-18],2,mean),col="darkgreen",lwd=3)
segments(c(2002:2018)-.15,apply((n1)/sn[,-18],2,mean),c(2002:2018)-.15,1,col=rgb(0,.5,1,.6),lwd=3)
axis(1,c(2001.5,2018.5),c("",""),tck=FALSE,pos=0)
axis(1,seq(2002,2018,4),seq(2002,2018,4),pos=0)
axis(2,c(0,0.5,1),las=T,pos=2001.5)
t2<-subset(t2mon,t2mon$month<6&t2mon$VC==1&is.na(t2mon$LengthSL)==FALSE)
mtext("Population structure",2,1.5,cex=0.7)
mtext("Year",1,2,cex=0.7)
legend(2002,1.4,"Spring - Age 1" ,fill=rgb(0,1,0,.7),bty="n",border=NA)
legend(2002,1.25,"Spring - Age 2+" ,fill="darkgreen",bty="n",border=NA)
legend(2008,1.4,"Spring - Augmented" ,fill=rgb(0,.5,1,.6),bty="n",border=NA)
legend(2008,1.25,"Fall - Age 0" ,fill="orange",bty="n",border=NA)
legend(2014.2,1.4,"Fall - Age 1+" ,fill="red",bty="n",border=NA)
text(2000.05,1.3,"B)",xpd=T)
# 3rd panel
phi0_nat<-exp(-124*extract(M2_1re,"M0")[[1]])
phi1_nat<-exp(-215*extract(M2_1re,"M1")[[1]])
phiW<-exp(-150*extract(M2_1re,"Mw")[[1]])
irphi<-extract(M2_1re,"irphi")[[1]]
phi_wdry1<-fN[,,,2]/N1[,-18,]
phi_wdry0<-fN[,,,1]/N0
par(mar=c(3,4,1,1))
xz<-c(1,1.2,2,2.2,3,3.2,5,5.2,6,6.2,7,7.2,9.25,10,10.75,12)
colz<-c("orange","orange","red","red","purple","purple","orange","orange","red","red","purple","purple","orange","red","purple","darkgrey")
plot(xz,c(mean(apply(phi0_nat[,1,],1,mean)),mean(apply(phi_wdry0[,,1],1,mean),xlim=c(0,13.5)),
	mean(apply(phi0_nat[,2,],1,mean)),mean(apply(phi_wdry0[,,2],1,mean)),
	mean(apply(phi0_nat[,3,],1,mean)),mean(apply(phi_wdry0[,,3],1,mean)),
	mean(apply(phi1_nat[,1,],1,mean)),mean(apply(phi_wdry1[,,1],1,mean)),
	mean(apply(phi1_nat[,2,],1,mean)),mean(apply(phi_wdry1[,,2],1,mean)),
	mean(apply(phi1_nat[,3,],1,mean)),mean(apply(phi_wdry1[,,3],1,mean)),
	mean(apply(phiW[,1,],1,mean)),mean(apply(phiW[,2,],1,mean)),
	mean(apply(phiW[,3,],1,mean)),mean(irphi)),pch=c(rep(c(1,19),6),19,19,19,19),col=colz,axes=FALSE,xlab="",ylab="",ylim=c(0.01,3.5),log="y")
segments(xz,c(q025(apply(phi0_nat[,1,],1,mean)),q025(apply(phi_wdry0[,,1],1,mean)),
	q025(apply(phi0_nat[,2,],1,mean)),q025(apply(phi_wdry0[,,2],1,mean)),
	q025(apply(phi0_nat[,3,],1,mean)),q025(apply(phi_wdry0[,,3],1,mean)),
	q025(apply(phi1_nat[,1,],1,mean)),q025(apply(phi_wdry1[,,1],1,mean)),
	q025(apply(phi1_nat[,2,],1,mean)),q025(apply(phi_wdry1[,,2],1,mean)),
	q025(apply(phi1_nat[,3,],1,mean)),q025(apply(phi_wdry1[,,3],1,mean)),
	q025(apply(phiW[,1,],1,mean)),q025(apply(phiW[,2,],1,mean)),
	q025(apply(phiW[,3,],1,mean)),q025(irphi)),
	xz,c(q975(apply(phi0_nat[,1,],1,mean)),q975(apply(phi_wdry0[,,1],1,mean)),
	q975(apply(phi0_nat[,2,],1,mean)),q975(apply(phi_wdry0[,,2],1,mean)),
	q975(apply(phi0_nat[,3,],1,mean)),q975(apply(phi_wdry0[,,3],1,mean)),
	q975(apply(phi1_nat[,1,],1,mean)),q975(apply(phi_wdry1[,,1],1,mean)),
	q975(apply(phi1_nat[,2,],1,mean)),q975(apply(phi_wdry1[,,2],1,mean)),
	q975(apply(phi1_nat[,3,],1,mean)),q975(apply(phi_wdry1[,,3],1,mean)),
	q975(apply(phiW[,1,],1,mean)),q975(apply(phiW[,2,],1,mean)),
	q975(apply(phiW[,3,],1,mean)),q975(irphi)),col=colz)
legend(1,5,c("San Acacia - w/o drying","San Acacia - w/ drying"),pch=c(1,19),col=c("orange","orange"),bty="n",border=NA)
legend(4.25,5,c("Isleta - w/o drying","Isleta - w/ drying"),pch=c(1,19),col=c("red","red"),bty="n",border=NA)
legend(6.75,5,c("Angostura - w/o drying","Angostura - w/ drying"),pch=c(1,19),col=c("purple","purple"),bty="n",border=NA)
legend(10,4,c("All river segments"),pch=19,col="darkgrey",bty="n",border=NA)
mtext("Survival",2,2.5,cex=0.7)
axis(1,c(0,2,6,10,12,13.5),c("","Age - 0 (Jul-Nov)","Age - 1+ (Apr-Nov)","Age - 1+ (Nov-Apr)","Stocking",""),pos=0.01,tck=FALSE)
axis(2,c(.01,.05,.1,.5,1),las=T,pos=0.75)
text(-.2,3.5,"C)",xpd=T)

# code to make calculations for Fig 5
A0<-extract(M2_1re,"A0_perpool")[[1]]
AtQ<-extract(M2_1re,"A0_perpool")[[1]]
a0_int<-extract(M2_1re,"alpha0_int")[[1]]
a0_max<-extract(M2_1re,"alpha0_max")[[1]]
a1_int<-extract(M2_1re,"alpha1_int")[[1]]
a1_max<-extract(M2_1re,"alpha1_max")[[1]]
sl_width<-extract(M2_1re,"sl_width")[[1]]
bankfull<-extract(M2_1re,"bankfull")[[1]]
p0<-extract(M2_1re,"p0")[[1]]
p1<-extract(M2_1re,"p1")[[1]]
sz<-extract(M2_1re,"sz")[[1]]
ppool<-array(NA,dim=c(1000,15000))
thab<-array(NA,dim=c(1000,3,15000))
tpool<-array(NA,dim=c(1000,3,15000))
trun<-array(NA,dim=c(1000,3,15000))
ta0<-array(NA,dim=c(1000,15000))
ta1<-array(NA,dim=c(1000,15000))
catch<-array(NA,dim=c(1000,6,15000))
cpe100<-array(NA,dim=c(1000,15000))
cpe100_lp<-array(NA,dim=c(1000,15000))
##assume 5K  m2 of effort in each segment with 20% in pools & 10K fish in each segment
for (i in 1:1000){
	ppool[i,]<-exp(A0+AtQ*(i/1000))
	thab[i,1,]<-200*sl_width[,1]*(i/1000)/(1+sl_width[,1]*(i/1000)/bankfull[,1])
	thab[i,2,]<-200*sl_width[,2]*(i/1000)/(1+sl_width[,2]*(i/1000)/bankfull[,2])
	thab[i,3,]<-200*sl_width[,3]*(i/1000)/(1+sl_width[,3]*(i/1000)/bankfull[,3])
	tpool[i,1,]<-1000*92.1*exp(A0+AtQ*(i/1000))*sl_width[,1]*(i/1000)/(1+sl_width[,1]*(i/1000)/bankfull[,1])
	trun[i,1,]<-1000*92.1*(1-exp(A0+AtQ*(i/1000)))*sl_width[,1]*(i/1000)/(1+sl_width[,1]*(i/1000)/bankfull[,1])
	tpool[i,2,]<-1000*85.5*exp(A0+AtQ*(i/1000))*sl_width[,2]*(i/1000)/(1+sl_width[,2]*(i/1000)/bankfull[,2])
    trun[i,2,]<-1000*85.5*(1-exp(A0+AtQ*(i/1000)))*sl_width[,2]*(i/1000)/(1+sl_width[,2]*(i/1000)/bankfull[,2])
    tpool[i,3,]<-1000*65*exp(A0+AtQ*(i/1000))*sl_width[,3]*(i/1000)/(1+sl_width[,3]*(i/1000)/bankfull[,3])
    trun[i,3,]<-1000*65*(1-exp(A0+AtQ*(i/1000)))*sl_width[,3]*(i/1000)/(1+sl_width[,3]*(i/1000)/bankfull[,3])
	ta0[i,]<-a0_int*(i/1000)/(1+a0_int*(i/1000)/a0_max)
	ta1[i,]<-a1_int*(i/1000)/(1+a1_int*(i/1000)/a1_max)
	#San Acacia - pool
	tmu<-.9*p0*5000*.2*ta0[i,]*10000/(ta0[i,]*tpool[i,1,]+trun[i,1,])+
		.1*p1*5000*.2*ta1[i,]*10000/(ta1[i,]*tpool[i,1,]+trun[i,1,])
	catch[i,1,]<-rnbinom(15000,tmu,sz)
	#San Acacia - run/riffle
	tmu<-.9*p0*5000*.8*10000/(ta0[i,]*tpool[i,1,]+trun[i,1,])+
		.1*p1*5000*.8*10000/(ta1[i,]*tpool[i,1,]+trun[i,1,])
	catch[i,2,]<-rnbinom(15000,tmu,sz)
	#Isleta
	tmu<-.9*p0*5000*.2*ta0[i,]*10000/(ta0[i,]*tpool[i,2,]+trun[i,2,])+
		.1*p1*5000*.2*ta1[i,]*10000/(ta1[i,]*tpool[i,2,]+trun[i,2,])
	catch[i,3,]<-rnbinom(15000,tmu,sz)
	#
	tmu<-.9*p0*5000*.8*10000/(ta0[i,]*tpool[i,2,]+trun[i,2,])+
		.1*p1*5000*.8*10000/(ta1[i,]*tpool[i,2,]+trun[i,2,])
	catch[i,4,]<-rnbinom(15000,tmu,sz)
	#Angostura
	tmu<-.9*p0*5000*.2*ta0[i,]*10000/(ta0[i,]*tpool[i,3,]+trun[i,3,])+
		.1*p1*5000*.2*ta1[i,]*10000/(ta1[i,]*tpool[i,3,]+trun[i,3,])
	catch[i,5,]<-rnbinom(15000,tmu,sz)
	#
	tmu<-.9*p0*5000*.8*10000/(ta0[i,]*tpool[i,3,]+trun[i,3,])+
		.1*p1*5000*.8*10000/(ta1[i,]*tpool[i,3,]+trun[i,3,])
	catch[i,6,]<-rnbinom(15000,tmu,sz)
	#
	cpe100[i,]<-apply(catch[i,,],2,sum)/150
	cpe100_lp[i,]<-((catch[i,1,]+catch[i,3,]+catch[i,5,])*.5+(catch[i,2,]+catch[i,4,]+catch[i,6,])*9/8)/150
	}
# code to plot figure 5
par(mfrow=c(2,2))
par(mar=c(3,4,2,2))
plot(1,col="white",xlim=c(1,1000),ylim=c(0,350),xlab="",ylab="",axes=FALSE,log="x")
mtext("Wetted area per site (x 100 m^2)",2,2.25,cex=0.7)
mtext("Discharge(cfs)",1,1.5,cex=0.7)
axis(1,c(1,10,100,1000),pos=0)
axis(2,c(0,100,200,300),las=T,pos=1)
polygon(c(1:1000,1000:1),c(apply(ppool[1:1000,],1,q025),apply(ppool[1000:1,],1,q975))*4000,col=rgb(0,0,1,.4),border=NA)
polygon(c(1:1000,1000:1),c(apply(thab[1:1000,1,],1,q025),apply(thab[1000:1,1,],1,q975))/100,col="orange",border="orange")
polygon(c(1:1000,1000:1),c(apply(thab[1:1000,2,],1,q025),apply(thab[1000:1,2,],1,q975))/100,col="red",border="red")
polygon(c(1:1000,1000:1),c(apply(thab[1:1000,3,],1,q025),apply(thab[1000:1,3,],1,q975))/100,col="purple",border="purple")
legend(1,375,c("% Pool (2nd y-axis)","San Acacia"),fill=c(rgb(0,0,1,.4),col="orange"),bty="n",border=NA)
legend(50,375,c("Isleta","Angostura"),fill=c("red","purple"),bty="n",border=NA)
axis(4,at=seq(0,300,100),labels=seq(0,7.5,2.5),las=T,pos=1000,col="blue",col.axis="blue")
mtext("% pool",4,1.5,col="blue",las=FALSE,cex=0.7)
text(.22,390,"A)",xpd=T)
## 2nd panel
plot(1,col="white",xlim=c(1,1000),ylim=c(0,15),xlab="",ylab="",axes=FALSE,log="x")
polygon(c(1:1000,1000:1),c(apply(ta0[1:1000,],1,q025),apply(ta0[1000:1,],1,q975)),col=rgb(0,1,0,.5),border=NA)
polygon(c(1:1000,1000:1),c(apply(ta1[1:1000,],1,q025),apply(ta1[1000:1,],1,q975)),col=rgb(1,0,1,.5),border=NA)
legend(5,15,c("Age-0","Age-1+"),fill=c(rgb(0,1,0,0.5),rgb(1,0,1,.5)),bty="n",border=NA)
mtext("Density in pools relative to runs/riffles",2,1.5,cex=0.7)
mtext("Discharge(cfs)",1,1.5,cex=0.7)
axis(1,c(1,10,100,1000),pos=0)
axis(2,las=T,pos=1)
segments(1,1,1000,1,lty=2)
text(.24,16.5,"B)",xpd=T)
# 3rd panel
plot(1,col="white",xlim=c(1,1000),ylim=c(0.1,100),xlab="",ylab="",axes=FALSE,log="xy")
polygon(c(1:1000,1000:1),c(apply(cpe100[1:1000,],1,q025),apply(cpe100[1000:1,],1,q975)),col=rgb(0,0,0,.4),border=NA)
mtext("Catch per 100 m^2",2,2.25,cex=0.7)
mtext("Discharge(cfs)",1,1.5,cex=0.7)
axis(1,c(1,10,100,1000),pos=0.1)
axis(2,c(0.1,1,10,100),las=T,pos=1)
text(.22,210,"C)",xpd=T)
# 4th panel
plot(c(7.5,10,15,30,60,90,120)/30*mean(cpe100[500,]),ylim=c(0,4),xlim=c(0.5,7.5),xlab="",ylab="",axes=FALSE,col="white")
mtext("Catch per 100 m^2",2,1.5,cex=0.7)
mtext("Overall abundance (x 1000)",1,1.5,cex=0.7)
for (i in 1:9){	
	polygon(c(i-.1,i+.1,i+.1,i-.1)-.15,c(7.5,10,15,30,60,90,120)[i]/30*c(q025(cpe100[500,]),q025(cpe100[500,]),q975(cpe100[500,]),q975(cpe100[500,])),col=rgb(0,0,0,.4),border=NA)
	polygon(c(i-.1,i+.1,i+.1,i-.1)+.15,c(7.5,10,15,30,60,90,120)[i]/30*c(q025(cpe100_lp[500,]),q025(cpe100_lp[500,]),q975(cpe100_lp[500,]),q975(cpe100_lp[500,])),col=rgb(0,1,1,.6),border=NA)
	}
legend(1.5,4,c("20% pool","10% pool"),fill=c(rgb(0,0,0,.4),col=rgb(0,1,1,.6)),bty="n",border=NA)
segments(0.5,q025(cpe100[500,]),9.5,q025(cpe100[500,]),lty=2)
segments(0.5,q975(cpe100[500,]),9.5,q975(cpe100[500,]),lty=2)
axis(1,at=c(-1,1:7,9),labels=c("",c(7.5,10,15,30,60,90,120),""),pos=0)
axis(2,las=T,pos=0.5)
text(-0.95,4.45,"D)",xpd=T)
#

# code to plot FigS4
par(mfrow=c(3,3))
par(mar=c(3,4,2,1))
plot(c(2002:2018)-0.1,apply(phi0_nat[,1,],2,median),pch=19,col=rgb(0,0,1,.4),ylim=c(0,1),axes=FALSE,xlim=c(2001.5,2018.5),main="San Acacia",ylab="",xlab="")
points(c(2002:2018)+.1,apply(phi_wdry0[,,1],2,median),pch=19,col=rgb(1,0,0,.4))
segments(c(2002:2018)+.1,apply(phi_wdry0[,,1],2,q025),c(2002:2018)+.1,apply(phi_wdry0[,,1],2,q975),col=rgb(1,0,0,.4))
segments(c(2002:2018)-.1,apply(phi0_nat[,1,],2,q025),c(2002:2018)-.1,apply(phi0_nat[,1,],2,q975),col=rgb(0,0,1,.4))
axis(1,seq(2002,2018,4),pos=0)
axis(1,c(2001.5,2018.15),c("",""),tck=0,pos=0)
axis(2,c(0,0.25,0.5,0.75,1),las=T,pos=2001.5)
mtext("Age-0 Survival (Jul-1  to Nov-1)",2,2.25,cex=0.7)
legend(2001.5,1.03,c("Survival without drying"),pch=19,col=rgb(0,0,1,.4),bty="n")
legend(2010,1.03,c("Survival with drying"),pch=19,col=rgb(1,0,0,.4),bty="n")
#
plot(c(2002:2018)-0.1,apply(phi0_nat[,2,],2,median),pch=19,col=rgb(0,0,1,.4),ylim=c(0,1),axes=FALSE,xlim=c(2001.5,2018.5),main="Isleta",ylab="",xlab="")
points(c(2002:2018)+.1,apply(phi_wdry0[,,2],2,median),pch=19,col=rgb(1,0,0,.4))
segments(c(2002:2018)+.1,apply(phi_wdry0[,,2],2,q025),c(2002:2018)+.1,apply(phi_wdry0[,,2],2,q975),col=rgb(1,0,0,.4))
segments(c(2002:2018)-.1,apply(phi0_nat[,2,],2,q025),c(2002:2018)-.1,apply(phi0_nat[,2,],2,q975),col=rgb(0,0,1,.4))
axis(1,seq(2002,2018,4),pos=0)
axis(1,c(2001.5,2018.15),c("",""),tck=0,pos=0)
axis(2,c(0,0.25,0.5,0.75,1),las=T,pos=2001.5)
mtext("Age-0 Survival (Jul-1  to Nov-1)",2,2.25,cex=0.7)
legend(2001.5,1.03,c("Survival without drying"),pch=19,col=rgb(0,0,1,.4),bty="n")
legend(2010,1.03,c("Survival with drying"),pch=19,col=rgb(1,0,0,.4),bty="n")
#
plot(c(2002:2018)-0.1,apply(phi0_nat[,3,],2,median),pch=19,col=rgb(0,0,1,.4),ylim=c(0,1),axes=FALSE,xlim=c(2001.5,2018.5),main="Angostura",ylab="",xlab="")
points(c(2002:2018)+.1,apply(phi_wdry0[,,3],2,median),pch=19,col=rgb(1,0,0,.4))
segments(c(2002:2018)+.1,apply(phi_wdry0[,,3],2,q025),c(2002:2018)+.1,apply(phi_wdry0[,,3],2,q975),col=rgb(1,0,0,.4))
segments(c(2002:2018)-.1,apply(phi0_nat[,3,],2,q025),c(2002:2018)-.1,apply(phi0_nat[,3,],2,q975),col=rgb(0,0,1,.4))
axis(1,seq(2002,2018,4),pos=0)
axis(1,c(2001.5,2018.15),c("",""),tck=0,pos=0)
axis(2,c(0,0.25,0.5,0.75,1),las=T,pos=2001.5)
mtext("Age-0 Survival (Jul-1  to Nov-1)",2,2.25,cex=0.7)
legend(2001.5,1.03,c("Survival without drying"),pch=19,col=rgb(0,0,1,.4),bty="n")
legend(2010,1.03,c("Survival with drying"),pch=19,col=rgb(1,0,0,.4),bty="n")
#
plot(c(2002:2018)-0.1,apply(phi1_nat[,1,],2,median),pch=19,col=rgb(0,0,1,.4),ylim=c(0,1),axes=FALSE,xlim=c(2001.5,2018.5),main="",ylab="",xlab="")
points(c(2002:2018)+.1,apply(phi_wdry1[,,1],2,median),pch=19,col=rgb(1,0,0,.4))
segments(c(2002:2018)+.1,apply(phi_wdry1[,,1],2,q025),c(2002:2018)+.1,apply(phi_wdry1[,,1],2,q975),col=rgb(1,0,0,.4))
segments(c(2002:2018)-.1,apply(phi1_nat[,1,],2,q025),c(2002:2018)-.1,apply(phi1_nat[,1,],2,q975),col=rgb(0,0,1,.4))
axis(1,seq(2002,2018,4),pos=0)
axis(1,c(2001.5,2018.15),c("",""),tck=0,pos=0)
axis(2,c(0,0.25,0.5,0.75,1),las=T,pos=2001.5)
mtext("Age-1+ Survival (April-1  to Nov-1)",2,2.25,cex=0.7)
legend(2001.5,1.03,c("Survival without drying"),pch=19,col=rgb(0,0,1,.4),bty="n")
legend(2010,1.03,c("Survival with drying"),pch=19,col=rgb(1,0,0,.4),bty="n")
#
plot(c(2002:2018)-0.1,apply(phi1_nat[,2,],2,median),pch=19,col=rgb(0,0,1,.4),ylim=c(0,1),axes=FALSE,xlim=c(2001.5,2018.5),main="",ylab="",xlab="")
points(c(2002:2018)+.1,apply(phi_wdry1[,,2],2,median),pch=19,col=rgb(1,0,0,.4))
segments(c(2002:2018)+.1,apply(phi_wdry1[,,2],2,q025),c(2002:2018)+.1,apply(phi_wdry1[,,2],2,q975),col=rgb(1,0,0,.4))
segments(c(2002:2018)-.1,apply(phi1_nat[,2,],2,q025),c(2002:2018)-.1,apply(phi1_nat[,2,],2,q975),col=rgb(0,0,1,.4))
axis(1,seq(2002,2018,4),pos=0)
axis(1,c(2001.5,2018.15),c("",""),tck=0,pos=0)
axis(2,c(0,0.25,0.5,0.75,1),las=T,pos=2001.5)
mtext("Age-1+ Survival (April-1  to Nov-1)",2,2.25,cex=0.7)
legend(2001.5,1.03,c("Survival without drying"),pch=19,col=rgb(0,0,1,.4),bty="n")
legend(2010,1.03,c("Survival with drying"),pch=19,col=rgb(1,0,0,.4),bty="n")
#
plot(c(2002:2018)-0.1,apply(phi1_nat[,3,],2,median),pch=19,col=rgb(0,0,1,.4),ylim=c(0,1),axes=FALSE,xlim=c(2001.5,2018.5),main="",ylab="",xlab="")
points(c(2002:2018)+.1,apply(phi_wdry1[,,3],2,median),pch=19,col=rgb(1,0,0,.4))
segments(c(2002:2018)+.1,apply(phi_wdry1[,,3],2,q025),c(2002:2018)+.1,apply(phi_wdry1[,,3],2,q975),col=rgb(1,0,0,.4))
segments(c(2002:2018)-.1,apply(phi1_nat[,3,],2,q025),c(2002:2018)-.1,apply(phi1_nat[,3,],2,q975),col=rgb(0,0,1,.4))
axis(1,seq(2002,2018,4),pos=0)
axis(1,c(2001.5,2018.15),c("",""),tck=0,pos=0)
axis(2,c(0,0.25,0.5,0.75,1),las=T,pos=2001.5)
mtext("Age-1+ Survival (April-1  to Nov-1)",2,2.25,cex=0.7)
legend(2001.5,1.03,c("Survival without drying"),pch=19,col=rgb(0,0,1,.4),bty="n")
legend(2010,1.03,c("Survival with drying"),pch=19,col=rgb(1,0,0,.4),bty="n")
#Over winter
plot(c(2002:2018)-0.1,apply(phiW[,1,1:17],2,median),pch=19,col=rgb(0,0,1,.4),ylim=c(0,1),axes=FALSE,xlim=c(2001.5,2018.5),main="",ylab="",xlab="")
segments(c(2002:2018)-.1,apply(phiW[,1,1:17],2,q025),c(2002:2018)-.1,apply(phiW[,1,1:17],2,q975),col=rgb(0,0,1,.4))
axis(1,seq(2002,2018,4),pos=0)
axis(1,c(2001.5,2018.15),c("",""),tck=0,pos=0)
axis(2,c(0,0.25,0.5,0.75,1),las=T,pos=2001.5)
mtext("Overwinter Survival (Nov-1 to April-1)",2,2.25,cex=0.7)
#
plot(c(2002:2018)-0.1,apply(phiW[,2,1:17],2,median),pch=19,col=rgb(0,0,1,.4),ylim=c(0,1),axes=FALSE,xlim=c(2001.5,2018.5),main="",ylab="",xlab="")
segments(c(2002:2018)-.1,apply(phiW[,2,1:17],2,q025),c(2002:2018)-.1,apply(phiW[,2,1:17],2,q975),col=rgb(0,0,1,.4))
axis(1,seq(2002,2018,4),pos=0)
axis(1,c(2001.5,2018.15),c("",""),tck=0,pos=0)
axis(2,c(0,0.25,0.5,0.75,1),las=T,pos=2001.5)
mtext("Overwinter Survival (Nov-1 to April-1)",2,2.25,cex=0.7)
#
plot(c(2002:2018)-0.1,apply(phiW[,3,1:17],2,median),pch=19,col=rgb(0,0,1,.4),ylim=c(0,1),axes=FALSE,xlim=c(2001.5,2018.5),main="",ylab="",xlab="")
segments(c(2002:2018)-.1,apply(phiW[,3,1:17],2,q025),c(2002:2018)-.1,apply(phiW[,3,1:17],2,q975),col=rgb(0,0,1,.4))
axis(1,seq(2002,2018,4),pos=0)
axis(1,c(2001.5,2018.15),c("",""),tck=0,pos=0)
axis(2,c(0,0.25,0.5,0.75,1),las=T,pos=2001.5)
mtext("Overwinter Survival (Nov-1 to April-1)",2,2.25,cex=0.7)

# do calculations for figure 6
# load back in t1,t2,t3,t4 used to calcuate preds_out
temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
t1<-numeric()
for (t in 1:6){
	se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
	t1[t]<-temp[t,6]-se
	se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
	t1[(t+7)]<-temp[(t+6),6]+se
	}
t1[7]<-1
t1[14]<-0
temp<-subset(ee[[3]],ee[[3]][,1]=="2")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t2<-temp[1,6]
temp<-subset(ee[[3]],ee[[3]][,1]=="4")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t4<-temp[1,6]-se
# use r<-1 san acacia
r<-1
	temp<-subset(ee[[3]],ee[[3]][,1]==riversegment[r]&is.na(ee[[3]][,5])==FALSE)
	t3<-numeric()
	for (t in 1:6){
		se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
		t3[t]<-temp[t,6]+se
		}	
	for (t in 8:15){
		se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
		t3[t]<-temp[t,6]-se
		}
	t3[7]<-temp[7,6]
#
EFFS<-extract(M2_1re,"effS")[[1]]
B_LBETA<-extract(M2_1re,"B_lbeta")[[1]]
M0<-extract(M2_1re,"M0")[[1]]
MOVE<-extract(M2_1re,"move")[[1]]
RP0<-extract(M2_1re,"rp0")[[1]]
A<-extract(M2_1re,"a")[[1]]
n<-array(NA,dim=c(17,5,15000))
eps<-extract(M2_1re,"lbeta_eps")[[1]]
efs<-10^seq(2,7,.1)
hab<-seq(0,1.25,0.01)
mn_a<-summary(M2_1re,"a")$summary[,1]
octNo<-array(NA,dim=c(length(efs),length(hab)))
meanphi0<-mean(phi_wdry0[,,1])
B<-summary(M2_1re,"B_lbeta")$summary[,1]
mu_lb<-summary(M2_1re,"mu_lbeta")$summary[1,1]+.5*(summary(M2_1re,"sd_lbeta")$summary[,1])^2
for (i in 1:length(hab)){
	for (j in 1:length(efs)){
			octNo[j,i]<-mn_a*meanphi0*(efs[j])/(1+mn_a*(efs[j])/(exp(mu_lb+B*(hab[i]-mean(preds2)))))
		}}
colz<-rep(rgb(0,0,0,.3),19)
colz[16]<-"blue3"
colz[18]<-"forestgreen"
# plot figure 6
par(mfrow=c(2,2))
par(mar=c(3,4,1,1))
plot(1000*c(simpflow[,1],simpflow_out[,1]),c(preds2[,1],preds_out[,1]),pch=c(rep(19,17),8,8),col=colz,axes=FALSE,ylab="",xlab="",xlim=c(0,5000),ylim=c(0,1.3))
axis(1,pos=0)
axis(2,seq(0,1.25,0.25),las=T,pos=0)
mtext("Mean May and June discharge (cfs)",1,2,cex=0.85)
mtext("Larval CC index",2,2.5,cex=0.85)
text(-1300,1.35,"A)",xpd=T)
#
contour(x=log10(efs),y=hab,octNo[,],levels=10^c(3:7),axes=FALSE,main="Age-0 production",cex.main=1)
axis(1,at=seq(1,7),c("","100","1000","10^4","10^5","10^6",""),pos=0)
axis(2,seq(0,1.25,0.25),las=T,pos=2)	
mtext("Number of effective spawners",1,2,cex=0.85)
mtext("Larval CC index",2,2.5,cex=0.85)
points(log10(c(summary(M2_1re,"effS")$summary[seq(1,49,3),6],apply(f1$effS[,,1],2,median))),c(preds2[,1],preds_out[,1]),pch=c(rep(19,17),8,8),col=colz)
text(0.7,1.27,"B)",xpd=T)
##
plot(c(32:155),sQ[32:155,16],type="l",axes=FALSE,ylab="",xlab="",ylim=c(0,6000),xlim=c(32,155),col=colz[16])
points(c(32:155),sQ_out[32:155,1],type="l",col=colz[18])
axis(2,c(0,3000,6000),las=T,pos=32)
axis(1,c(-10,32,62,93,124,155),c("","Apr-1","May-1","Jun-1","Jul-1","Aug-1"),pos=0)
mtext("Mean daily flow (cfs)",2,2.5,cex=0.85)
mtext("Month-Day",1,2,cex=0.85)
text(1,6600,"C)",xpd=T)
# 
curve(mn_a*meanphi0*(10^x)/(1+mn_a*(10^x)/(exp(mu_lb+B*(preds2[16,1]-mean(preds2))))),from=3,to=6,axes=FALSE,ylab="",xlab="",log="xy",ylim=c(1000000,8000000),col=colz[16])
curve(mn_a*meanphi0*(10^x)/(1+mn_a*(10^x)/(exp(mu_lb+B*(preds_out[1,1]-mean(preds2))))),add=T,col=colz[18])
points(log10(summary(M2_1re,"effS")$summary[46,6]),mn_a*meanphi0*(summary(M2_1re,"effS")$summary[46,6])/(1+mn_a*(summary(M2_1re,"effS")$summary[46,6])/(exp(mu_lb+B*(preds2[16,1]-mean(preds2))))),col=colz[16],pch=19)
points(log10(median(f1$effS[,1,1])),mn_a*meanphi0*(median(f1$effS[,1,1]))/(1+mn_a*(median(f1$effS[,1,1]))/(exp(mu_lb+B*(preds_out[1,1]-mean(preds2))))),col=colz[18],pch=19)
axis(2,las=T,at=c(1000000,2000000,4000000,8000000),c(1,2,4,8),pos=3)
axis(1,at=c(2:6),c("","10^3","10^4","10^5","10^6"),pos=1000000)
mtext("Number of effective spawners",1,2,cex=0.85)
mtext("Age-0 production",2,2.5,cex=0.85)
mtext("(x 10^6)",2,1.5,cex=0.85)
text(2.5,9000000,"D)",xpd=T)
#


######################
#### STEP 6: Evaluate management effectiveness
######################	


#algorithm to calculate near optimal Q by optimizing for san acacia
buildQ<-function(steps=25){
	stepQ<-cbind(diag(flexdays)*totflexQ/steps,matrix(0,nrow=flexdays,ncol=length(fixQ)-flexdays))
	for (s in 2:(2*ceiling(t4))){
		tstepQ<-matrix(0,nrow=(flexdays-s+1),ncol=length(fixQ))
		for (i in 1:(flexdays-s+1)){
			for (j in i:(i+s-1)){
				tstepQ[i,j]<-totflexQ/steps/s}}
		stepQ<-rbind(stepQ,tstepQ)
		}
	Noptions<-length(stepQ[,1])
	setQ<-fixQ
	temp<-numeric()
	for (j in 1:steps){
		temp<-numeric()
		for (i in 1:Noptions){
			temp[i]<-calc_cov(stepQ[i,]+setQ,t3,t2,t4,t1)
			}
		setQ<-stepQ[which(temp==max(temp))[1],]+setQ}
	rfQ<-calc_cov(setQ,t3,t2,t4,t1)
	trfQ<-rfQ+.0000001
	t<-0
	tsetQ<-setQ
	while(trfQ>rfQ&t<100){
		setQ<-tsetQ
		rfQ<-trfQ
		t<-t+1
		temp<-numeric()
		for (i in 1:Noptions){
			if (length(which(setQ-stepQ[i,]<fixQ))==0){
			temp[i]<-calc_cov(setQ-stepQ[i,],t3,t2,t4,t1)} else {temp[i]<-0}
			}
		tsetQ<-setQ-stepQ[which(temp==max(temp))[1],]
		temp<-numeric()
		for (i in 1:Noptions){
			temp[i]<-calc_cov(stepQ[i,]+tsetQ,t3,t2,t4,t1)
			}
		tsetQ<-tsetQ+stepQ[which(temp==max(temp))[1],]
		trfQ<-calc_cov(tsetQ,t3,t2,t4,t1)
		}
	trfQ<-rfQ+.0000001
	t<-0
	while(trfQ>rfQ&t<100){
		setQ<-tsetQ
		rfQ<-trfQ
		t<-t+1
		temp<-numeric()
		for (i in 1:Noptions){
			if (length(which(setQ-stepQ[i,]/2<fixQ))==0){
			temp[i]<-calc_cov(setQ-stepQ[i,]/2,t3,t2,t4,t1)} else {temp[i]<-0}
			}
		tsetQ<-setQ-stepQ[which(temp==max(temp))[1],]/2
		temp<-numeric()
		for (i in 1:Noptions){
			temp[i]<-calc_cov(stepQ[i,]/2+tsetQ,t3,t2,t4,t1)
			}
		tsetQ<-tsetQ+stepQ[which(temp==max(temp))[1],]/2
		trfQ<-calc_cov(tsetQ,t3,t2,t4,t1)
		}
	trfQ<-rfQ+.0000001
	t<-0
	while(trfQ>rfQ&t<100){
		setQ<-tsetQ
		rfQ<-trfQ
		t<-t+1
		temp<-numeric()
		for (i in 1:Noptions){
			if (length(which(setQ-stepQ[i,]/4<fixQ))==0){
			temp[i]<-calc_cov(setQ-stepQ[i,]/4,t3,t2,t4,t1)} else {temp[i]<-0}
			}
		tsetQ<-setQ-stepQ[which(temp==max(temp))[1],]/4
		temp<-numeric()
		for (i in 1:Noptions){
			temp[i]<-calc_cov(stepQ[i,]/4+tsetQ,t3,t2,t4,t1)
			}
		tsetQ<-tsetQ+stepQ[which(temp==max(temp))[1],]/4
		trfQ<-calc_cov(tsetQ,t3,t2,t4,t1)
		}
	trfQ<-rfQ+.0000001
	t<-0
	while(trfQ>rfQ&t<100){
		setQ<-tsetQ
		rfQ<-trfQ
		t<-t+1
		temp<-numeric()
		for (i in 1:Noptions){
			if (length(which(setQ-stepQ[i,]/10<fixQ))==0){
			temp[i]<-calc_cov(setQ-stepQ[i,]/10,t3,t2,t4,t1)} else {temp[i]<-0}
			}
		tsetQ<-setQ-stepQ[which(temp==max(temp))[1],]/10
		temp<-numeric()
		for (i in 1:Noptions){
			temp[i]<-calc_cov(stepQ[i,]/10+tsetQ,t3,t2,t4,t1)
			}
		tsetQ<-tsetQ+stepQ[which(temp==max(temp))[1],]/10
		trfQ<-calc_cov(tsetQ,t3,t2,t4,t1)
		}
	return(list(setQ=setQ,tsetQ=tsetQ,t=t,cv=calc_cov(setQ,t3,t2,t4,t1)))}
# function assumes fixed Q vector, totflexQ, and flexdays are in workspace
add10_array<-array(NA,dim=c(183,17,3))
add10_cv<-matrix(NA,nrow=17,ncol=3)
add30_array<-array(NA,dim=c(183,17,3))
add30_cv<-matrix(NA,nrow=17,ncol=3)
## this loop can take many hours - BEWARE
for (i in 1:17){
	fixQ<-sQ[1:183,i]
	flexdays<-154
	totflexQ<-5042 #10000 acre feet
	temp<-buildQ()
	temp2<-buildQ(steps=10)
	temp3<-buildQ(steps=3)
	add10_array[,i,1]<-temp$setQ
	add10_array[,i,2]<-temp2$setQ
	add10_array[,i,3]<-temp3$setQ
	add10_cv[i,]<-c(temp$cv,temp2$cv,temp3$cv)
	totflexQ<-5042*3 #30000 acre feet
	temp<-buildQ()
	temp2<-buildQ(steps=10)
	temp3<-buildQ(steps=3)
	add30_array[,i,1]<-temp$setQ
	add30_array[,i,2]<-temp2$setQ
	add30_array[,i,3]<-temp3$setQ
	add30_cv[i,]<-c(temp$cv,temp2$cv,temp3$cv)
	}
wmax<-function(x){which(x==max(x))[1]}
add10<-matrix(NA,ncol=17,nrow=183)
add30<-matrix(NA,ncol=17,nrow=183)
for (i in 1:17){
	add10[,i]<-add10_array[,i,wmax(add10_cv[i,])]
	add30[,i]<-add30_array[,i,wmax(add30_cv[i,])]
		}

#now calculate resulting abundance
base_cov<-matrix(NA,ncol=3,nrow=17)
add10_cov<-matrix(NA,ncol=3,nrow=17)
add30_cov<-matrix(NA,ncol=3,nrow=17)
add10_mj<-matrix(NA,ncol=3,nrow=17)
add30_mj<-matrix(NA,ncol=3,nrow=17)
#
t3mat<-matrix(NA,nrow=3,ncol=15)
for (r in 1:3){
	temp<-subset(ee[[3]],ee[[3]][,1]==riversegment[r]&is.na(ee[[3]][,5])==FALSE)
	for (t in 1:6){
		se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
		t3mat[r,t]<-temp[t,6]+se
		}	
	for (t in 8:15){
		se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
		t3mat[r,t]<-temp[t,6]-se
		}
	t3mat[r,7]<-temp[7,6]
	}
#add just to may 16 to june 9 - from Walsworth paper
mj<-rep(0,183)
mj[77:101]<-1/25

for (i in 1:17){
	base_cov[i,]<-c(calc_cov(sQ[1:183,i],t3mat[1,],t2,t4,t1)-mean(preds2),
		calc_cov(sQ[1:183,i],t3mat[2,],t2,t4,t1)-mean(preds2),
		calc_cov(aQ[1:183,i],t3mat[3,],t2,t4,t1)-mean(preds2)) 
	add10_cov[i,]<-c(calc_cov(add10[1:183,i],t3mat[1,],t2,t4,t1)-mean(preds2),
		calc_cov(add10[1:183,i],t3mat[2,],t2,t4,t1)-mean(preds2),
		calc_cov(aQ[1:183,i]+add10[1:183,i]-sQ[1:183,i],t3mat[3,],t2,t4,t1)-mean(preds2)) 
	add30_cov[i,]<-c(calc_cov(add30[1:183,i],t3mat[1,],t2,t4,t1)-mean(preds2),
		calc_cov(add30[1:183,i],t3mat[2,],t2,t4,t1)-mean(preds2),
		calc_cov(aQ[1:183,i]+add30[1:183,i]-sQ[1:183,i],t3mat[3,],t2,t4,t1)-mean(preds2))
	add10_mj[i,]<-c(calc_cov(sQ[1:183,i]+5042*mj,t3mat[1,],t2,t4,t1)-mean(preds2),
		calc_cov(sQ[1:183,i]+5042*mj,t3mat[2,],t2,t4,t1)-mean(preds2),
		calc_cov(aQ[1:183,i]+5042*mj,t3mat[3,],t2,t4,t1)-mean(preds2)) 
	add30_mj[i,]<-c(calc_cov(sQ[1:183,i]+5042*mj*3,t3mat[1,],t2,t4,t1)-mean(preds2),
		calc_cov(sQ[1:183,i]+5042*mj*3,t3mat[2,],t2,t4,t1)-mean(preds2),
		calc_cov(aQ[1:183,i]+5042*mj*3,t3mat[3,],t2,t4,t1)-mean(preds2)) 
	#
	tn1<-A*EFFS[,i,1]/(1+A*EFFS[,i,1]/exp(eps[,1,i]+B_LBETA*base_cov[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1])+(1-MOVE)*RP0*(cum_phiR[i,215,1]-cum_phiR[i,91,1]))
	tn2<-A*EFFS[,i,2]/(1+A*EFFS[,i,2]/exp(eps[,2,i]+B_LBETA*base_cov[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2])+(1-MOVE)*RP0*(cum_phiR[i,215,2]-cum_phiR[i,91,2]))
	tn3<-A*EFFS[,i,3]/(1+A*EFFS[,i,3]/exp(eps[,3,i]+B_LBETA*base_cov[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3])+(1-MOVE)*RP0*(cum_phiR[i,215,3]-cum_phiR[i,91,3]))
	n[i,1,]<-tn1+tn2+tn3
	#
	tn1<-A*EFFS[,i,1]/(1+A*EFFS[,i,1]/exp(eps[,1,i]+B_LBETA*add10_cov[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1])+(1-MOVE)*RP0*(cum_phiR[i,215,1]-cum_phiR[i,91,1]))
	tn2<-A*EFFS[,i,2]/(1+A*EFFS[,i,2]/exp(eps[,2,i]+B_LBETA*add10_cov[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2])+(1-MOVE)*RP0*(cum_phiR[i,215,2]-cum_phiR[i,91,2]))
	tn3<-A*EFFS[,i,3]/(1+A*EFFS[,i,3]/exp(eps[,3,i]+B_LBETA*add10_cov[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3])+(1-MOVE)*RP0*(cum_phiR[i,215,3]-cum_phiR[i,91,3]))
	n[i,2,]<-tn1+tn2+tn3
	#
	tn1<-A*EFFS[,i,1]/(1+A*EFFS[,i,1]/exp(eps[,1,i]+B_LBETA*add30_cov[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1])+(1-MOVE)*RP0*(cum_phiR[i,215,1]-cum_phiR[i,91,1]))
	tn2<-A*EFFS[,i,2]/(1+A*EFFS[,i,2]/exp(eps[,2,i]+B_LBETA*add30_cov[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2])+(1-MOVE)*RP0*(cum_phiR[i,215,2]-cum_phiR[i,91,2]))
	tn3<-A*EFFS[,i,3]/(1+A*EFFS[,i,3]/exp(eps[,3,i]+B_LBETA*add30_cov[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3])+(1-MOVE)*RP0*(cum_phiR[i,215,3]-cum_phiR[i,91,3]))
	n[i,3,]<-tn1+tn2+tn3
	#
	tn1<-A*EFFS[,i,1]/(1+A*EFFS[,i,1]/exp(eps[,1,i]+B_LBETA*add10_mj[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1])+(1-MOVE)*RP0*(cum_phiR[i,215,1]-cum_phiR[i,91,1]))
	tn2<-A*EFFS[,i,2]/(1+A*EFFS[,i,2]/exp(eps[,2,i]+B_LBETA*add10_mj[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2])+(1-MOVE)*RP0*(cum_phiR[i,215,2]-cum_phiR[i,91,2]))
	tn3<-A*EFFS[,i,3]/(1+A*EFFS[,i,3]/exp(eps[,3,i]+B_LBETA*add10_mj[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3])+(1-MOVE)*RP0*(cum_phiR[i,215,3]-cum_phiR[i,91,3]))
	n[i,4,]<-tn1+tn2+tn3
	#
	tn1<-A*EFFS[,i,1]/(1+A*EFFS[,i,1]/exp(eps[,1,i]+B_LBETA*add30_mj[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1])+(1-MOVE)*RP0*(cum_phiR[i,215,1]-cum_phiR[i,91,1]))
	tn2<-A*EFFS[,i,2]/(1+A*EFFS[,i,2]/exp(eps[,2,i]+B_LBETA*add30_mj[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2])+(1-MOVE)*RP0*(cum_phiR[i,215,2]-cum_phiR[i,91,2]))
	tn3<-A*EFFS[,i,3]/(1+A*EFFS[,i,3]/exp(eps[,3,i]+B_LBETA*add30_mj[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3])+(1-MOVE)*RP0*(cum_phiR[i,215,3]-cum_phiR[i,91,3]))
	n[i,5,]<-tn1+tn2+tn3
	}

# plot figure 7
CX<-0.7
layout(matrix(c(1,1,2,3,4,5,6,6,7,7),byrow=T,ncol=2),heights=c(1,4,4,4,4))
par(mar=c(0,0,0,0))
plot(0,xlim=c(0,1),ylim=c(0,1),axes=FALSE,ylab="",xlab="",col="white")
legend(0.08,1,c("10 kaf added May 16th to June 9","10 kaf added near optimally"),pch=19,col=c(rgb(0,0,1,.4),"darkblue"),bty="n",cex=CX/0.7)
legend(0.62,1,c("30 kaf added May 16th to June 9","30 kaf added near optimally"),pch=19,col=c(rgb(0,1,0,.4),"darkgreen"),bty="n",cex=CX/0.7)
#
par(mar=c(3,4,0.5,0.5))
plot(preds2[,1],add10_cov[,1]+mean(preds2),pch=19,col="darkblue",xlim=c(0,1.25),ylim=c(0,1.25),axes=FALSE,ylab="",xlab="")
points(preds2[,1],add10_mj[,1]+mean(preds2),pch=19,col=rgb(0,0,1,.4))
curve(1*x,add=T,lty=2)
mtext("Larval CC index without added water",1,1.25,cex=CX)
mtext("Larval CC index",2,2,cex=CX)
axis(2,c(0,.5,1,1.5),las=T,pos=0,cex.axis=CX/0.7,hadj=0.5,tck=-0.03)
axis(1,c(0,.5,1,1.5),las=T,pos=0,cex.axis=CX/0.7,padj=-0.5,tck=-0.03)
text(-0.24,1.31,"A)",xpd=T)
#
plot(preds2[,1],add30_cov[,1]+mean(preds2),pch=19,col="forestgreen",xlim=c(0,1.25),ylim=c(0,1.25),axes=FALSE,ylab="",xlab="")
points(preds2[,1],add30_mj[,1]+mean(preds2),pch=19,col=rgb(0,1,0,.4))
curve(1*x,add=T,lty=2)
mtext("Larval CC index without added water",1,1.25,cex=CX)
mtext("Larval CC index",2,2,cex=CX)
axis(2,c(0,.5,1,1.5),las=T,pos=0,cex.axis=CX/0.7,hadj=0.5,tck=-0.03)
axis(1,c(0,.5,1,1.5),las=T,pos=0,cex.axis=CX/0.7,padj=-0.5,tck=-0.03)
text(-0.24,1.31,"B)",xpd=T)
#
plot(apply(n[,1,],1,median),apply(n[,2,],1,median),pch=19,col="darkblue",xlim=c(1000,10000000),log="xy",ylim=c(1000,10000000),axes=FALSE,ylab="",xlab="")
points(apply(n[,1,],1,median),apply(n[,4,],1,median),pch=19,col=rgb(0,0,1,.4))
curve(1*x,add=T,lty=2)
mtext("Age-0 abundance without added water",1,1.5,cex=CX)
mtext("Age-0 abundance",2,2,cex=CX)
axis(2,at=c(1000,100000,10000000),labels=c("10^3","10^5","10^7"),las=T,pos=1000,cex.axis=CX/0.7,hadj=0.7,tck=-0.03)
axis(1,at=c(1000,100000,10000000),labels=c("10^3","10^5","10^7"),pos=1000,cex.axis=CX/0.7,padj=-0.5,tck=-0.03)
text(170,15000000,"C)",xpd=T)
#
plot(apply(n[,1,],1,median),apply(n[,3,],1,median),pch=19,col="forestgreen",xlim=c(1000,10000000),ylim=c(1000,10000000),log="xy",axes=FALSE,ylab="",xlab="")
points(apply(n[,1,],1,median),apply(n[,5,],1,median),pch=19,col=rgb(0,1,0,.4))
points(preds2[,1],apply(n[,1,],1,median),pch=19,col=rgb(1,0,0,.2))
curve(1*x,add=T,lty=2)
mtext("Age-0 abundance without added water",1,1.5,cex=CX)
mtext("Age-0 abundance",2,2,cex=CX)
axis(2,at=c(1000,100000,10000000),labels=c("10^3","10^5","10^7"),las=T,pos=1000,cex.axis=CX/0.7,hadj=0.7,tck=-0.03)
axis(1,at=c(1000,100000,10000000),labels=c("10^3","10^5","10^7"),pos=1000,cex.axis=CX/0.7,padj=-0.5,tck=-0.03)
text(170,15000000,"D)",xpd=T)
#
par(mar=c(3,3.5,0.5,0.5))
plot(sQ[1:155,14],type="l",axes=FALSE,ylab="",xlab="",ylim=c(0,4000),xlim=c(1,155))
polygon(c(77:101,101:77),c(sQ[77:101,14],sQ[101:77,14]+5042/25),col=rgb(0,0,1,.4),border=NA)
polygon(c(77:101,101:77),c(sQ[77:101,14]+5042/25,sQ[101:77,14]+3*5042/25),col=rgb(0,1,0,.4),border=NA)
polygon(c(1:155,155:1),c(sQ[1:155,14],add30[155:1,14]),col="forestgreen",border=NA)
polygon(c(1:155,155:1),c(sQ[1:155,14],add10[155:1,14]),col="darkblue",border=NA)
axis(2,c(0,2000,4000),las=T,pos=1,cex.axis=CX/0.7,hadj=0.7,tck=-0.03)
axis(1,c(1,32,62,93,124,155),c("Mar-1","Apr-1","May-1","Jun-1","Jul-1","Aug-1"),pos=0,cex.axis=CX/0.7,padj=-.5,tck=-0.03)
axis(1,c(0,160),c("",""),pos=0,cex.axis=CX/0.7,padj=-.5,tck=0)
mtext("Mean daily flow (cfs)",2,1.5,cex=CX)
mtext("Month-Day",1,1,cex=CX)
text(4,4000,"Larval CC index values",cex=CX/0.7,adj=0)
text(4,3500,"baseline (2015): 0.67",cex=CX/0.7,adj=0)
text(4,3000,"10 kaf (5/16-6/19): 0.71",cex=CX/0.7,adj=0,col=rgb(0,0,1,.4))
text(4,2500,"10 kaf (near optimal): 0.91",cex=CX/0.7,adj=0,col="darkblue")
text(4,2000,"30 kaf (5/16-6/19): 0.74",cex=CX/0.7,adj=0,col=rgb(0,1,0,.4))
text(4,1500,"30 kaf (near optimal): 0.95",cex=CX/0.7,adj=0,col="forestgreen")
text(-13,4200,"E)",xpd=T)
#
plot(sQ[1:155,9],type="l",axes=FALSE,ylab="",xlab="",ylim=c(0,4000),xlim=c(1,155))
polygon(c(77:101,101:77),c(sQ[77:101,9],sQ[101:77,9]+5042/25),col=rgb(0,0,1,.4),border=NA)
polygon(c(77:101,101:77),c(sQ[77:101,9]+5042/25,sQ[101:77,9]+3*5042/25),col=rgb(0,1,0,.4),border=NA)
polygon(c(1:155,155:1),c(sQ[1:155,9],add30[155:1,9]),col="forestgreen",border=NA)
polygon(c(1:155,155:1),c(sQ[1:155,9],add10[155:1,9]),col="darkblue",border=NA)
axis(2,c(0,2000,4000),las=T,pos=1,cex.axis=CX/0.7,hadj=0.7,tck=-0.03)
axis(1,c(1,32,62,93,124,155),c("Mar-1","Apr-1","May-1","Jun-1","Jul-1","Aug-1"),pos=0,cex.axis=CX/0.7,padj=-.5,tck=-0.03)
axis(1,c(0,160),c("",""),pos=0,cex.axis=CX/0.7,padj=-.5,tck=0)
mtext("Mean daily flow (cfs)",2,1.5,cex=CX)
mtext("Month-Day",1,1,cex=CX)
text(110,4000,"Larval CC index values",cex=CX/0.7,adj=0)
text(110,3500,"baseline (2010): 0.77",cex=CX/0.7,adj=0)
text(110,3000,"10 kaf (5/16-6/19): 0.77",cex=CX/0.7,adj=0,col=rgb(0,0,1,.4))
text(110,2500,"10 kaf (near optimal): 0.79",cex=CX/0.7,adj=0,col="darkblue")
text(110,2000,"30 kaf (5/16-6/19): 0.77",cex=CX/0.7,adj=0,col=rgb(0,1,0,.4))
text(110,1500,"30 kaf (near optimal): 0.81",cex=CX/0.7,adj=0,col="forestgreen")
text(-13,4200,"F)",xpd=T)

# estimated impacts of different management scenarios for figure 8
M1<-extract(M2_1re,"M1")[[1]]
RP1<-extract(M2_1re,"rp1")[[1]]
SP_N<-extract(M2_1re,"Sp_N")[[1]]
BETA_2<-extract(M2_1re,"beta_2")[[1]]
# restored 
t3r_array<-array(NA,dim=c(3,15,15000))
riversegment2<-c("7c","7b","7a")
res_unc<-rnorm(15000)
for (r in 1:3){
	temp2<-subset(ee[[3]],ee[[3]][,1]==riversegment2[r]&is.na(ee[[3]][,5])==FALSE&ee[[3]][,2]!=500)
	for (i in 1:15){
		se<-ifelse(temp2[i,4]==temp2[i,3],0,calcsig(temp2[i,4],temp2[i,3],temp2[i,6],temp2[i,5]/100))
		t3r_array[r,i,]<-t3mat[r,i]*(.1*(ifelse(temp2[i,6]+se*res_unc<1,1,temp2[i,6]+se*res_unc))+.9)
		}}
#10% less drying
lessdry<-function(x){
	if (sum(x)<=0.1) (rep(0,215)) else {
		tx<-x[215:1]
		wx<-which(cumsum(tx)>0.1)[1]
		tt2<-sum(tx[1:(wx-1)])
		tx[1:(wx-1)]<-0
		tx[wx]<-tx[wx]-(0.1-tt2)
		tx[215:1]}
		}
prop_nd_10less<-prop_nd
for (i in 1:17){
	for (j in 1:3){
		prop_nd_10less[i,,j]<-lessdry(prop_nd[i,,j])}}
#
cum_nd_10less<-cum_nd
cum_phiR_10less<-cum_phiR
for (t in 1:Nyears){
	cum_nd_10less[t,,1]<-cumsum(prop_nd_10less[t,,1])
	cum_nd_10less[t,,2]<-cumsum(prop_nd_10less[t,,2])
	cum_nd_10less[t,,3]<-cumsum(prop_nd_10less[t,,3])
	cum_phiR_10less[t,,1]<-cumsum(prop_nd_10less[t,,1]*pred_phiR)
	cum_phiR_10less[t,,2]<-cumsum(prop_nd_10less[t,,2]*pred_phiR)
	cum_phiR_10less[t,,3]<-cumsum(prop_nd_10less[t,,3]*pred_phiR)
	}
#Effective spawners without augmentation
EFFS2<-EFFS
for (i in 1:17){
	for (j in 1:3){
		EFFS2[,i,j]<-SP_N[,i,j,1]+BETA_2*SP_N[,i,j,2]}}
#
scenarios<-array(NA,dim=c(17,6,3,15000))
# baseline is no water added, no augmentation, no restoration, no salvage, and no restoration
# this loop takes a while, but not hours
for (i in 1:17){
	restore_cov<-matrix(NA,ncol=3,nrow=15000)
	for (j in 1:15000){
			restore_cov[j,]<-c(calc_cov(sQ[1:183,i],t3r_array[1,,j],t2,t4,t1)-mean(preds2),
		calc_cov(sQ[1:183,i],t3r_array[2,,j],t2,t4,t1)-mean(preds2),
		calc_cov(aQ[1:183,i],t3r_array[3,,j],t2,t4,t1)-mean(preds2))}
	#scenario 1 - base case: no augmentation, no added water, no restoration, no salvage, all drying
	tn1<-A*EFFS2[,i,1]/(1+A*EFFS2[,i,1]/exp(eps[,1,i]+B_LBETA*base_cov[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1]))
	tn2<-A*EFFS2[,i,2]/(1+A*EFFS2[,i,2]/exp(eps[,2,i]+B_LBETA*base_cov[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2]))
	tn3<-A*EFFS2[,i,3]/(1+A*EFFS2[,i,3]/exp(eps[,3,i]+B_LBETA*base_cov[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3]))
	tn4<-(SP_N[,i,1,1]+SP_N[,i,1,2])*exp(-M1[,1,i]*215)*(1+(MOVE-1)*cum_nd[i,215,1])
	tn5<-(SP_N[,i,2,1]+SP_N[,i,2,2])*exp(-M1[,2,i]*215)*(1+(MOVE-1)*cum_nd[i,215,2])
	tn6<-(SP_N[,i,3,1]+SP_N[,i,3,2])*exp(-M1[,3,i]*215)*(1+(MOVE-1)*cum_nd[i,215,3])
	scenarios[i,1,1,]<-tn1+tn4
	scenarios[i,1,2,]<-tn2+tn5
	scenarios[i,1,3,]<-tn3+tn6
	#scenario 2 - augmentation - a
	tn1<-A*EFFS[,i,1]/(1+A*EFFS[,i,1]/exp(eps[,1,i]+B_LBETA*base_cov[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1]))
	tn2<-A*EFFS[,i,2]/(1+A*EFFS[,i,2]/exp(eps[,2,i]+B_LBETA*base_cov[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2]))
	tn3<-A*EFFS[,i,3]/(1+A*EFFS[,i,3]/exp(eps[,3,i]+B_LBETA*base_cov[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3]))
	scenarios[i,2,1,]<-tn1+tn4
	scenarios[i,2,2,]<-tn2+tn5
	scenarios[i,2,3,]<-tn3+tn6
	# scenario 3 - restore 10% of habitat in each river segment only -r
	tn1<-A*EFFS2[,i,1]/(1+A*EFFS2[,i,1]/exp(eps[,1,i]+B_LBETA*restore_cov[,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1]))
	tn2<-A*EFFS2[,i,2]/(1+A*EFFS2[,i,2]/exp(eps[,2,i]+B_LBETA*restore_cov[,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2]))
	tn3<-A*EFFS2[,i,3]/(1+A*EFFS2[,i,3]/exp(eps[,3,i]+B_LBETA*restore_cov[,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3]))
	scenarios[i,3,1,]<-tn1+tn4
	scenarios[i,3,2,]<-tn2+tn5
	scenarios[i,3,3,]<-tn3+tn6
	#scenario 4 - add10kaf to base - f
	tn1<-A*EFFS2[,i,1]/(1+A*EFFS2[,i,1]/exp(eps[,1,i]+B_LBETA*add10_cov[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1]))
	tn2<-A*EFFS2[,i,2]/(1+A*EFFS2[,i,2]/exp(eps[,2,i]+B_LBETA*add10_cov[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2]))
	tn3<-A*EFFS2[,i,3]/(1+A*EFFS2[,i,3]/exp(eps[,3,i]+B_LBETA*add10_cov[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3]))
	scenarios[i,4,1,]<-tn1+tn4
	scenarios[i,4,2,]<-tn2+tn5
	scenarios[i,4,3,]<-tn3+tn6
    # scenario 5 - salvage only - s
	tn1<-A*EFFS2[,i,1]/(1+A*EFFS2[,i,1]/exp(eps[,1,i]+B_LBETA*base_cov[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd[i,215,1]-cum_nd[i,91,1])+(1-MOVE)*RP0*(cum_phiR[i,215,1]-cum_phiR[i,91,1]))
	tn2<-A*EFFS2[,i,2]/(1+A*EFFS2[,i,2]/exp(eps[,2,i]+B_LBETA*base_cov[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd[i,215,2]-cum_nd[i,91,2])+(1-MOVE)*RP0*(cum_phiR[i,215,2]-cum_phiR[i,91,2]))
	tn3<-A*EFFS2[,i,3]/(1+A*EFFS2[,i,3]/exp(eps[,3,i]+B_LBETA*base_cov[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd[i,215,3]-cum_nd[i,91,3])+(1-MOVE)*RP0*(cum_phiR[i,215,3]-cum_phiR[i,91,3]))
	tn4<-(SP_N[,i,1,1]+SP_N[,i,1,2])*exp(-M1[,1,i]*215)*(1+(MOVE-1)*cum_nd[i,215,1]+(1-MOVE)*RP1*cum_phiR[i,215,1])
	tn5<-(SP_N[,i,2,1]+SP_N[,i,2,2])*exp(-M1[,2,i]*215)*(1+(MOVE-1)*cum_nd[i,215,2]+(1-MOVE)*RP1*cum_phiR[i,215,2])
	tn6<-(SP_N[,i,3,1]+SP_N[,i,3,2])*exp(-M1[,3,i]*215)*(1+(MOVE-1)*cum_nd[i,215,3]+(1-MOVE)*RP1*cum_phiR[i,215,3])
	scenarios[i,5,1,]<-tn1+tn4
	scenarios[i,5,2,]<-tn2+tn5
	scenarios[i,5,3,]<-tn3+tn6
	#scenario 6 - decrease drying by 10% only - d
	tn1<-A*EFFS2[,i,1]/(1+A*EFFS2[,i,1]/exp(eps[,1,i]+B_LBETA*base_cov[i,1]))*exp(-M0[,1,i]*124)*
		((1-cum_nd[i,91,1])+(MOVE-1)*(cum_nd_10less[i,215,1]-cum_nd_10less[i,91,1]))
	tn2<-A*EFFS2[,i,2]/(1+A*EFFS2[,i,2]/exp(eps[,2,i]+B_LBETA*base_cov[i,2]))*exp(-M0[,2,i]*124)*
		((1-cum_nd[i,91,2])+(MOVE-1)*(cum_nd_10less[i,215,2]-cum_nd_10less[i,91,2]))
	tn3<-A*EFFS2[,i,3]/(1+A*EFFS2[,i,3]/exp(eps[,3,i]+B_LBETA*base_cov[i,3]))*exp(-M0[,3,i]*124)*
		((1-cum_nd[i,91,3])+(MOVE-1)*(cum_nd_10less[i,215,3]-cum_nd_10less[i,91,3]))
	tn4<-(SP_N[,i,1,1]+SP_N[,i,1,2])*exp(-M1[,1,i]*215)*(1+(MOVE-1)*cum_nd_10less[i,215,1])
	tn5<-(SP_N[,i,2,1]+SP_N[,i,2,2])*exp(-M1[,2,i]*215)*(1+(MOVE-1)*cum_nd_10less[i,215,2])
	tn6<-(SP_N[,i,3,1]+SP_N[,i,3,2])*exp(-M1[,3,i]*215)*(1+(MOVE-1)*cum_nd_10less[i,215,3])
	scenarios[i,6,1,]<-tn1+tn4
	scenarios[i,6,2,]<-tn2+tn5
	scenarios[i,6,3,]<-tn3+tn6
	}
per_benefit<-array(NA,dim=c(17,5,15000))
ben_mns<-numeric()
ben_up<-numeric()
ben_down<-numeric()
for (i in 1:17){
	for (j in 1:5){per_benefit[i,j,]<-apply(scenarios[i,(j+1),,],2,sum)/apply(scenarios[i,1,,],2,sum)}
	ben_mns[(i-1)*5+c(1:5)]<-100*(apply(per_benefit[i,,],1,median))
	ben_up[(i-1)*5+c(1:5)]<-100*(apply(per_benefit[i,,],1,q90))
	ben_down[(i-1)*5+c(1:5)]<-100*(apply(per_benefit[i,,],1,q10))
	}
ben_up<-ifelse(ben_mns>140,NA,ben_up)
ben_down<-ifelse(ben_mns>140,NA,ben_down)
ben_mns<-ifelse(ben_mns>140,NA,ben_mns)
xs<-rep(c(2002:2018),each=5)+rep(c(-2:2)/10,17)
#plot figure 8
colorz<-c("forestgreen","mediumorchid1","darkblue","red2","orange")
layout(matrix(c(1:3),ncol=1),heights=c(2,1,1))
par(mar=c(3,4,1,3))
plot(xs,ben_mns,pch=c(21,23,16,24,25),ylim=c(100,140),log="y",xlab="",ylab="",axes=FALSE,col=colorz,lwd=2,xlim=c(2001.5,2018.5))
segments(xs,ben_up,xs,ben_down,col=colorz)
axis(1,at=c(2001.4,2002:2018,2018.5),labels=c("",c(2002:2018),""),pos=100,tck=0)
axis(2,las=T,pos=2001.4)
mtext("Percent of baseline",2,2,cex=0.7)
arrows(xs[which(is.na(ben_mns)==T)],138,xs[which(is.na(ben_mns)==T)],140,lwd=2,col="darkblue",length=0.05)
text(2000.05,142,"A)",xpd=T)
legend(2010.4,130,c("Augmentation","Restoration","Flow","Flow: >140%","Fish Rescue","Reduce drying"),pch=c(21,23,16,16,24,25),col=c(colorz[1:3],"white",colorz[4:5]),bty="n")
arrows(2010.67,121.5,2010.67,123,lwd=2,col="darkblue",length=0.05)
#panel 2 - spring abundance going into year and amount augmented
aug<-c(rowSums(w_um),rowSums(w_m))[1:17]
sN<-extract(M2_1re,"Sp_N")[[1]]
N1<-extract(M2_1re,"N1")[[1]]
N0<-extract(M2_1re,"N0")[[1]]
n1<-N1[,-18,1]+N1[,-18,2]+N1[,-18,3]
sn<-sN[,,1,1]+sN[,,2,1]+sN[,,3,1]+sN[,,1,2]+sN[,,2,2]+sN[,,3,2]+sN[,,1,3]+sN[,,2,3]+sN[,,3,3]
plot(c(2002:2018),apply(sn,2,mean)[-18],pch=19,col="white",axes=FALSE,ylab="",xlab="",log="y",ylim=c(100,100000000),xlim=c(2001.5,2018.5))
segments(c(2002:2018),apply(sn[,-18]-n1,2,q975),c(2002:2018),apply(sn[,-18]-n1,2,q025),col=rgb(0,.5,1,.4),lwd=6)
segments(c(2002:2018),apply(n1,2,q975),c(2002:2018),apply(n1,2,q025),col="green",lwd=6)
axis(2,c(100,1000,10000,100000,1000000,10000000),c("10^2","10^3","10^4","10^5","10^6","10^7"),las=T,pos=2001.4)
axis(1,at=c(2001.4,2002:2018,2018.5),labels=c("",c(2002:2018),""),pos=100,tck=0)
mtext("Abundance",2,2,cex=0.7)
legend(2002.25,100000000,"Spring - wild born" ,fill="green",bty="n",border=NA)
legend(2012.25,100000000,"Spring - augmented" ,fill=rgb(0,.5,1,.4),bty="n",border=NA)
text(2000.05,199526231,"B)",xpd=T)
#panel 3 - flow during reproductive season + % drying
plot(c(2002:2018),preds2[,1],pch=19,col="blue",ylim=c(0,1.5),ylab="",xlab="",axes=FALSE,xlim=c(2001.5,2018.5))
rmdry<-(cum_nd[,215,1]*60.6+cum_nd[,215,2]*54)/154.8
points(c(2002:2018),rmdry*3,col="red",pch=19)
axis(4,c(0,.5,1,1.5),col="blue",las=T,pos=2018.4,col.axis="blue")
axis(2,col="red",at=c(0,0.6,1.2),labels=c(0,0.2,0.4),las=T,pos=2001.4,col.axis="red")
axis(1,at=c(2001.4,2002:2018,2018.5),labels=c("",c(2002:2018),""),pos=0,tck=0)
mtext("Larval carry capacity index",4,1,cex=0.7,col="blue")
mtext("Drying (%)",2,1.75,col="red",cex=0.7)
text(2000.05,1.55,"C)",xpd=T)

# now calculate by riversegment benefits for figure 9
byriversegment<-array(NA,dim=c(17,5,3))
for (i in 1:17){
	for (j in 1:5){
		for (k in 1:3){
			byriversegment[i,j,k]<-mean(scenarios[i,(j+1),k,]/scenarios[i,1,k,])
			}}}
byriversegment<-ifelse(byriversegment<1.05,1,byriversegment)
#
wmax<-function(x){which(x==max(x))[1]}
w2<-function(x){
	tx<-which(x==sort(x)[4])[1]
	ifelse(x[tx]<1.1,6,tx)}
w3<-function(x){
	tx<-which(x==sort(x)[3])[1]
	ifelse(x[tx]<1.1,6,tx)}
w4<-function(x){
	tx<-which(x==sort(x)[2])[1]
	ifelse(x[tx]<1.1,6,tx)}
# plot figure 9
layout(matrix(c(1:4),ncol=1),heights=c(1,2,2,2))
par(mar=c(0,0,0,0))
plot(0.1,0.5,xlim=c(0,1),ylim=c(0,1),axes=FALSE,col=grey(0.5),pch=16,cex=12,xlab="",ylab="")
#points(0.1,0.5,col=grey(0.55),pch=16,cex=12)
points(0.1,0.5,col=grey(0.7),pch=16,cex=9)
points(0.1,0.5,col=grey(0.9),pch=16,cex=5)
segments(c(0.1,0.142,0.15),c(0.58,0.5,0.34),c(0.25,0.25,0.25),c(0.75,0.53,0.3),lwd=1.5)
text(0.27,0.8,"Most effective",adj=0,cex=1.2)
text(0.27,0.55,"2nd most effective",adj=0,cex=1.2)
text(0.27,0.3,"3rd most effective",adj=0,cex=1.2)
#text(0.27,0.15,"4th most effective",adj=0,cex=1.2)
legend(0.57,1.1,c("Augment Spawners","Restore larval habitat","Increase flows by 10kaf","Fish rescue","Reduce drying"),col=colorz,pch=16,bty="n",cex=1.2,y.intersp=1.2)
#
par(mar=c(3.5,4,0.5,0.5))
plot(preds2[,1],apply(apply(EFFS2[,,1],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,1],1,w4)],pch=16,cex=5,log="y",ylab="",xlab="",axes=FALSE,xlim=c(0,1.25),ylim=c(10,2000000))
points(preds2[,1],apply(apply(EFFS2[,,1],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,1],1,w3)],pch=16,cex=4)
points(preds2[,1],apply(apply(EFFS2[,,1],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,1],1,w2)],pch=16,cex=3)
points(preds2[,1],apply(apply(EFFS2[,,1],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,1],1,wmax)],pch=16,cex=1.5)
axis(2,at=c(10,1000,1000000),labels=c("10","10^3","10^6"),las=T,pos=0)
axis(1,at=c(0,0.5,1,1.5),pos=10)
text(4.5,100,"San Acacia",cex=1.2)
mtext("Larval carrying capacity index",1,1.75,cex=0.7)
mtext("Effective Spawners",2,2,cex=0.7)
text(-.17,3000000,"A)",xpd=T)
#
plot(preds2[,2],apply(apply(EFFS2[,,2],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,2],1,w4)],pch=16,cex=5,log="y",ylab="",xlab="",axes=FALSE,xlim=c(0,1.25),ylim=c(10,2000000))
points(preds2[,2],apply(apply(EFFS2[,,2],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,2],1,w3)],pch=16,cex=4)
points(preds2[,2],apply(apply(EFFS2[,,2],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,2],1,w2)],pch=16,cex=3)
points(preds2[,2],apply(apply(EFFS2[,,2],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,2],1,wmax)],pch=16,cex=1.5)
axis(2,at=c(10,1000,1000000),labels=c("10","10^3","10^6"),las=T,pos=0)
axis(1,at=c(0,0.5,1,1.5),pos=10)
text(4.5,100,"Isleta",cex=1.2)
mtext("Larval carrying capacity index",1,1.75,cex=0.7)
mtext("Effective Spawners",2,2,cex=0.7)
text(-.17,3000000,"B)",xpd=T)
#
plot(preds2[,3],apply(apply(EFFS2[,,3],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,3],1,w4)],pch=16,cex=5,log="y",ylab="",xlab="",axes=FALSE,xlim=c(0,1.25),ylim=c(10,2000000))
points(preds2[,3],apply(apply(EFFS2[,,3],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,3],1,w3)],pch=16,cex=4)
points(preds2[,3],apply(apply(EFFS2[,,3],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,3],1,w2)],pch=16,cex=3)
points(preds2[,3],apply(apply(EFFS2[,,3],c(1,2),sum),2,mean),col=c(colorz,NA)[apply(byriversegment[,,3],1,wmax)],pch=16,cex=1.5)
axis(2,at=c(10,1000,1000000),labels=c("10","10^3","10^6"),las=T,pos=0)
axis(1,at=c(0,0.5,1,1.5),pos=10)
text(4.5,100,"Angostura",cex=1.2)
mtext("Larval carrying capacity index",1,1.75,cex=0.7)
mtext("Effective Spawners",2,2,cex=0.7)
text(-.17,3000000,"C)",xpd=T)

