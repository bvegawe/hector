##==============================================================================
## Script to calculate the normalizations for BRICK subcomponent inputs, 
## required because each model was designed to take in temperature relative 
## to some time period, and DAIS requires SLR relative to a time period. 
##
## Because we can't normalize to future values in Hector (it runs sequentially 
## in time), we approximate the offsets with the observed offsets.
##
## The values in the output csv are hard-coded into src/brick_slr.f90 (currently
## input the simple temperature offset in the simple_step_forward call;
## input the ais slr offset when setting sea_level_noAIS_previous and
## sea_level_noAIS).
##==============================================================================

t0.hector = 1880. # Normalize obs to this year to ~match Hector offset.
		  # Even if Hector goes back further, the change should be fairly small
dat = read.csv("data/GISTEMP_2018.csv",header=FALSE,na.strings="***",skip=4)
obs.temp = dat[,14]
obs.temp.time = dat[,1]
obs.temp = obs.temp[!is.nan(obs.temp)]
obs.temp.time = obs.temp.time[!is.nan(obs.temp)]
oitmp = which(obs.temp.time == t0.hector) : which(obs.temp.time == t0.hector + 20)
obs.temp      = obs.temp - mean(obs.temp[oitmp])

dat = read.table("data/GMSL_ChurchWhite2011_yr_2015.txt")
obs.sl        = dat[,2]/1000. # data are in mm
obs.sl.time   = dat[,1]-0.5   # they give half-year so the year is unambiguous
t0.hector = max(obs.sl.time[1],t0.hector) # SL obs only go back to 1880
oitmp = which(obs.sl.time == t0.hector) : which(obs.sl.time == t0.hector + 20)
obs.sl        = obs.sl - mean(obs.sl[oitmp])

#Normalization years for each component
t0.gsic    = 1880; tf.gsic    = 1900
t0.te      = 1880; tf.te      = 1900
t0.simple  = 1960; tf.simple  = 1990
t0.dais    = 1880; tf.dais    = 1900
t0.sl.dais = 1961; tf.sl.dais = 1990

#Mean.obs contains the normalizations for each component
mean.obs = vector("list", 5)
names(mean.obs)=as.character(c("temp.gsic","temp.te","temp.simple","temp.dais","sl.dais"))
oitmp = which(obs.temp.time == t0.gsic)    : which(obs.temp.time == tf.gsic)
mean.obs$temp.gsic[1]   = t0.gsic   ; mean.obs$temp.gsic[2]   = tf.gsic
mean.obs$temp.gsic[3]   = mean(obs.temp[oitmp])
oitmp = which(obs.temp.time == t0.te)      : which(obs.temp.time == tf.te)
mean.obs$temp.te[1]     = t0.te     ; mean.obs$temp.te[2]     = tf.te
mean.obs$temp.te[3]     = mean(obs.temp[oitmp])
oitmp = which(obs.temp.time == t0.simple)  : which(obs.temp.time == tf.simple)
mean.obs$temp.simple[1] = t0.simple ; mean.obs$temp.simple[2] = tf.simple
mean.obs$temp.simple[3] = mean(obs.temp[oitmp])
oitmp = which(obs.temp.time == t0.dais)    : which(obs.temp.time == tf.dais)
mean.obs$temp.dais[1]   = t0.dais   ; mean.obs$temp.dais[2]   = tf.dais
mean.obs$temp.dais[3]   = mean(obs.temp[oitmp])
oitmp = which(obs.sl.time   == t0.sl.dais) : which(obs.sl.time   == tf.sl.dais)
mean.obs$sl.dais[1]     = t0.sl.dais; mean.obs$sl.dais[2]     = tf.sl.dais
mean.obs$sl.dais[3]     = mean(obs.sl[oitmp]  )

print(mean.obs)
write.csv(mean.obs,"BRICK_norm_gistemp.csv",row.names=c("t0","tf","avg"))
