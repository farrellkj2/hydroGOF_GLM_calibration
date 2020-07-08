#### Example script for GLM and hydroGOF ####

# Install/load libraries, set directory and date ranges ####
#install.packages('pacman')
# devtools::install_github("CareyLabVT/GLMr")
# devtools::install_github("CareyLabVT/glmtools") 
pacman::p_load(broom, GLMr,  glmtools,  hydroGOF,  lubridate,  tidyverse)
options(scipen=999)

sim_folder  <- ('./Sunapee_GLM')

# Set model start and end dates
startDate = ymd('2003-11-08') # start date for model runs
endDate = ymd('2008-12-31') # end date for model runs

# Set calibration period start end dates (within model run dates)
calStart = ymd('2004-04-01') # calibration start
calEnd = ymd('2005-12-31') # calibration end

# Define goodness of fit function
calculate_GOF <- function(x, y, z){
  RMSE = round(rmse(x, y), 2)
  R2 = round((cor(x, y, method="pearson", use= 'na.or.complete'))^2, 2)
  NMAE = round((mae(x, y) * count(z)) / (count(z)* mean(y)),2)
  rho = round(cor(x, y, method="spearman", use= 'na.or.complete'),2)
  N = count(z)
  
  return(data.frame(N, R2, RMSE, rho, NMAE))
  
}

#### Run GLM model (note: WQ modules off in this example) ####
nml_file <- (paste0(sim_folder, "/glm2.nml"))
nml <- read_nml(nml_file)  # Read nml file
run_glm(sim_folder,  verbose=TRUE)
SimFile <- file.path(sim_folder, 'output.nc')

#### GOF metrics: Manual temperatures ####
# Load observational data
manTemp <- read_csv(paste0(sim_folder, '/observed_data/manual_temps.csv'))

# Match GLM water temperatures with field observations (resample_to_field in glmtools)
manTempResampled <- resample_to_field(SimFile, './Sunapee_GLM/observed_data/manual_temps.csv', method = 'interp') %>% 
  mutate(period = if_else((DateTime >= calStart & DateTime <= calEnd), "cal", 
                          if_else(DateTime > calEnd, "val", "NA"))) %>% filter(period != "NA")

# Calculate GOF between observed, modeled
GOFmetrics <- bind_rows(
  # Select surface depths, calculate GOF based on daily max. temperature in that depth range
  (manTempResampled %>% filter(Depth <= 4) %>% 
                           mutate(Metric = "Manual Temp Daily Max 0-4m") %>% 
                           group_by(Metric, period, DateTime) %>% 
                           summarise(MaxObs = max(Observed_temp), MaxMod = max(Modeled_temp)) %>% 
                           do(calculate_GOF(x=.$MaxMod, y= .$MaxObs, z = .))),
  # Select bottom depths, calculate GOF based on daily max. temperature in that depth range                   
  (manTempResampled %>% filter(Depth >= 16 & Depth <= 20) %>% 
                           mutate(Metric = "Manual Temp Daily Max 16-20m") %>% 
                           group_by(Metric, period, DateTime) %>% 
                           summarise(MaxObs = max(Observed_temp), MaxMod = max(Modeled_temp)) %>% 
                           do(calculate_GOF(x=.$MaxMod, y= .$MaxObs, z = .))))

# Plot sequence of graphs: observed vs. modeled along a sequence of depths
depths<-seq(1,4,1)

for(i in 1:length(depths)){
  tempdf <- subset(manTempResampled, manTempResampled$Depth==depths[i])
  plot(tempdf$DateTime, tempdf$Observed_temp, type='b', col='red',
       ylab="temperature", xlab="time", main=paste0("Obs=red,Mod=black,Depth=",depths[i]),
       ylim=c(0,30))
  points(tempdf$DateTime,tempdf$Modeled_temp, type="b", col='black')
  
}

#### GOF metrics: Buoy temperatures ####
buoyTemp <- read_csv(paste0(sim_folder, '/observed_data/buoy_temps.csv'))

buoyTempResampled <- resample_to_field(SimFile, './Sunapee_GLM/observed_data/buoy_temps.csv',  method = 'interp') %>% 
  mutate(period = if_else((DateTime >= calStart & DateTime <= calEnd), "cal", 
                          if_else(DateTime > calEnd, "val", "NA"))) %>% filter(period != "NA")
# Note: in this subset example, no buoy data during "cal" period

GOFmetrics <- bind_rows(GOFmetrics, 
                        (buoyTempResampled %>% filter(Depth <= 4) %>% 
                           mutate(Metric = "Buoy Daily Max. Temp 0-4m") %>%
                           group_by(Metric, period, DateTime) %>% 
                           summarise(MaxObs = max(Observed_temp), MaxMod = max(Modeled_temp)) %>% 
                           do(calculate_GOF(x=.$MaxMod, y= .$MaxObs, z = .))),
                        (buoyTempResampled %>% filter(Depth >= 12) %>% 
                           mutate(Metric = "Buoy Temp Daily Max 12-14m") %>% 
                           group_by(Metric, period, DateTime) %>% 
                           summarise(MaxObs = max(Observed_temp), MaxMod = max(Modeled_temp)) %>% 
                           do(calculate_GOF(x=.$MaxMod, y= .$MaxObs, z = .))))

# Plot sequence of graphs: observed vs. modeled along a sequence of depths
depths<-c(1,2,4)

for(i in 1:length(depths)){
  tempdf <- subset(buoyTempResampled, buoyTempResampled$Depth==depths[i])
  plot(tempdf$DateTime, tempdf$Observed_temp, type='b', col='red',
       ylab="temperature", xlab="time", main=paste0("Obs=red,Mod=black,Depth=",depths[i]),
       ylim=c(0,30))
  points(tempdf$DateTime,tempdf$Modeled_temp, type="b", col='black')
  
}

#### Write GOF table ####
GOFmetrics <- GOFmetrics %>% rename(NMAE = n.1) 
write_csv(GOFmetrics, paste('./GOF_metrics_Sunapee_', format(Sys.Date(), "%Y%m%d"),'.csv', sep=""), append=F)

