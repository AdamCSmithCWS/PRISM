
# PRISM population size estimates -----------------------------------------


# Packages ----------------------------------------------------------------

library(tidyverse)
library(jagsUI)
library(lubridate)
library(dclone)
library(ggforce)
library(ggmcmc)

# Data --------------------------------------------------------------------

data.dir <- "data/"

int.i = read.csv(paste0(data.dir,"IntensivePlots_Int.csv"),stringsAsFactors = F)
int.r = read.csv(paste0(data.dir,"IntensivePlots_Rapid.csv"),stringsAsFactors = F)
rap = read.csv(paste0(data.dir,"RapidPlots.csv"),stringsAsFactors = F)
rap$month <- month(mdy(rap$survey_date))
rap$day <- day(mdy(rap$survey_date))
earliest_survey <- min(rap[which(rap$month == 6),"day"])
rap$days_s_June4 <- as.integer(mdy(rap$survey_date)-mdy(paste("June",earliest_survey-1,year(mdy(rap$survey_date)))))


areas = read.csv(paste0(data.dir,"StratumAreas.csv"),stringsAsFactors = F)


# 
 table(aggregate(plot ~ X11.UniqueZone, data = rap, FUN = function(x) length(unique(x)))$plot)


# Data Prep ---------------------------------------------------------------


# complete list of Rapid surveys conducted --------------------------------

survs = unique(rap[,c("X1.Region","plot","survey_date",
                      "X7.Hab","X11.UniqueZone","year",
                      "type","survey_method","start_time1",
                      "area","surveyID")])
nrapid = nrow(survs)

# species list ------------------------------------------------------------

shorebirds = unique(rap[which(rap$group == "Shorebirds"),"species"])

n_non_zeros = sort(table(rap$species)[shorebirds],decreasing = T)

p_non_zeros = n_non_zeros/nrapid

minimum_proportion = 0.03 #change this to add less prevalent species
species_to_inc = names(p_non_zeros[p_non_zeros > minimum_proportion]) 
# including only species present in > 10% of the surveys



# all surveys and species -------------------------------------------------

survs_all = expand_grid(survs,species = species_to_inc)
rap_counts_inc = rap[which(rap$species %in% species_to_inc),c("X1.Region","plot","survey_date",
                                                              "X7.Hab","X11.UniqueZone","year",
                                                              "type","survey_method","start_time1",
                                                              "area","surveyID","species","number")]
all_data = full_join(survs_all,rap_counts_inc)
all_data[which(is.na(all_data$number)),"number"] <- 0

all_data$j_count <- as.integer(all_data$number)
all_data$j_species <- as.integer(factor(all_data$species))
all_data$j_habitat <- all_data$X7.Hab
all_data$j_sub_region <- as.integer(factor(all_data$X1.Region))
all_data$j_cluster <- as.integer(factor(all_data$X11.UniqueZone))
nplots = NA
for(k in 1:max(all_data$j_cluster)){
  wl <- which(all_data$j_cluster == k)
  all_data[wl,"j_plot"] <- as.integer(factor(all_data$plot[wl]))
  nplots[k] <- max(all_data[wl,"j_plot"])
}






# Data objects required for the rapid survey portion of the model ---------


#nrapid = number of counts
nrapid = nrow(all_data)
### following vectors all nrapid in length
#rapid_count = counts at rapid plots by species
rapid_count = all_data$j_count
#species = species integer ID
species = all_data$j_species
#habitat = habitat integer ID
habitat = all_data$j_habitat
#sub_region = subregion integer ID
sub_region = all_data$j_sub_region
#cluster = cluster of plots integer ID
cluster = all_data$j_cluster
#plot = plot integer ID
plot = all_data$j_plot
#log_area = log(surveyed area of the plot)
log_area = log(all_data$area)

#nhabitats = number of habitat types (3)
nhabitats = 3
#nsub_regions = number of sub_regions
nsub_regions = max(sub_region)
#nplots[cluster] = number of rapidly surveyed plots with each cluster
nplots = nplots
#nclusters = total number of rapidly surveyed plot-clsuters
nclusters = max(all_data$j_cluster)
#nspecies = number of species to analyse
nspecies = max(all_data$j_species)







# Intenstive survey prep --------------------------------------------------

int.i = int.i[,c("plotID","number","species","subregion_code","year","plot_area")] 


int.r = int.r[,c("plotID","number","species","subregion_code","year","plot_area","surveyID")] 

intsurvs = unique(int.r[,c("plotID","subregion_code","year","plot_area","surveyID")]) 

intsurvs = expand_grid(intsurvs,species = species_to_inc)

int.i2 = int.i[which(int.i$species %in% species_to_inc),]
int.r2 = int.r[which(int.r$species %in% species_to_inc),]


intsurvs = full_join(intsurvs,int.r2)
intsurvs[which(is.na(intsurvs$number)),"number"] <- 0
intsurvs = rename(intsurvs,count_r = number)  
 

intsurvs = full_join(intsurvs,int.i2)
intsurvs[which(is.na(intsurvs$number)),"number"] <- 0
intsurvs = rename(intsurvs,count_i = number)  

intsurvs[which(is.na(intsurvs$count_r)),"count_r"] <- 0

sp_names_indices = unique(all_data[,c("species","j_species")])

intsurvs <- full_join(intsurvs,sp_names_indices)


  
# data objects required for hte visibility portion of the model -----------


# nintensive = number of intensive and rapid survey pairs, each survey on a given plot treated as independent - i.e., simple random survey of plots following Bart and Johnston
nintensive = nrow(intsurvs)
# count_i
count_i = ceiling(as.integer(intsurvs$count_i))
# count_r
count_r = intsurvs$count_r
# species_r
species_r = intsurvs$j_species





# Area scaling factors ----------------------------------------------------

all_sub_regions <- unique(all_data[,c("X1.Region","j_sub_region")])
names(all_sub_regions)[1] <- "Region"
areas <- full_join(all_sub_regions,areas)

areas_final <- areas[which(!is.na(areas$j_sub_region)),]
areas_final <- areas_final[order(areas_final$j_sub_region),]

areas_mat = matrix(0,nrow = nsub_regions,ncol = nhabitats)
for(i in 1:nsub_regions){
  for(h in 1:nhabitats){
    areas_mat[i,h] <- as.numeric(areas_final[which(areas_final$j_sub_region == i & areas_final$Habitat == h),"Area"])
  }
}
areas_mat[is.na(areas_mat)] <- 0
# compile jags data object ------------------------------------------------


jags_data = list(
  nrapid = nrapid,
  rapid_count = rapid_count,
  species = species,
  habitat = habitat,
  sub_region = sub_region,
  #cluster = cluster,
  #plot = plot,
  log_area = log_area,
  nhabitats = 3,
  nsub_regions = nsub_regions,
  nplots = nplots,
  nclusters = nclusters,
  nspecies = nspecies,
  nintensive = nintensive,
  count_i = count_i,
  count_r = count_r,
  species_r = species_r,
  areas_mat = areas_mat
)



# Model Description -------------------------------------------------------

### investigate:
####  the necessity of the over-dispersion parameter which is correlated with site-effects
####  the convergence of the habitat parameters (particularly among regions)
####  

sink("PRISM_model.txt")
cat("
model
{

#### likelihood for rapid surveys
for(i in 1:nrapid){

rapid_count[i] ~ dpois(lambda[i])

log(lambda[i]) <- elambda[i] 

#elambda[i] <- log_area[i] + alpha[species[i]] + strat[habitat[i],sub_region[i],species[i]] + site[plot[i],cluster[i],species[i]] + noise[i]
elambda[i] <- log_area[i] + alpha[species[i]] + strat[habitat[i],sub_region[i],species[i]] + noise[i]

#noise[i] ~ dnorm(0,taunoise)
noise[i] ~ dt(0,taunoise,nu) #alternate heavy-tailed

}#i



### precision priors
nu ~ dgamma(2, 0.1) #alternative degrees of freedom (i.e., the heavy-tail component of the t-distribution), if nu is large (infinite) this is equivalent to the normal distribution noise

#taunoise ~ dgamma(0.001,0.001) #extra Poisson variance on counts
taunoise ~ dscaled.gamma(1,100) #implicit prior on sigma of a half-t dist: sigma = 1*t(df = 100) , i.e., 95% prob sd < 2.3
sd_noise <- taunoise^-0.5
#varnoise <- 1/taunoise
adj <- (1.422*pow(nu,0.906))/(1+(1.422*pow(nu,0.906)))
varnoise <- (sd_noise/adj)^2
hab[1] <- 0  
for(h in 2:nhabitats){
  hab[h] ~ dnorm(0,0.1) #mean habitat effects for a species
}

for(s in 1:nspecies){

alpha[s] ~ dnorm(0,0.1)#### species means - fixed effects


#### main effects of habitat by stratum and species

tau_hab_region[s] ~ dscaled.gamma(1,100) #implicit prior on sigma of a half-t dist: sigma = 1*t(df = 100) , i.e., 95% prob sd < 2.3
sd_hab_region[s] <- tau_hab_region[s]^-0.5

#tau_hab_region[s] ~ dgamma(0.001,0.001) #prior on species specific variance

for(h in 1:nhabitats){
  #hab[h,s] ~ dnorm(0,0.1) #mean habitat effects for a species
  for(r in 1:nsub_regions){
    #strat[h,r,s] ~ dnorm(hab[h,s],tau_hab_region) ## alternative without species specific variance
    strat[h,r,s] ~ dnorm(hab[h],tau_hab_region[s]) #stratum, habitat, and species specific intercepts, centered on the species intercepts

  }#r
}#h




#### random effects of site within clusters (zones)
#tau_cluster[s] ~ dgamma(0.001,0.001) #
#tau_site_cluster[s] ~ dgamma(0.001,0.001)
#tau_cluster[s] ~ dscaled.gamma(1,100) #implicit prior on sigma of a half-t dist: sigma = 1*t(df = 100) , i.e., 95% prob sd < 2.3
#sd_cluster[s] <- tau_cluster[s]^-0.5
#tau_site_cluster[s] ~ dscaled.gamma(1,100) #implicit prior on sigma of a half-t dist: sigma = 1*t(df = 100) , i.e., 95% prob sd < 2.3
#sd_site_cluster[s] <- tau_site_cluster[s]^-0.5




#v_cluster[s] <- 1/tau_cluster[s]
#v_site_cluster[s] <- 1/tau_site_cluster[s]

# for(k in 1:nclusters){
#   #clust[k,s] ~ dnorm(0,tau_cluster[s])
#   for(p in 1:nplots[k]){
#     site[p,k,s] ~ dnorm(0,tau_site_cluster[s]) #stratum, habitat, and species specific intercepts, centered on the species intercepts
#       #site[p,k,s] ~ dnorm(clust[k,s],tau_site_cluster[s]) #stratum, habitat, and species specific intercepts, centered on the species intercepts
#   }#p
# }#h



#### Derived parameters


  for(r in 1:nsub_regions){
for(h in 1:nhabitats){

    #n[s,h,r] <- areas_mat[r,h] * exp(alpha[s] + strat[h,r,s] + sr[s] + 0.5*v_cluster[s] + 0.5*v_site_cluster[s] + 0.5*varnoise)
    #n_uncor[s,h,r] <- exp(alpha[s] + strat[h,r,s] + sr[s] + 0.5*v_cluster[s] + 0.5*v_site_cluster[s] + 0.5*varnoise)
     n[s,h,r] <- areas_mat[r,h] * exp(alpha[s] + strat[h,r,s] + sr[s] + 0.5*varnoise)
    n_uncor[s,h,r] <- exp(alpha[s] + strat[h,r,s] + sr[s] + 0.5*varnoise)
    
}#h
Nr[s,r] <- sum(n[s,1:nhabitats,r])
}#r
N[s] <- sum(Nr[s,1:nsub_regions])

}#s











#### intensive survey parameters for visibility


for(i in 1:nintensive){
 ## raw counts
 count_i[i] ~ dpois(lambda_i[i])
 count_r[i] ~ dpois(lambda_r[i])
 ## poisson counts at sites
 
 ## species specific counts at each intensive survey
  mu_i[i] <- si[species_r[i]] #species mean expected count
	loglambda_i[i] ~ dnorm(mu_i[i],taunoise_i)
	log(lambda_i[i]) <- loglambda_i[i]
	
 ## species specific counts at each rapid survey of intensive plots
   mu_r[i] <- loglambda_i[i] - sr[species_r[i]] #expected count for that species in the intensive survey minus a correction value
	loglambda_r[i] ~ dnorm(mu_r[i],taunoise_r)
	log(lambda_r[i]) <- loglambda_r[i]


}#i



## priors
taunoise_i ~ dgamma(0.001,0.001) # count variance on intensive surveys
taunoise_r ~ dgamma(0.001,0.001) # count variance on rapid surveys
#taustratspecies ~ dgamma(0.001,0.001) #variance among strata on species true abundance in intensive surveys 

tau_sr ~ dgamma(0.001,0.001) # variance among species in the detection ratio
#taustrat_ratio ~ dgamma(0.001,0.001) #variance among strata on species level detection ratio
SR ~ dnorm(0,0.1) #mean detection ratio across species 

## visibility corrections
	for(s in 1:nspecies){
		si[s] ~ dnorm(0.0,0.1)
		sr[s] ~ dnorm(SR,tau_sr)

		  # si[s,sr] ~ dnorm(mu_si[s],taustratspecies)
		  # sr[s,sr] ~ dnorm(mu_sr[s],taustrat_ratio)
		  # 
		  vis[s] <- exp(sr[s]) #species by region visibility adjustments on the multiplicative scale
	}#s
	all_muvis <- exp(SR)


}",fill = TRUE)
sink()



# Run the model -----------------------------------------------------------


nburn = 2000
nsave = 5000
nthin = 10

mod = jags(data = jags_data,
           model.file = "PRISM_model.txt",
           parameters.to.save = c("vis","all_muvis","sd_sr",
                                  "n","n_uncor",
                                  "N","Nr",
                                  "nu",
                                  #"sd_cluster",
                                  "hab",
                                  "strat",
                                  #"sd_site_cluster",
                                  "sd_noise",
                                  "sd_hab_region"),
           n.chains = 3,
           n.iter = (nsave*nthin)+nburn,
           n.burnin = nburn,
           n.thin = nthin,
           parallel = T)




out = data.frame(mod$summary)  
out$node = row.names(out)
names(out)[3:7] <- c("lci","lqrt","med","uqrt","uci")



# visibility plot ---------------------------------------------------------


out_vis = filter(out,node %in% paste0("vis[",1:nspecies,"]"))
out_vis$j_species = 1:nrow(out_vis)
out_vis = left_join(out_vis,sp_names_indices)


visib = ggplot(data = out_vis,aes(x = species,y = med))+
  geom_linerange(aes(x = species,ymin = lci,ymax = uci),alpha = 0.5)+
  geom_point()+
  ylab("Multiplicative visibility correction factor")+
  geom_hline(yintercept = 1)+
  coord_cartesian(ylim = c(0.25,1.5))+
  labs(title = "Detectability ratios all indicate more birds observed during rapid surveys")

pdf("visibility correction factors.pdf")
print(visib)  
  dev.off()
  
  

# regional abundance plot by habitat -------------------------------------------------

  
  
  out_n = filter(out,grepl(node,pattern = "n[",fixed = T))
  out_n$j_species = as.integer(str_match(out_n$node,"(?<=n[:punct:])[:digit:]+"))
  out_n$habitat = as.integer(str_match(out_n$node,"(?<=n[:punct:][:digit:]{1,2}[:punct:])[:digit:]+"))
  out_n$j_sub_region = as.integer(str_match(out_n$node,"(?<=n[:punct:][:digit:]{1,2}[:punct:][:digit:]{1,2}[:punct:])[:digit:]+"))
  
  out_n = left_join(out_n,sp_names_indices)
  out_n = left_join(out_n,all_sub_regions)
  out_n$Region = factor(out_n$Region)
  out_n$habitat = factor(out_n$habitat)
  
  pdf("sub_regional population size estimates by habitat.pdf",
      width = 11,
      height = 8.5)
for(pg in 1:ceiling(nspecies/6)){
  nbyreg_s = ggplot(data = out_n,aes(x = Region,y = med,colour = habitat))+
    geom_linerange(aes(x = Region,ymin = lci,ymax = uci),alpha = 0.5,position = position_dodge(width = 0.6))+
    geom_point(position = position_dodge(width = 0.6))+
    facet_wrap_paginate(facets = ~species,scales = "free",ncol = 2,nrow = 3,page = pg)
  print(nbyreg_s)
  }


  dev.off()
  
  
  

# regional density plot uncorrected for area ----------------------
# add mean observed
  
  

  out_n = filter(out,grepl(node,pattern = "n_uncor[",fixed = T))
  out_n$j_species = as.integer(str_match(out_n$node,"(?<=n_uncor[:punct:])[:digit:]+"))
  out_n$habitat = as.integer(str_match(out_n$node,"(?<=n_uncor[:punct:][:digit:]{1,2}[:punct:])[:digit:]+"))
  out_n$j_sub_region = as.integer(str_match(out_n$node,"(?<=n_uncor[:punct:][:digit:]{1,2}[:punct:][:digit:]{1,2}[:punct:])[:digit:]+"))
  
  out_n = left_join(out_n,sp_names_indices)
  out_n = left_join(out_n,all_sub_regions)
  out_n$Region = factor(out_n$Region)
  out_n$habitat = factor(out_n$habitat)
  for(j in 1:nrow(out_n)){
    ww = which(all_data$X1.Region == out_n[j,"Region"] &
                 all_data$X7.Hab == out_n[j,"habitat"] &
                 all_data$species == out_n[j,"species"])
    out_n[j,"obs_mean"] <- mean(as.numeric(unlist(all_data[ww,"number"]/all_data[ww,"area"])),na.rm = T)
      }

  pdf("sub_regional density estimates not corrected for area.pdf",
      width = 11,
      height = 8.5)
  for(pg in 1:ceiling(nspecies/6)){
    nbyreg_s = ggplot(data = out_n,aes(x = Region,y = med,colour = habitat))+
    geom_linerange(aes(x = Region,ymin = lci,ymax = uci),alpha = 0.5,position = position_dodge(width = 0.6))+
    geom_point(position = position_dodge(width = 0.6))+
      geom_point(data = out_n,aes(x = Region,y = obs_mean,colour = habitat),shape = 2,alpha = 0.5,position = position_dodge(width = 0.6))+
      facet_wrap_paginate(facets = ~species,scales = "free",ncol = 2,nrow = 3,page = pg)
    print(nbyreg_s)
  }
    dev.off()
  
  
  
  
  
  
  
  # regional abundance plot  ----------------------
  
  
  out_n = filter(out,grepl(node,pattern = "Nr[",fixed = T))
  out_n$j_species = as.integer(str_match(out_n$node,"(?<=Nr[:punct:])[:digit:]+"))
  out_n$j_sub_region = as.integer(str_match(out_n$node,"(?<=Nr[:punct:][:digit:]{1,2}[:punct:])[:digit:]+"))
  
  out_n = left_join(out_n,sp_names_indices)
  out_n = left_join(out_n,all_sub_regions)
  out_n$Region = factor(out_n$Region)
 
  nbyreg_s = ggplot(data = out_n,aes(x = Region,y = med))+
    geom_linerange(aes(x = Region,ymin = lci,ymax = uci),alpha = 0.5)+
    geom_point()+
    facet_wrap(facets = ~species,scales = "free")
  
  pdf("sub_regional population size estimates.pdf",
      width = 11,
      height = 8.5)
  for(pg in 1:ceiling(nspecies/6)){
    nbyreg_s = ggplot(data = out_n,aes(x = Region,y = med))+
    geom_linerange(aes(x = Region,ymin = lci,ymax = uci),alpha = 0.5)+
    geom_point()+
      facet_wrap_paginate(facets = ~species,scales = "free",ncol = 2,nrow = 3,page = pg)
    print(nbyreg_s)
  }
  dev.off()
  
  
  
  
  # Total estimated populations ---------------------------------------------------------
  
  
  out_N = filter(out,node %in% paste0("N[",1:nspecies,"]"))
  out_N$j_species = 1:nrow(out_N)
  out_N = left_join(out_N,sp_names_indices)
  
  
  visib = ggplot(data = out_N,aes(x = species,y = med))+
    geom_linerange(aes(x = species,ymin = lci,ymax = uci),alpha = 0.5)+
    geom_point()+
    ylab("Population size")+
    geom_hline(yintercept = 1)+
    #coord_cartesian(ylim = c(0.25,1.5))+
    labs(title = "Total estimated populations")
  
  pdf("Total estimated populations.pdf")
  print(visib)  
  dev.off()
  
  

# MCMC explore ------------------------------------------------------------

conv = ggs(mod$samples)  
  
  
  
  conv_sd = filter(conv,grepl(pattern = "sd",Parameter))
  
ggmcmc(conv_sd,file="full convergence on sd parameters.pdf", param_page=10)

  
  
  
habs = filter(conv,grepl(pattern = "hab[",Parameter,fixed = T))
ggmcmc(habs,file="full convergence on hab parameters.pdf", param_page=10)

  

strats = filter(conv,grepl(pattern = "strat[",Parameter,fixed = T))
ggmcmc(strats,file="full convergence on strats parameters.pdf", param_page=10)

# raw data distributions --------------------------------------------------

all_data$Habitat = factor(all_data$X7.Hab)
all_data_nonz = filter(all_data,number > 5)

  raw.p = ggplot()+
    geom_violin(data = all_data,aes(x = Habitat, y = number))+
    #geom_dotplot(data = all_data_nonz,aes(x = Habitat, y = number),binaxis = "y",binwidth = 1,stackdir = "center",dotsize = 1,alpha = 0.3)+
    facet_wrap(facets = ~species,nrow = 3,ncol = 5,scales = "free")
  
  
  pdf(file = "violin plots of the raw counts by species and habitat.pdf",
      width = 11,
      height = 8.5)
  print(raw.p)
  dev.off()
  
  
  
  
  
  

  # #### intensive survey parameters for visibility
  # 
  # 
  # for(i in 1:nintensive){
  #   ## raw counts
  #   count_i[i] ~ dpois(lambda_i[i])
  #   count_r[i] ~ dpois(lambda_r[i])
  #   ## poisson counts at sites
  #   
  #   ## species specific counts at each intensive survey
  #   mu_i[i] <- si[species_r[i],sub_region_r[i]] #species mean expected count in that subregion
  #   loglambda_i[i] ~ dnorm(mu_i[i],taunoise_i)
  #   log(lambda_i[i]) <- loglambda_i[i]
  #   
  #   ## species specific counts at each rapid survey of intensive plots
  #   mu_r[i] <- loglambda_i[i] - sr[species_r[i],sub_region_r[i]] #expected count for that species in the intensive survey minus a correction value
  #   loglambda_r[i] ~ dnorm(mu_r[i],taunoise_r)
  #   log(lambda_r[i]) <- loglambda_r[i]
  #   
  #   
  # }#i
  # 
  # 
  # 
  # ## priors
  # taunoise_i ~ dgamma(0.001,0.001) # count variance on intensive surveys
  # taunoise_r ~ dgamma(0.001,0.001) # count variance on rapid surveys
  # taustratspecies ~ dgamma(0.001,0.001) #variance among strata on species true abundance in intensive surveys 
  # 
  # tau_sr ~ dgamma(0.001,0.001) # variance among species in the detection ratio
  # taustrat_ratio ~ dgamma(0.001,0.001) #variance among strata on species level detection ratio
  # SR ~ dnorm(0,0.1) #mean detection ratio across species 
  # 
  # ## visibility corrections
  # for(s in 1:nspecies){
  #   mu_si[s] ~ dnorm(0.0,0.1)
  #   mu_sr[s] ~ dnorm(SR,tau_sr)
  #   for(sr in 1:nsub_region){
  #     si[s,sr] ~ dnorm(mu_si[s],taustratspecies)
  #     sr[s,sr] ~ dnorm(mu_sr[s],taustrat_ratio)
  #     
  #     vis[s,sr] <- exp(sr[s,sr]) #species by region visibility adjustments on the multiplicative scale
  #   }#sr
  #   muvis[s] <- exp(mu_sr[s])
  # }#s
  # all_muvis <- exp(SR)
  # 









