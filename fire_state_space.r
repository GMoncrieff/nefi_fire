library(ggplot2)
library(coda)
library(rjags)
library(dplyr)
library(tidyr)
library(ecoforecastR)

data_raw <- read.csv('nefi_fire.csv')

#lets only use one region for now
data_raw <- data_raw[data_raw$region=="CB",]

#reshape to a matrix with ncols = timesteps and nrow = sites
data_mat <- data_raw %>%
  select(id,year,NDVI) %>% 
  spread(year,NDVI) 

#get df to match rows to site names for later
data_id <- data_mat %>%
  mutate(site=as.numeric(row.names(.))) %>%
  select(site,id)

#convert to mat for jags
data_mat <- data_mat %>%
  select(-id) %>%
  as.matrix(.)

#n sites
NS = nrow(data_mat)
#n time
NT = ncol(data_mat)

###JAGS model

random_site_walk = "
model{

#### Priors
tau_obs ~ dgamma(0.01,0.00000001) ## Observation error precision
tau_add ~ dgamma(0.01,0.00000001) ## Process errror precision 

#### initial conditions: t = time, s = site
for(s in 1:NS){
x[s,1] ~ dnorm(0.5,1)           ## prior on Initial condition
}

#### process model: random walk

for(t in 2:NT){
for(s in 1:NS){
x[s,t]~dnorm(x[s,t-1],tau_add)
}
}

### observation model
for(t in 1:NT){
for(s in 1:NS){
y[s,t] ~ dnorm(x[s,t],tau_obs)
}
}
}
"

# Create new data object, and set initial conditions using random samples from this data:
data <- list(y=data_mat,
             NS = NS,
             NT = NT
)
# Send the model, the data, and the initial conditions to JAGS, which will return the JAGS model object:
j.model.rw   <- jags.model (file = textConnection(random_site_walk),
                            data = data,
                            n.chains = 3)

var.out   <- coda.samples (model = j.model.rw,
                           variable.names = "x",
                           n.iter = 5000)

#var.mat      <- as.matrix(var.out)
#GBR <- gelman.plot(var.out)
#pairs(var.mat)

jags.burn <- window(var.out,start=2000) #discard burnin
#calc ci
ci <- apply(as.matrix(jags.burn),2,quantile,c(0.025,0.5,0.975))

#tranformn results to long
cidf <- as.data.frame(ci) %>%
  mutate(ci = row.names(.)) %>%
  gather("lab","x",-ci) %>%
  mutate(lab = substr(lab, 1, nchar(lab)-1)) %>%
  mutate(lab = substr(lab, 3, nchar(lab))) %>%
  separate(lab,into=c("site","time"),sep=',') %>%
  mutate(site=as.numeric(site),time=as.numeric(time)) %>%
  spread(ci,x) %>%
  rename(lower=`2.5%`,middle=`50%`,upper=`97.5%`)

#transform data to long
datadf <- data_mat %>%
  as.data.frame(.) %>%
  mutate(site=rownames(.)) %>%
  gather("year","obs",-site) %>%
  mutate(time = as.numeric(as.factor(year)),
         year=as.numeric(year),
         site=as.numeric(site))

#join to get years and sites
cijoin <- cidf %>% 
  left_join(datadf,by=c('time','site')) %>%
  left_join(data_id,by='site')

#select subset
ciplot <- cijoin %>%
  filter(site<=10)

#plot!
ggplot(ciplot, aes(x = time,y = obs)) +
  geom_point() +
  geom_line(aes(x = time,y = middle),color="red",alpha=0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper),alpha=0.2,fill="red") +
  xlab("time") +
  ylab("NDVI") +
  facet_wrap(~site)+
  theme_bw()
