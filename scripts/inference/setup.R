#library(readr)
library(tidyverse)
library(data.table)
library(Matrix)
library(purrr)
library(parallel)
library(kableExtra)
library(splines)

#### 1. Setting up data ####

divis_per = 6

#### 1.1 Import and repurpose ####

# table 1 : population (correspondance between keys and names)
all_names = data.table(nom_prenom_std = unique(c(fread("data/liste_corps/liste_corps_complet.csv")$nom_prenom_std, 
                     fread("data/intermediate/correspondance_names.csv")$nom_prenom_alt_std)))
all_names$ident = 1:nrow(all_names)

population_names1 = unique(fread("data/liste_corps/liste_corps_complet.csv")[,c("nom_prenom_std")])
setnames(population_names1,"nom_prenom_std","name1")
population_names1$name2 = population_names1$name1

population_names2a = fread("data/intermediate/correspondance_names2.csv")[,c("nom_prenom_std","nom_prenom_alt_std")]
setnames(population_names2a, "nom_prenom_std","name1");setnames(population_names2a, "nom_prenom_alt_std","name2")
population_names2b = fread("data/intermediate/correspondance_names2.csv")[,c("nom_prenom_std","nom_prenom_alt_std")]
setnames(population_names2b, "nom_prenom_std","name2");setnames(population_names2b, "nom_prenom_alt_std","name1")
population_names2 = rbind(population_names2a,population_names2b)

all_connectors = rbind(population_names1, population_names2)
all_connectors = merge(all_connectors, all_names, by.x = "name2", by.y = "nom_prenom_std")
population_names = all_connectors %>% group_by(name1) %>% filter(ident == min(ident)) %>% ungroup() %>% distinct(name1, ident) %>% as.data.table()
setnames(population_names, "name1","nom_prenom_std")

# table 2 : population by public organization
population_corps = fread("data/liste_corps/liste_corps_complet_date.csv") %>%
  
  left_join(population_names, by = "nom_prenom_std") %>%
  select(-nom_prenom_std) %>%
  filter(poste%in%c("ena","dgfip","insee","cdiplo","cprefet","ce","ccomptes","dgtresor","igf","igas","iga")) %>%
  mutate(source_date_y = as.integer(substr(source_date, 1, 4))-1990,
         source_date_m = as.integer(substr(source_date, 6, 7))-1) %>%
  mutate(source_date = (12*source_date_y + source_date_m)%/%divis_per+1) %>%
  
  # Exception student ena 
  mutate(source_date = ifelse(poste == "ena" & source_date >= 56, source_date+1,source_date)) %>%

  arrange(ident, source_date) %>% distinct() %>% drop_na() %>% 
  filter(source_date <= 65) %>%
  
  as.data.table() 

# table intermédiaire: noms alternatifs (to do once)

# alternative_names <- fread("data/intermediate/jorf_orders3.csv") %>%
#   merge(population_names, by = "nom_prenom_std") %>%
#   filter(ident %in% population_corps$ident) %>%
#   select(nom_prenom_std, nom, prenom, nom_alternatif) %>%
#   filter(nom_alternatif != "") %>% unique() %>%
#   mutate(nom_prenom_alt_std = paste0(nettoyage_nom(nom), " ", nettoyage_nom(prenom)),
#          full_name_alt = paste0(str_to_lower(prenom)," ", str_to_lower(nom))) %>%
#   select(nom_prenom_std, nom_prenom_alt_std)
# 
# write_csv(alternative_names, "temp.csv")
# 
# alternative_names2 = fread("data/intermediate/correspondance_names.csv")[,c("nom_prenom_std","nom_prenom_alt_std")]
# write_csv(rbind(alternative_names, alternative_names2),"data/intermediate/correspondance_names2.csv")


# table 3 : jorf orders 
jorf_orders <- fread("data/intermediate/jorf_orders3.csv") %>%
  merge(population_names, by = "nom_prenom_std") %>%
  select(-"nom_prenom_std") %>%
  filter(ident %in% population_corps$ident) %>%

  # Definition date 
  mutate(trace_date_y = as.integer(substr(source_date, 1, 4))-1990,
         trace_date_m = as.integer(substr(source_date, 6, 7))-1) %>%
  mutate(trace_date = (12*trace_date_y + trace_date_m)%/%divis_per+1) %>%
  
  # Exception student ENA #
  mutate(trace_date = ifelse(substr(eleve_ena,1,1)=="1" & trace_date >= 56, trace_date+1, trace_date)) %>%
  
  # Filtering
  filter(trace_date <= 65) %>%

  arrange(trace_date) %>%
  drop_na(trace_date) %>%
  
  mutate(trace_type = case_when(
    depart_retraite != "" ~ 3,
    type_ordre %in% c("nomination", "admission", "affectation", "réintégration",
                      "renouvellement", "intégration", "inscription", "désignation") ~ 1,
    type_ordre %in% c("délégation de signature", "promotion", "titularisation", "habilitation",
                      "composition", "prime") ~ 2,
    TRUE ~ 0)) %>%
  
  select(ident, trace_type, trace_date) %>%
  unique()

NPER = max(jorf_orders$trace_date)

# table 4 : linkedin profiles

linkedin_profiles = fread(file = "data/intermediate/linkedin_profiles.csv") %>%
  left_join(population_names, by = "nom_prenom_std") %>%
  select(-nom_prenom_std) %>%
  
  arrange(ident, date_start_year, date_start_month) %>%
  distinct() %>%
  
  # filtering crucial missing data
  filter(!is.na(date_start_year)) %>%
  filter(str_length(company_name)>1) %>%
  
  # missing data: last date not described (so we take the last one available)
  mutate(
    date_end_month = ifelse((is.na(date_end_year) & date_start_year >= 2012), 5, date_end_month),
    date_end_year = ifelse((is.na(date_end_year) & date_start_year >= 2012), 2022, date_end_year)) %>%

  # missing data: no month described 
  mutate(date_start_month = ifelse(is.na(date_start_month), 1, date_start_month),
         date_end_month = case_when(is.na(date_end_month) & date_end_year < 2022 ~ 12, 
                                    is.na(date_end_month) & date_end_year == 2022 ~ 5,
                                    T ~ date_end_month)) %>%
  
  # Conversion to quantitative format 
  mutate(date_start_y = date_start_year - 1990, date_end_y = date_end_year - 1990,
         date_start_m = date_start_month - 1, date_end_m = date_end_month - 1) %>%
  mutate(date_start = (date_start_y*12 + date_start_m)%/%divis_per+1,
         date_end = (date_end_y*12 + date_end_m)%/%divis_per+1) %>%
  
  # Correction de date 
  filter(date_start <= NPER) %>%
  mutate(date_end = ifelse(date_end >= NPER, NPER, date_end)) %>%
  
  # Simplifying # 
  mutate(date_end = ifelse(is.na(date_end), NPER, date_end)) %>% # handling the censoring 
  mutate(check = T) %>%
  distinct(ident, nature, company_name, date_start, date_end, check, profile_id) %>%
  
  # Filter on org type & date
  filter(nature %in% c("public","private")) %>% 
  filter(date_end >= 1) %>% 
  mutate(date_start = ifelse(date_start < 1, 1, date_start)) %>%
  
  mutate(nature = case_when(
    nature == "public" ~ 0,
    nature == "private" ~ 1))

# Profile information
population_corps = left_join(population_corps, as.data.table(distinct(linkedin_profiles, ident))[,profile:=1], by = "ident")
population_corps$profile[is.na(population_corps$profile)] = 0

#### 1.2 Functions ####

aggregating_intervals = function(x, edges = c(1,NPER)){
  
    rp = rep(98, NPER)
    
    for (nature_ in c(1,0,98)){ # priority to public sector job
      x2 = x[nature == nature_,]
      if (nrow(x2) == 0){next}
      for (i in 1:nrow(x2)){
        sel = x2[[i,"date_start"]]:x2[[i,"date_end"]]
        rp[sel] = nature_
      }
    }
    
    rp[setdiff((1:NPER), (edges[1]:edges[2]))] = 99
    rp2 = rle(rp)
    
    return(list("full"=rp, "agg"=rp2))
  }
  
  # pre-computing
  spans_max = jorf_orders[,.(begin = min(trace_date), end = min(ifelse(trace_type==3, trace_date, NPER))),.(ident)]
  traces = jorf_orders[trace_type %in% c(1,2),.(N=ifelse(.N>0,1,0)), by = .(ident, trace_date)]
  
  generate_data = function(id){
  
  # Determining main organization
  x1 = population_corps[ident==id,]
  x2 = rep(0, NPER)
  for (i in 1:nrow(x1)){
    x2[x1[[i,"source_date"]]] = x1[[i,"poste"]]
  }
  
  k = "NA"
  x3 = rle(x2)
  for (i in 1:length(x3$values)){
    if (x3$values[i]=="0"){
      x3$values[i]=k
    } else {
      k = x3$values[i]
    }
  }
  main_org = inverse.rle(x3)
  
  # Spans
  span_max = spans_max[ident==id,]$begin:spans_max[ident==i,]$end
  span_min = max(which.min(main_org=="NA"), min(span_max)):max(span_max)
  span_ena = max(which.max(main_org=="ena"), min(span_max)):max(span_max)
  
  # Main org def
  main_org = main_org[span_min]
  
  # transition between main org
  transition_group = as.double(lag(main_org)!=main_org)
  transition_group[1]=1
  
  # JO orders
  trace = rep(0, NPER)
  x = traces[ident == id,]
  trace[x$trace_date] = as.character(x$N)
  trace_obs = as.character(trace[span_min])
  
  # Trajectories 
  traj_obs = aggregating_intervals(linkedin_profiles[ident==id,])$full[span_min]
  
  # Period covariate
  period = (span_min-1)/(NPER-1)
  
  # Age covariates : maximum, minimum, ENA
  age_max = (span_min-min(span_max))/(NPER-1)
  age_ena = (span_min-min(span_ena))/(NPER-1)
  age_ena[age_ena<0]=NA
  if (!('ena'%in% main_org)) age_ena = rep(NA, length(age_ena))
  
  # profile info 
  profile_present = unique(population_corps[ident==id,profile])
  
  # correct for incoherence
  trace_obs[transition_group == 1]=1
  trace_obs[traj_obs==1 & trace_obs==1]=0
  
  # sex (fill_na = M)
  sex = x1$sexe[1]
  
  w = (data.table(
    ident = id,
    group = main_org,
    traj_obs = traj_obs,
    trace_obs = as.double(trace_obs),
    period = period,
    sex = sex,
    age_max = age_max,
    age_ena = age_ena,
    profile_present = profile_present))
  
  return(w[period<=1,])
}

#### 1.3 Applying functions ####

# V2
ncl = detectCores()-1
cl = makeCluster(ncl, type = "FORK")
clusterSetRNGStream(cl)

names_ = as.character(spans_max$ident[order(spans_max$ident)])
names_ = names_[!is.na(names_)]

#for (id in names_) {print(id);generate_data(id)}

res = parLapply(cl,names_,fun=generate_data)
res2 = list()
for (id in 1:length(res)){res2[[names_[id]]]=res[[id]]} 

saveRDS(res2, file = "data/clean/dataTot.rds")

