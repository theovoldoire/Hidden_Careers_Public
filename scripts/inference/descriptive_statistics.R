
library(purrr)
library(data.table)
library(Matrix)
library(HMM)
library(parallel)
library(mgcv)
library(kableExtra)
library(dplyr)
library(splines)

ncl = detectCores()-1
cl = makeCluster(ncl, type = "FORK")
clusterSetRNGStream(cl)

source("scripts/inference/inference_functions.R")
clusterEvalQ(cl, {source("scripts/inference/inference_functions.R")})

dataTot = readRDS("data/clean/dataTot.rds")

#### Table 1 : cardinality pop 1 ####

NPER = 65
dataTot2 = parLapply(cl, dataTot, clean1, 
                     groups = c("ena"), ena_pop = T, remove_early = F,
                     linkedin_correct = 1, data_augmentation = 0.01)
dataPop = list()
for (x in names(dataTot)) if (nrow(dataTot2[[x]])>1) dataPop[[x]] = dataTot2[[x]]

x = rbindlist(dataPop)[,c("ident","period","age_ena")]
x$cohort_ena = (x$period-(x$age_ena+0.0625))*(NPER-1)

tab1 = table(x[age_ena==0,]$cohort_ena)
names(tab1) = as.character(1990:2019)
tab1 = as.data.frame(tab1)
colnames(tab1) = c("Year","Count")
kbl(cbind(tab1[1:10,],tab1[11:20,],tab1[21:30,]), format = "latex", booktabs = T,
    linesep = "",caption="Number of students by year of admission.") |>
  save_kable(file = "documents/figures/population_ena_students_year.tex", booktabs = TRUE)

#### Table 2 : cardinality pop 2 ####

# OK, donc population finale: on filtre pour <= 54

dataTot2 = parLapply(cl, dataTot, clean1, 
                     groups = c("cprefet","ccomptes","iga","ce","cdiplo","igf","dgtresor",
                                "dgfip","igas","insee"), ena_pop = F, remove_early = F,
                     linkedin_correct = 1, data_augmentation = 0.01)
dataPop = list()
for (x in names(dataTot)) if (nrow(dataTot2[[x]])>1) dataPop[[x]] = dataTot2[[x]]

tab2 = as.data.table(table(unique(rbindlist(dataPop)[,c("group","ident")])$group))
tab2 = tab2[order(-N)]
colnames(tab2) = c("group","count")

kbl(cbind(tab2[1:5],tab2[6:10]), format = "latex",booktabs = T,
    caption = "Number of individuals affiliated by service (1990-2022).")|>
  save_kable(file = "documents/figures/population_important_groups.tex")


#### Table 3 : proportion of profiles ####

dataTot2 = parLapply(cl, dataTot, clean1, 
                     groups = c("ena"), ena_pop = T, remove_early = F,
                     linkedin_correct = 1, data_augmentation = 0.01)
dataPop = list()
for (x in names(dataTot)) if (nrow(dataTot2[[x]])>1) dataPop[[x]] = dataTot2[[x]]

rr1 =  unique(rbindlist(dataPop)[,c("ident","group","profile_present")])

dataTot2 = parLapply(cl, dataTot, clean1, 
                     groups = c("cprefet","ccomptes","iga","ce","cdiplo","igf","dgtresor",
                                "dgfip","igas","insee"), ena_pop = F, remove_early = F,
                     linkedin_correct = 1, data_augmentation = 0.01)
dataPop = list()
for (x in names(dataTot)) if (nrow(dataTot2[[x]])>1) dataPop[[x]] = dataTot2[[x]]

rr2 =  unique(rbindlist(dataPop)[,c("ident","group","profile_present")])

rr = rbind(rr1,rr2)
tab3 = rr %>% group_by(group) %>% summarise(prop = mean(profile_present)) %>% arrange(desc(prop)) %>%
  mutate(prop=as.character(round(prop*100,1)))

tab3a = tab3[1:6,]
tab3b = tab3[7:11,] %>% add_row("group"="",prop="")
kbl(cbind(tab3a,tab3b), format = "latex",booktabs = T,linesep = "",
    caption = "Proportion of individuals by group with a matched profile.")|>
  save_kable(file = "documents/figures/proportion_profiles.tex")


#### Table 6: Men/ENA for study 2 #### 


dataTot2 = parLapply(cl, dataTot, clean1, 
                     groups = c("cprefet","ccomptes","iga","ce","cdiplo","igf","dgtresor",
                                "dgfip","igas","insee"), ena_pop = F, remove_early = F,
                     linkedin_correct = 1, data_augmentation = 0.01)
dataPop = list()
for (x in names(dataTot)) if (nrow(dataTot2[[x]])>1) dataPop[[x]] = dataTot2[[x]]

dd = rbindlist(dataPop)
dd$is_ena = as.double(!is.na(dd$age_ena))
dd2 = unique(dd[,c("ident","group","sex","is_ena")])

tab6 = dd2 |> 
  group_by(group) |> 
  summarise(prop_men = round(mean(sex=="M")*100,1), 
            prop_ena = round(mean(is_ena)*100,1)) %>%
  arrange(desc(prop_men))

kbl(cbind(tab6[1:5,],tab6[6:10,]), format = "latex",booktabs = T,linesep = "",
    caption = "Proportion of men and ENA graduates by organization.")|>
  save_kable(file = "documents/figures/descriptives_men_ena.tex")


#### Table : codex ####


df = readxl::read_xlsx("data/liste_corps/codex.xlsx")
df = df %>% select(abrégé, affectation, nomination, intégration, titularisation) %>% arrange(abrégé)
df[is.na(df)]=0
setnames(df, "abrégé","short_name")
kbl(df, format = "latex", booktabs=T, linesep = "",
    caption = "Type of traces used to define membership by organization") |>
  save_kable(file = "documents/figures/condex_sum.tex")
