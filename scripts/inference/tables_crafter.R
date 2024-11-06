library(data.table)
library(kableExtra)
library(ggplot2)
library(stringr)
library(ggplot2)
library(latex2exp)
library(ggrepel)
library(extrafont)
library("gridExtra")
library(splines)

"
Careful, obsburn parameter in config was renamed from burnin, as it is a more
appropriate name (we only burn the observations, not the parameters).
This was not updated for models trained before the R&R.
"

#### Study 1 ####

configs = expand.grid("model_index"=1:5,"linkedin_correct"=1,
                      "N"=1000, "obsburn"=900, "size_simu"=1,"iter_metropolis" = 30,
                      "data_augmentation"=0.01,"ena_pop"=T)

dt = list()
for (i in 1:nrow(configs)){
  print(i)
  dt[[i]] = readRDS(paste0("data/results/ena_config_",i,".rds"))
}
span = 201:5000

#### Table 1: Parameter posteriors for models 1-2 ####

# Importing 

span = 201:5000
pars = c("gamma0","gamma1","lambda")
dt2 = list()
for (i in c(1)){
  for (par in pars){
    dt2[[length(dt2)+1]] = data.table(est = quantile(dt[[i]]$params[[par]][,span], p = 0.05), 
                                      p = "q05", par = par, mod = dt[[i]]$config$model_index, coef = "(Intercept)", correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)+1]] = data.table(est = quantile(dt[[i]]$params[[par]][,span], p = 0.50), 
                                      p = "q50", par = par, mod = dt[[i]]$config$model_index, coef = "(Intercept)", correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)+1]] = data.table(est = quantile(dt[[i]]$params[[par]][,span], p = 0.95), 
                                      p = "q95", par = par, mod = dt[[i]]$config$model_index, coef = "(Intercept)", correct = dt[[i]]$config$linkedin_correct)
  }
}

for (i in c(2)){
  for (par in c("gamma0","gamma1")){
    dt2[[length(dt2)+1]] = data.table(est = quantile(dt[[i]]$params[[par]][,span], p = 0.05), 
                                      p = "q05", par = par, mod = dt[[i]]$config$model_index, coef = "(Intercept)", correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)+1]] = data.table(est = quantile(dt[[i]]$params[[par]][,span], p = 0.50), 
                                      p = "q50", par = par, mod = dt[[i]]$config$model_index, coef = "(Intercept)", correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)+1]] = data.table(est = quantile(dt[[i]]$params[[par]][,span], p = 0.95), 
                                      p = "q95", par = par, mod = dt[[i]]$config$model_index, coef = "(Intercept)", correct = dt[[i]]$config$linkedin_correct)
  }
}

for (i in c(2)){
par = "lambda"
dt2[[length(dt2)+1]] = data.frame(est = apply(dt[[i]]$params[[par]][,span], MARGIN = 1, FUN = quantile, p = 0.05), 
                                  p = "q05", par = par, mod = dt[[i]]$config$model_index, correct = dt[[i]]$config$linkedin_correct)
dt2[[length(dt2)]]$coef = row.names(dt2[[length(dt2)]]);rownames(dt2[[length(dt2)]]) = 1:nrow(dt2[[length(dt2)]])
dt2[[length(dt2)+1]] = data.frame(est = apply(dt[[i]]$params[[par]][,span], MARGIN = 1, FUN = quantile, p = 0.50), 
                                  p = "q50", par = par, mod = dt[[i]]$config$model_index, correct = dt[[i]]$config$linkedin_correct)
dt2[[length(dt2)]]$coef = row.names(dt2[[length(dt2)]]);rownames(dt2[[length(dt2)]]) = 1:nrow(dt2[[length(dt2)]])
dt2[[length(dt2)+1]] = data.frame(est = apply(dt[[i]]$params[[par]][,span], MARGIN = 1, FUN = quantile, p = 0.95), 
                                  p = "q95", par = par, mod = dt[[i]]$config$model_index, correct = dt[[i]]$config$linkedin_correct)
dt2[[length(dt2)]]$coef = row.names(dt2[[length(dt2)]]);rownames(dt2[[length(dt2)]]) = 1:nrow(dt2[[length(dt2)]])
}
for (i in c(3:5)){
  for (par in pars){
    dt2[[length(dt2)+1]] = data.frame(est = apply(dt[[i]]$params[[par]][,span], MARGIN = 1, FUN = quantile, p = 0.05), 
                                      p = "q05", par = par, mod = dt[[i]]$config$model_index, correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)]]$coef = row.names(dt2[[length(dt2)]]);rownames(dt2[[length(dt2)]]) = 1:nrow(dt2[[length(dt2)]])
    dt2[[length(dt2)+1]] = data.frame(est = apply(dt[[i]]$params[[par]][,span], MARGIN = 1, FUN = quantile, p = 0.50), 
                                      p = "q50", par = par, mod = dt[[i]]$config$model_index, correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)]]$coef = row.names(dt2[[length(dt2)]]);rownames(dt2[[length(dt2)]]) = 1:nrow(dt2[[length(dt2)]])
    dt2[[length(dt2)+1]] = data.frame(est = apply(dt[[i]]$params[[par]][,span], MARGIN = 1, FUN = quantile, p = 0.95), 
                                      p = "q95", par = par, mod = dt[[i]]$config$model_index, correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)]]$coef = row.names(dt2[[length(dt2)]]);rownames(dt2[[length(dt2)]]) = 1:nrow(dt2[[length(dt2)]])
  }
}

#  Renaming #

dt3 = rbindlist(dt2, use.name = T)

dt3[coef == "I(period - age_ena)",]$coef = "Cohort"
dt3[coef == "lagged1",]$coef = "$A(1)$"
dt3[coef == "lagged2",]$coef = "$A(0.8)$"
dt3[coef == "sqrt(lagged3)",]$coef = "$L$"
dt3[coef == "as.double(lagged1 == 0)",]$coef = "$I(A(1)=0)$"
dt3[coef == "I(lagged1 == 0)TRUE",]$coef = "$I(A(1)=0)$"
dt3[coef == "sqrt(lagged3):as.double(lagged1 == 0)",]$coef  = "$I(A(1)=0):L$"
dt3[coef == "as.double(round(period * 64) %in% c(10, 11, 24, 25, 34, 35, 44, 45, 54, 55, 64, 65))",]$coef = "Election"
dt3[coef == "as.double(round(period * 64) %in% c(9, 10, 23, 24, 33, 34, 43, 44, 53, 54, 63, 64))",]$coef = "Election"


dt3$coef = factor(dt3$coef, levels = c("(Intercept)","$L$","$A(1)$","$A(0.8)$","$I(A(1)=0)$","$I(A(1)=0):L$",
                                       "Cohort","Election"))

dt3$mod = paste0("$M_",dt3$mod,"$")

dt3[par=="gamma0",]$par = "$\\gamma_0$"
dt3[par=="gamma1",]$par = "$\\gamma_1$"
dt3[par=="lambda",]$par = "$\\lambda$"

dt3$est = round(dt3$est,3)

setnames(dt3,"par","parameter")
setnames(dt3, "p", "quantile")
setnames(dt3, "mod","model")

# Actual table

dt4 = dcast(dt3[correct==1 & !is.na(coef),], model + quantile ~ parameter+coef, value.var = "est",sep=", ")
dt4[is.na(dt4)] = ""

dt5 = t(dt4) |> as.data.frame()
rownames(dt5) = colnames(dt4)
colnames(dt5) = dt5[1,]
dt5 = dt5[2:nrow(dt5),]  

sel1 = !(str_detect(rownames(dt5), "Election")|str_detect(rownames(dt5), "Cohort"))
sel2 = str_detect(rownames(dt5), "Cohort")|str_detect(rownames(dt5), "Election")

kbl(dt5[sel1,1:6], format = "latex", booktabs=T, longtable = T, align = "r",
  row.names = T, escape = F, linesep = "", col.names = colnames(dt5)[1:6],
  caption = "Coefficients for Models 1 and 2.") |>
  kable_styling(font_size = 7) |>
  save_kable(file = "documents/figures/table_parameters1-1.tex") 

#### Figure 2: Marginal probability to leave (model 3) ####

y0 = bs((0:100)/100, df=5, intercept=T, degree=2,Boundary.knots = c(0,1)) %*% dt[[3]]$params[["gamma0"]][1:5,span]

#for (i in 1:101) y0[i,] = y0[i,]+dt[[3]]$params[["gamma0"]][1,span]

y1 = bs((0:100)/100, df=5, intercept=T, degree=2,Boundary.knots = c(0,1)) %*% dt[[3]]$params[["gamma1"]][1:5,span]
#for (i in 1:101) y1[i,] = y1[i,]+dt[[3]]$params[["gamma1"]][1,span]

res = rbind(data.frame(time = (0:100)/100*35, type = "gamma0",q50 = apply(y0, 1, median),
                       q05 = apply(y0, 1, quantile, probs=0.05), q95 = apply(y0, 1, quantile, probs=0.95)),
            data.frame(time = (0:100)/100*35, type = "gamma1", q50 = apply(y1, 1, median),
                       q05 = apply(y1, 1, quantile, probs=0.05), q95 = apply(y1, 1, quantile, probs=0.95))) |>
  as.data.table() |> melt(id.vars = c("time","type"))


ggplot(res[plogis(res$value)<=0.15,]) + geom_line(aes(x = time, y = plogis(value), linetype = variable)) + 
  scale_x_continuous("Years since admission") + 
  scale_y_continuous("Avg. prob. in 6 months") + 
  scale_linetype_manual("Quantile", values =c("q50"="solid","q05"="dashed","q95"="dashed")) +
  theme_bw() + facet_grid(type~.,scales="free", labeller = labeller) + guides(linetype = "none")

ggsave(units = "cm", width = 11, height = 6, filename = "documents/figures/time_prob.png")


#### Figure 3: Parameters of interest for models 4/5 ####

rr1 = readRDS("data/results/ena_config_4.rds")
rr2 = readRDS("data/results/ena_config_5.rds")

var1 = "I(period - age_ena)"
var2a = "as.double(round(period * 64) %in% c(9, 10, 23, 24, 33, 34, 43, 44, 53, 54, 63, 64))"
var2b = "as.double(round(period * 64) %in% c(10, 11, 24, 25, 34, 35, 44, 45, 54, 55, 64, 65))"
res1 = data.frame(
  param = c(rep("gamma0",length(span)),rep("gamma1",length(span)),rep("lambda",length(span))),
  value = c(rr1$params$gamma0[var1,span], rr1$params$gamma1[var1,span], rr1$params$lambda[var1,span]),
  mod = rep("Model 4: Cohort",3*length(span)))

res2 = data.frame(
  param = c(rep("gamma0",length(span)),rep("gamma1",length(span)),rep("lambda",length(span))),
  value = c(rr2$params$gamma0[var2a,span], rr2$params$gamma1[var2a,span], rr2$params$lambda[var2b,span]),
  mod = rep("Model 5: Election time",3*length(span)))

res = rbind(res1,res2) %>%
  filter(param != "lambda")

d2 <- res %>%
  group_by(mod, param) %>%
  summarize(lower = quantile(value, probs = .05),
            mid = quantile(value, probs = .50),
            upper = quantile(value, probs = .95))
#res$param[res$param=="gamma0"] = ("$\\gamma_0$")
#res$param[res$param=="gamma1"] = ("$\\gamma_1$")
table2 = ggplot(res) + 
  geom_histogram(aes(x = value),col = "black",fill="white") + 
  facet_grid(rows=vars(param),cols = vars(mod),scales="free") + theme_bw() +
  geom_vline(data=d2,aes(xintercept = lower), linetype="dashed") + 
  geom_vline(data=d2,aes(xintercept = upper), linetype = "dashed") +
  geom_vline(aes(xintercept=0),col = "black") + 
  theme(text = element_text(size = 9)) + 
  xlab("Parameter value") + ylab("Count after burn-in") 

ggsave(table2, 
       units = "cm", 
       width = 12, 
       height = 6, filename = "documents/figures/figure_models45.png",scale = 1)

#### Figure X: proportion 10 years after ####

fun_marginal_ena = function(x, span){
  
  y=merge(x$obs[it%in%span],x$given_data[,c("ident","age_max","age_ena","period")],
          by=c("ident","age_max"))[,-"age_max"][order(ident,it,age_ena),c("it","ident","traj","age_ena","period","avgA0")]
  y$cohort_ena = floor(((y$period - y$age_ena)*(NPER-1)+1)/2+1990)
  y$traj = as.double(y$traj)
  
  # we average within each iteration 
  rr2 = y[cohort_ena>=1990 & cohort_ena<=2010 & age_ena==20/(NPER-1), #20/(NPER-1)
          .(leaveO = mean(traj, na.rm = T), leaveE = 1-mean(avgA0, na.rm = T)),.(cohort_ena,it)]
  
  b = rr2[,.(
    leaveQ05 = quantile(leaveO,0.05, na.rm = F),
    leaveQ50 = quantile(leaveO,0.5, na.rm = F),
    leaveQ95 = quantile(leaveO,0.95, na.rm = F)),.(cohort_ena)]
  b$type = "Observation"
  
  return(b)
}

NPER = 65
res = fun_marginal_ena(dt[[3]],span = 4700:5000)

ggplot() + 
  geom_line(data = res[type == "Observation",], aes(x = cohort_ena, y = leaveQ50)) + 
  geom_errorbar(data = res, aes(x = cohort_ena, ymin = leaveQ05, ymax = leaveQ95, linetype = type), 
                position=position_dodge(width = 0.1), linetype = "dashed",
                linewidth = 0.5, width = 0.3) + 
  scale_y_continuous("Proportion") +
  scale_x_continuous("Year of admission") +
  theme_bw() + theme(text = element_text(size = 8),plot.margin=grid::unit(c(0,0,0,0), "mm"))  

ggsave(units = "cm", width = 13, height = 7, filename = "documents/figures/figure1.png",scale = 1)

#### Figure Rouban : average year for people leaving ####

NPER = 65
result = readRDS(paste0("data/results/ena_config_",3,".rds"))
y=merge(result$obs,result$given_data[,c("ident","age_max","age_ena","period")],
        by=c("ident","age_max"))[,-"age_max"][order(ident,it,age_ena),
                       c("it","ident","traj","age_ena","period","avgA0")]
y$cohort_ena = floor(((y$period - y$age_ena)*(NPER-1)+1)/2+1990)
y$traj = as.double(y$traj)

y2 = merge(y, y[,any(traj == 1 & age_ena <= 0.8),.(it,ident)][V1==T,c("it","ident")], on = c("it","ident"))[order(it,ident,age_ena)]

y3 = y2[,.(rle(traj)$lengths[1]/2),.(cohort_ena,it,ident)]

y3[,.(mean(V1)),.(cohort_ena)][order(V1) & cohort_ena,][order(cohort_ena)]

# 1997 et 1998 
table(y$cohort_ena)








#### Study 2 ####

configs = expand.grid("model_index"=1:3,"linkedin_correct"=1,
                      "N"=5000, "obsburn"=4700, "size_simu"=1,"iter_metropolis" = 30,
                      "data_augmentation"=0.1,"ena_pop"=F, "remove_early"=F)
dt = list()
for (i in 1:2){
  print(i)
  dt[[i]] = readRDS(paste0("data/results/groups_config_",i,".rds"))
}


#### Table X ####

# Importing # 
pars = c("gamma0","gamma1","lambda")
span = 21:50
dt2 = list()

for (i in c(1:4)){
  for (par in pars){
    dt2[[length(dt2)+1]] = data.frame(est = apply(dt[[i]]$params[[par]][,span], MARGIN = 1, FUN = quantile, p = 0.05), 
                                      p = "q05", par = par, mod = dt[[i]]$config$model_index, correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)]]$coef = row.names(dt2[[length(dt2)]]);rownames(dt2[[length(dt2)]]) = 1:nrow(dt2[[length(dt2)]])
    dt2[[length(dt2)+1]] = data.frame(est = apply(dt[[i]]$params[[par]][,span], MARGIN = 1, FUN = quantile, p = 0.50), 
                                      p = "q50", par = par, mod = dt[[i]]$config$model_index, correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)]]$coef = row.names(dt2[[length(dt2)]]);rownames(dt2[[length(dt2)]]) = 1:nrow(dt2[[length(dt2)]])
    dt2[[length(dt2)+1]] = data.frame(est = apply(dt[[i]]$params[[par]][,span], MARGIN = 1, FUN = quantile, p = 0.95), 
                                      p = "q95", par = par, mod = dt[[i]]$config$model_index, correct = dt[[i]]$config$linkedin_correct)
    dt2[[length(dt2)]]$coef = row.names(dt2[[length(dt2)]]);rownames(dt2[[length(dt2)]]) = 1:nrow(dt2[[length(dt2)]])
  }
}

# Renaming #

groups = c("cprefet","ccomptes","iga","ce","cdiplo","igf","dgfip","insee",
           "igas","dgtresor")
dt3 = rbindlist(dt2, use.name = T)
dt3$est = round(dt3$est,3)

#dt3[coef == "bs(age_max, df = 5, degree = 2, Boundary.knots = c(0, 1))1",]$coef = "$Age_1$"
#dt3[coef == "bs(age_max, df = 5, degree = 2, Boundary.knots = c(0, 1))2",]$coef = "$Age_2$"
#dt3[coef == "bs(age_max, df = 5, degree = 2, Boundary.knots = c(0, 1))3",]$coef = "$Age_3$"
#dt3[coef == "bs(age_max, df = 5, degree = 2, Boundary.knots = c(0, 1))4",]$coef = "$Age_4$"
#dt3[coef == "bs(age_max, df = 5, degree = 2, Boundary.knots = c(0, 1))5",]$coef = "$Age_5$"

#dt3[coef == "I(period - age_ena)",]$coef = "Cohort"
dt3[coef == "lagged1",]$coef = "$A(1)$"
dt3[coef == "lagged2",]$coef = "$A(0.8)$"
dt3[coef == "sqrt(lagged3)",]$coef = "$L$"
dt3[coef == "I(lagged1 == 0)TRUE",]$coef = "$I(A(1)=0)$"
dt3[coef == "sqrt(lagged3):as.double(lagged1 == 0)",]$coef  = "$I(A(1)=0):L$"


dt3[coef == "I(!is.na(age_ena))TRUE",]$coef = "ena-student"

org = groups[1]
for (org in groups){
  dt3[coef == paste0("group",org),]$coef = paste0(org,":(Intercept)")
  dt3[coef == paste0("group",org,":I(!is.na(age_ena))TRUE"),]$coef = paste0(org,":comparison") #ena-student
  dt3[coef == paste0("group",org,":sexF"),]$coef = paste0(org,":comparison") # women
  dt3[coef == paste0("group",org,":I(period - age_max)"),]$coef = paste0(org,":comparison") #cohort
}



dt3$coef = factor(dt3$coef, levels = 
                    c(paste0(groups,(":(Intercept)")),paste0(groups,(":comparison")),
                      paste0(groups,(":cohort")),"$L$","$A(1)$","$A(0.8)$","$I(A(1)=0)$","$I(A(1)=0):L$",
                      paste0("$Age_",1:5,"$")))


dt3$mod = paste0("$M_",dt3$mod,"$")

dt3[par=="gamma0",]$par = "$\\gamma_0$"
dt3[par=="gamma1",]$par = "$\\gamma_1$"
dt3[par=="lambda",]$par = "$\\lambda$"

setnames(dt3,"par","parameter")
setnames(dt3, "p", "quantile")
setnames(dt3, "mod","model")



dt4 = dcast(dt3[!is.na(coef) & correct==1 & parameter != "$\\lambda$",], model + quantile ~ parameter+coef, value.var = "est",sep=", ")
dt4[is.na(dt4)] = ""

dt5 = t(dt4) |> as.data.frame()
rownames(dt5) = colnames(dt4)
colnames(dt5) = dt5[1,]
dt5 = dt5[2:nrow(dt5),]

sel1 = !(str_detect(rownames(dt5),"cohort")|str_detect(rownames(dt5),"Age")|str_detect(rownames(dt5),"women"))
sel2 = !(str_detect(rownames(dt5), "ena-student"))

kbl(dt5[,4:9], format = "latex", digits = 4, booktabs = T, linesep = "",
    col.names = colnames(dt5)[4:9], row.names = T, escape = F, longtable = T,
    caption = "Coefficients for Models 2-4") |>
  kable_styling(font_size = 6) |>
  save_kable(file = "documents/figures/table_parameters2.tex")
  
#### Figure 5 Comparison across organizations ####

span = 201:5000
rr = data.table(group = substr(rownames(dt[[1]]$params$gamma0),6,nchar(rownames(dt[[1]]$params$gamma0))),
                
                gamma0_05 = plogis(apply(dt[[1]]$params$gamma0[,span], 1, quantile, p = 0.05)),
                gamma0_50 = plogis(apply(dt[[1]]$params$gamma0[,span], 1, median)),
                gamma0_95 = plogis(apply(dt[[1]]$params$gamma0[,span], 1, quantile, p = 0.95)),
                
                gamma1_05 = plogis(apply(dt[[1]]$params$gamma1[,span], 1, quantile, p = 0.05)),
                gamma1_50 = plogis(apply(dt[[1]]$params$gamma1[,span], 1, median)),
                gamma1_95 = plogis(apply(dt[[1]]$params$gamma1[,span], 1, quantile, p = 0.95)))

groups = c("ccomptes" ,"cdiplo", "ce", "cprefet", "dgtresor", "iga", "igas","igf", "insee")

ggplot(rr[group%in%groups,]) + 
  geom_point(aes(x = gamma0_50, y = gamma1_50)) +
  geom_text_repel(aes(x = gamma0_50, y = gamma1_50, label = group)) +
  geom_errorbar(aes(xmin = gamma0_05, xmax = gamma0_95, y = gamma1_50), linewidth = 0.25) + 
  geom_errorbar(aes(ymin = gamma1_05, ymax = gamma1_95, x = gamma0_50), linewidth = 0.25) +
  
  scale_y_continuous(TeX("Probability of coming back $\\gamma_1$")) + scale_x_continuous(TeX("Probability of leaving $\\gamma_0$")) + 
  
  theme_bw() + theme(text = element_text(size = 10),plot.margin=grid::unit(c(0,0,0,0), "mm"))  
  

ggsave(units = "cm", width = 8, height = 8, filename = "documents/figures/model1_coeffs.png",scale = 1)

#### Figure 6 : Comparison HF & ENA ####

# que faire pour dgfip? 
groups = c("ccomptes","ce","iga","igas","igf","dgtresor","cprefet","cdiplo","insee")
par = 'gamma0'
result = 0

r2a = list()
for (org in groups) {
  p1 = plogis(dt[[2]]$params[[par]][paste0("group",org,":sexF"),201:1000]+dt[[2]]$params[[par]][paste0("group",org),201:1000])
  p2 = plogis(dt[[2]]$params[[par]][paste0("group",org),201:1000])
  r2a[[org]] = data.table(org = org, p = p1*(1-p2)/(p2*(1-p1)))
}
r2 = rbindlist(r2a)
r2$org = factor(r2$org, levels=r2[,.(mean(p)),.(org)][order(V1)]$org)


r3a = list()
for (org in setdiff(groups,"insee")) {
  p1 = plogis(dt[[3]]$params[[par]][paste0("group",org,":I(!is.na(age_ena))TRUE"),201:1000]+dt[[3]]$params[[par]][paste0("group",org),201:1000])
  p2 = plogis(dt[[3]]$params[[par]][paste0("group",org),201:1000])
  r3a[[org]] = data.table(org = org, p = p1*(1-p2)/(p2*(1-p1)))
}
r3 = rbindlist(r3a)
r3$org = factor(r3$org, levels=setdiff(groups,"insee"))

plot2 = ggplot(r2,aes(p)) + 
  facet_wrap(.~org, ncol = 5) + theme_bw() +
  geom_histogram(bins = 30, col = "black", fill = "white") +
  geom_vline(xintercept = 1, col = "black") +
  scale_x_continuous("",limits = c(0.25,1.75)) + scale_y_continuous("") + 
  theme(text = element_text(size = 9))     

ggsave(plot = plot2, filename = "documents/figures/model2_coeffs_RR.png",
       width = 18, height = 8, units = "cm")


#### Les énarques démissionnaires ? ####

x = rbindlist(dataTot)[!is.na(age_ena),] 
setkey(x, "ident")
x = x[x[, .(has_true = any(age_ena==0.0625 & trace_obs==0)), by = ident][has_true==T,ident]]
x = x[x[, .(has_true = any(age_ena>0 & trace_obs == 1)), by = "ident"][has_true==F,ident]]
x$ident = as.integer(x$ident)


View(population_corps[ident%in%x$ident])
x=merge(x,population_corps, by ="ident")
write_csv(unique(population_corps[ident%in%x$ident,c("nom_prenom_std","ident","sexe","source_date")])[order(source_date),],file="names_to_check.csv")


#### Figure Rouban 2 : Passage cabinet ena ####

X = rbindlist(dataPop)
X$tet = X$group == "ena" & lead(X$group) == "cab"
X = X[ident == lead(ident),]
summary(glm(tet ~ age_ena+period,X,family=binomial("logit")))


summary(glm(tet ~ bs(age_ena,df=5,degree=2)+period,X,family=binomial("logit")))

#### ESS ####

readRDS(paste0("data/results/ena_config_",2,".rds"))
dt=readRDS(paste0("data/results/groups_config_",2,".rds"))
x = acf(dt$params$gamma0[21,200:5000])
4800 /(1+2*sum(x$acf[2:length(x$acf)]))



