rm(list=ls())
library(dplyr)
library(tidyr)

# load calibaration data
calib <- read.table("./cruises_calib.txt",header = T)

# data cleaning and manipulation
calib$Year <- calib$Year-2009
recency <- calib %>% group_by(ID) %>% summarise(recency = max(Year))
calib <- cbind(calib, value = rep(1,nrow(calib)))
calib <- calib %>%
  spread(key="Year",value="value")
calib[is.na(calib)] <- 0
calib <- merge(cbind(calib, rowSums(calib[,3:7])),recency, by = "ID")
colnames(calib) <- c("customerID","2009","2010","2011","2012","2013","2014","frequency","recency")
NumTot <- nrow(calib)
calib_combine <- calib %>% count(frequency,recency)

################################## BB and BG-BB estimation
# BB model building and parameter estimating
log_bb <- function(pars, x, m) {
  a <- pars[1]
  b <- pars[2]
  lchoose(m, x) + lbeta(a+x, b+m-x) - lbeta(a,b)
}

LL_BB <-function(pars, x, m, n) {
  p <- exp(pars)
  LL_ind <- log_bb(p, x, m)
  return(-sum(LL_ind * n))
}

pairs <- c(0,0)

opt_BB <- optim(pairs, fn=LL_BB, x = calib_combine$frequency, m = 5, n=calib_combine$n)
params_bb <- exp(opt_BB[["par"]])
LL_bb <- -opt_BB[["value"]]
params_bb
LL_bb

# BG/BB model building and paramter estimating
library(BTYD)
LL_BGBB <- function(pars, x, tx, m, n) {
  p <- exp(pars)
  k <- pmax(0,lchoose(tx-1,x-1))
  LL_ind <- bgbb.LL(p, x, tx, m) + k
  return(-sum(LL_ind * n))
}

pairs <- c(0,0,0,0)
opt_BGBB <- optim(pairs, fn=LL_BGBB, x=calib_combine$frequency, tx=calib_combine$recency, m=5, n=calib_combine$n, control = list(maxit = 100000))
params_bgbb <- exp(opt_BGBB[["par"]])
LL_bgbb <- -opt_BGBB[["value"]]
params_bgbb
LL_bgbb

################################## Models fit 
# observation vs. prediction
nbr_calib <- calib %>% count(frequency)
colnames(nbr_calib) <- c("x", "N")
nbr_calib <- cbind(nbr_calib, Pred_N_BB = round(exp(log_bb(params_bb,0:5,5))*NumTot,0), 
             Pred_N_BGBB = round(bgbb.pmf(params_bgbb,5,0:5)*NumTot,0))
nbr_calib

library(ggplot2)
select(nbr_calib,x,Obs=N,Pred_N_BB, Pred_N_BGBB) %>%
  gather(Model,Count,c(Obs,Pred_N_BB, Pred_N_BGBB)) %>%
  ggplot(aes(x=x,y=Count,fill=Model,color=Model)) %>%
  + geom_bar(position="dodge",stat="identity") %>%
  + ggtitle("Calibration Data: Obs vs. Pred") %>%
  + scale_fill_manual(values=c('red','green','darkblue')) %>%
  + scale_color_manual(values=c('red','green','darkblue')) %>%
  + scale_y_continuous("Respondents") %>%
  + theme_bw() %>%
  + theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=15),
          legend.position="bottom")


dis_calib <- nbr_calib %>%
  mutate(Dis_BB=(nbr_calib$N-nbr_calib$Pred_N_BB)/nbr_calib$N,Dis_BGBB=(nbr_calib$N-nbr_calib$Pred_N_BGBB)/nbr_calib$N)

library(ggplot2)
select(dis_calib,x,Dis_BB,Dis_BGBB) %>%
  gather(Model,Count,c(Dis_BB,Dis_BGBB)) %>%
  ggplot(aes(x=x,y=Count,fill=Model,color=Model)) %>%
  + geom_line(stat="identity") %>%
  + ggtitle("Discrepancy Percentage") %>%
  + scale_fill_manual(values=c('green','darkblue')) %>%
  + scale_color_manual(values=c('green','darkblue')) %>%
  + scale_y_continuous("Respondents") %>%
  + theme_bw() %>%
  + theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=15),
          legend.position="bottom")

################################## let’s see how well the models capture the distribution of
################################## cruise counts during a holdout (forecast) period
# load holdout data
hold <- read.table("./cruises_hold.txt",header = T)

# data cleaning and manipulation
hold$Year <- hold$Year-2014
recency <- hold %>% group_by(ID) %>% summarise(recency = max(Year))
hold <- cbind(hold, value = rep(1,nrow(hold)))
hold <- hold %>%
  spread(key="Year",value="value")
id <- data.frame(ID = seq(length(calib$customerID)))
hold <- merge(id, hold, by = "ID", all.x = T)
hold[is.na(hold)] <- 0
hold <- merge(cbind(hold, rowSums(hold[,2:5])),recency, by = "ID", all.x = T)
hold[is.na(hold)] <- 0
colnames(hold) <- c("customerID","2015","2016","2017","2018","frequency","recency")
hold_agg <- hold %>% count(frequency,recency)

# observation vs. prediction
nbr_hold <- hold %>% count(frequency)
colnames(nbr_hold) <- c("x", "N")
nbr_hold <- cbind(nbr_hold, Pred_N_BB = round(exp(log_bb(params_bb,0:4,4))*NumTot,0), 
                  Pred_N_BGBB = round(bgbb.pmf.General(params_bgbb,5,4,0:4)*NumTot,0))
nbr_hold

select(nbr_hold,x,Obs=N,Pred_N_BB, Pred_N_BGBB) %>%
  gather(Model,Count,c(Obs,Pred_N_BB, Pred_N_BGBB)) %>%
  ggplot(aes(x=x,y=Count,fill=Model,color=Model)) %>%
  + geom_bar(position="dodge",stat="identity") %>%
  + ggtitle("Holdout Data: Obs vs. Pred") %>%
  + scale_fill_manual(values=c('red','green','darkblue')) %>%
  + scale_color_manual(values=c('red','green','darkblue')) %>%
  + scale_y_continuous("Respondents") %>%
  + theme_bw() %>%
  + theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=15),
          legend.position="bottom")

################################## Prediction and Analysis 
# merge calib and hold
df <- merge(calib[,1:7], hold[,1:5], by = "customerID")
df <- cbind(df, rowSums(df[,3:11]))
df <- df %>% mutate(df_recn = apply(df[2:11], 1, function(x) {
                   max(x * seq(0,9))
                 }))
#df <- cbind(df, df_recn)
colnames(df)[c(12,13)] <- c("frequency","recency")
df_agg <- df %>% count(frequency,recency)

# expected number of alive customers 
alive_exp <- round(sum(bgbb.PAlive(params_bgbb, df_agg$frequency, df_agg$recency, 9) * df_agg$n),0)
alive_exp


# total expected residual transactions
delta <- 0.13
DERT_exp <- sum(bgbb.DERT(params_bgbb, df_agg$frequency, df_agg$recency, 9, delta) * df_agg$n)
DERT_exp

# total expected residual lifetime value 
ERLV_exp <- DERT_exp * 241
ERLV_exp

##################################  how many cruises can the company expect from this new cohort in each year from 2020 to 2024?
# incremental transactions by year 
new_cust <- 11993
incre_trans <- round(diff(bgbb.Expectation(params_bgbb, 0:5) * new_cust),0)
year_cnt <- data.frame(year = c("2020","2021","2022","2023","2024"), incremental = incre_trans)
year_cnt
ggplot(year_cnt,aes(x = year,y = incremental,group = 1)) %>%
  + geom_line(color = "blue",size = 1) %>%
  + geom_point(color = "darkblue",size = 4) %>%
  + ggtitle("Incremental Transactions by Year") %>%
  + theme_bw() %>%
  + expand_limits(y = 700) %>%
  + geom_text(aes(label = incremental),hjust = 0.5,vjust = 2) %>%
  + theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size = 15),
          legend.position = "bottom")

# the maximum spend on the campaigns, in other words, the ECLV
max_camp <- (bgbb.DERT(params_bgbb, 0, 0, 0, delta) + 1) * new_cust*241
max_camp


