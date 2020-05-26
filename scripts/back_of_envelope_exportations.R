new_model_func <- function(pars, t, tmax=100){
  
  # 1. How much growth is there in Wuhan?
  r <- pars["growth_rate"]
  K <- pars["K"]
  t0 <- pars["t0"]
  tmax <- 100
  
  I_t <- numeric(tmax)
  t <- seq(0, tmax-t0, by=1)
  
  I_t_lead <- K / ( 1 + ( K - 1 ) * exp( -r * (t+1) ) )
  I_t_lag <- K / ( 1 + ( K - 1 ) * exp ( -r * t) )
  
  I_t[t0:tmax] <- I_t_lead - I_t_lag
  I_t
  
}

pars <- c("growth_rate"=0.3,"K"=200000,"t0"=37,
          "wuhan_pop"=11080000, "total_travellers"=4000000,
          "incu_par1" = 1.62,"incu_par2"=0.418,
          "prob_report"=0.5)


tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax,by="1 day")
duration <- as.numeric(tmax - tmin)

solve <- function(pars){
  # a) Pre-compute the incubation period
  incu_par1 <- pars["incu_par1"]
  incu_par2 <- pars["incu_par2"]
  
  ## Probability of remaining pre-symptomatic t days after infection
  presymptom_probs <- calculate_probs_presymptomatic_lnorm(duration, incu_par1, incu_par2)
  ## Probability of symptom onset t days after infection
  onset_probs <- calculate_onset_probs_lnorm(duration, incu_par1, incu_par2)
  
  # 1. Logistic growth of cumulative incidence in Wuhan
  r <- pars["growth_rate"]
  K <- pars["K"]
  t0 <- pars["t0"]
  prob_report <- pars["prob_report"]
  times[t0]
  I_wuhan <- daily_sigmoid_interval_cpp(r, K, duration, t0)
  
  # 2. Figure out how many of these infected individuals are going to leave
  wuhan_pop_ini <- pars["wuhan_pop"]
  total_travellers <- pars["total_travellers"]
  tmp <- create_export_prob_matrix(total_travellers, wuhan_pop_ini, export_prob_dat, tmin, tmax,
                                            index_date_end=index_date_end)
  ## This is the daily probability of someone leaving Wuhan
  export_probs <- tmp$probs
  
  ## Given infection on day i (row), what's the probability of leaving on each day in the future, j (col)?
  leave_matrix <- prob_leave_on_day(export_probs, duration)
  daily_prob_leaving <- prob_leave_pre_symptoms_vector(leave_matrix, presymptom_probs)
  cases_left <- daily_prob_leaving*I_wuhan*prob_report
  
  ret <- list(signif(sum(I_wuhan)/wuhan_pop_ini,3),
              sum(floor(cases_left)),
              t0+1,
              which.max(I_wuhan))
  ret
}
res <- solve(pars)

growth_rates <- seq(0.05,0.5,by=0.01)
Ks <- seq(100000,3000000,by=100000)
total_travellers <- seq(5000000,5000000,by=500000)
start_times <- c(37)
report_rates <- seq(0.1,1,by=0.1)

all_pars <- expand.grid("r"=growth_rates,"K"=Ks,
                        "total_travellers"=total_travellers,"t0"=start_times,
                        "report_rates"=report_rates)

final_sizes <- numeric(nrow(all_pars))
exported_cases <- numeric(nrow(all_pars))
start_times <- numeric(nrow(all_pars))
peak_times <- numeric(nrow(all_pars))

for(i in 1:nrow(all_pars)){
  if(i %% 1000 == 0) print(i)
  pars1 <- pars
  pars1["prob_report"] <- all_pars[i,5]
  pars1["t0"] <- all_pars[i,4]
  pars1["growth_rate"] <- all_pars[i,1]
  pars1["K"] <- all_pars[i,2]
  pars1["total_travellers"] <- all_pars[i,3]
  res <- solve(pars1)
  
  final_sizes[i] <- res[[1]]
  exported_cases[i] <- res[[2]]
  start_times[i] <- res[[3]]
  peak_times[i] <- res[[4]]
}

final <- tibble(growth_rate=all_pars[,1], K=all_pars[,2], travellers=all_pars[,3], start_time=all_pars[,4],report_rate=all_pars[,5],
                final_sizes, exported_cases, start_times, peak_times) %>% 
  mutate(start_times=as.Date(times[start_times]), peak_times=as.Date(times[peak_times]))

works <- final %>% filter(peak_times >= as.Date("2020-01-18") & peak_times <= as.Date("2020-01-23") &
                   exported_cases < 8400*1.1 & exported_cases > 8400*0.9)
p1 <- ggplot(works) + geom_tile(aes(x=report_rate,y=final_sizes,fill=growth_rate)) + scale_fill_viridis_c(limits=c(0.25,0.35)) +
  theme_bw() +
  scale_x_continuous(limits=c(0.05,1.05),breaks=seq(0.1,1,by=0.1)) +
  scale_y_continuous(limits=c(0,0.3),breaks=seq(0,0.3,by=0.05)) +
  xlab("Ascertainment rate outside of Wuhan") +
  ylab("Final cumulative incidence in Wuhan")

png("fig1.png",width=7,height=5,res=300,units="in")
p1
dev.off()


## 12753 is number of confirmed cases outside of Wuhan
print(paste0("Total incidence: ", signif(sum(I_wuhan)/wuhan_pop_ini,3)))
print(paste0("Total exported and detected infections: ", sum(floor(cases_left))))
print(paste0("Seed date: ", times[t0+1]))
print(paste0("Peak infection incidence: ", times[which.max(I_wuhan)]))

zhang_dat <- read_csv("data/real/zhang_appendix_lancet.csv")

zhang_dat

table(zhang_dat$`exposure_Wuhan/Hubei`)
res <- zhang_dat %>% 
  mutate(Location = ifelse(Location %in% c("Shenzhen","Other cities in Guangdong"), "Guangdong", Location)) %>%
  group_by(Location, `exposure_Wuhan/Hubei`) %>% 
  count() %>% 
  rename(exposure=`exposure_Wuhan/Hubei`) %>% 
  pivot_wider(names_from=exposure,values_from = n) %>% 
  ungroup() %>%
  group_by(Location) %>%
  mutate(
    n = (`0`+`1`+`2`),
    n_in_hubei=`1` + `2`,
    n_in_wuhan = `1`,
    prop_from_wuhan=n_in_wuhan/n,
    lower_confint = ifelse(n > 1, prop.test(n_in_wuhan, n)$conf.int[1], 0),
    upper_confint = ifelse(n > 1, prop.test(n_in_wuhan, n)$conf.int[2], 1),
    prop_from_hubei=n_in_hubei/n,
    lower_confint_hubei = ifelse(n > 1, prop.test(n_in_hubei, n)$conf.int[1], 0),
    upper_confint_hubei = ifelse(n > 1, prop.test(n_in_hubei, n)$conf.int[2], 1)
    ) %>%
  select(Location, prop_from_wuhan, lower_confint, upper_confint,
         prop_from_hubei, lower_confint_hubei, upper_confint_hubei)
res

p2 <- res %>% ggplot() + 
  geom_pointrange(aes(x=Location, y=prop_from_wuhan,ymin=lower_confint,ymax=upper_confint), col="blue") +
  geom_pointrange(aes(x=Location, y=prop_from_hubei,ymin=lower_confint_hubei,ymax=upper_confint_hubei), col="red") +
  ylab("Proportion of confirmed cases that were exposed in Wuhan") +
  theme_bw() +
  scale_y_continuous(limits=c(0,1)) + theme(axis.text.x=element_text(angle=45, hjust=1))

png("fig2.png",width=7,height=5,res=300,units="in")
p2
dev.off()


confirmed_data1 <- read_csv("data/real/midas_data_final.csv")
confirmed_data1 %>% group_by(province_raw) %>% summarise(n=sum(n)) %>% rename(Location=province_raw) %>%
  full_join(res) %>% mutate(n_imported=n*prop_from_wuhan) %>% summarise(sum(n_imported, na.rm=TRUE))

tmp <- zhang_dat %>% 
  group_by(Location) %>%
  mutate(report_date = as.Date(report_date, origin="2019/01/11"),
                     onset_date= as.Date(onset_date, oprigin="2019/01/11"),
                     report_delay = report_date - onset_date,
                     stage=onset_date < as.Date("2020/01/27")) %>% 
  select(ID, onset_date, report_date,report_delay, stage) %>% drop_na()

p3 <- tmp %>% group_by(onset_date) %>% 
  summarise(n=n(),mean_delay=mean(report_delay),sd_delay=sd(report_delay),
            lower_confint = mean_delay - 1.96*(sd_delay/sqrt(n)),
            upper_confint = mean_delay + 1.96*(sd_delay/sqrt(n))) %>% 
  filter(n > 10) %>%
  ggplot() +
  #geom_pointrange(aes(x=onset_date, ymin=lower_confint,ymax=upper_confint,y=mean_delay)) +
  scale_y_continuous(limits=c(0,25)) +
  geom_ribbon(aes(x=onset_date,ymin=lower_confint,ymax=upper_confint),alpha=0.25,fill="red") +
  theme_bw() +
  ylab("Delay from onset to report") +
  xlab("Date of onset") +
  geom_line(aes(x=onset_date,y=mean_delay), col="red") +
  geom_jitter(data=tmp, aes(x=onset_date,y=report_delay),size=0.1) +
  geom_smooth(data=tmp, aes(x=onset_date,y=report_delay),fill="blue") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(~Location)

png("fig3.png",width=8,height=6,res=300,units="in")
p3
dev.off()

tmp %>% ggplot() + geom_histogram(aes(x=report_delay)) + facet_wrap(~stage, scales="free_y")


