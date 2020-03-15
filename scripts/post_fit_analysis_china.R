

## Check convergence
chain <- read.csv(output$file)
pdf("tmp.pdf")
plot(coda::as.mcmc(chain[,c("shape","scale","weibull_alpha","weibull_sigma","t0","K","lnlike")]))
dev.off()

chain <- chain[chain$sampno > mcmcPars["adaptive_period"],]
quants <- generate_prediction_intervals(chain, parTab, confirmed_data, 
                                        daily_import_probs = import_probs, daily_export_probs = export_probs,
                                        nsamp=100,return_draws = FALSE,model_ver=2)
#samps <- quants$samp_ids
#samps <- sort(samps)
chain_save <- chain[chain$sampno %in% samps,]
chain_save$sampno <- match(chain_save$sampno, samps)

#quants <- quants$draws

quants$province <- row.names(import_probs)[quants$province]
confirmed_data2 <- confirmed_data

confirmed_data2$province <- row.names(import_probs)[confirmed_data2$province]
confirmed_data2$date <- as.Date(confirmed_data2$date, origin="2019-11-01")
quants$date <- as.Date(quants$date, origin="2019-11-01")

#quants %>% filter(var %in% c("total_prevalence","infection_prev","onset_prev") & sampno == 1) %>% ggplot() +
#  geom_line(aes(x=date,y=n,col=var)) + facet_wrap(~province,scales="free_y")
#write_csv(quants, "prevalence_estimates_summary_logistic_growth.csv")
#write_csv(quants, "prevalence_estimates_logistic_growth.csv")
#write_csv(chain_save, "mcmc_chain_thinned_logistic_growth.csv")

hubei_plot <- ggplot(quants[quants$province == "Hubei" & #quants$date <= "2020-01-25" & 
                              quants$date >= "2019-12-01" &
                              quants$var %in% c("confirmations","onsets","total_prevalence"),]) + 
  geom_vline(xintercept=(as.Date("2020-01-23", format="%Y-%m-%d", tz="LMT", origin="2019-11-01")),linetype="dashed") +
  geom_line(aes(x=date,y=median,col=var)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  scale_x_date(breaks="2 days") +
  theme_bw() +
  ylab("Daily prevalence") +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position=c(0.2,0.2),
        panel.grid.minor=element_blank())

hubei_plot_inset <- ggplot(quants[quants$province == "Hubei" & quants$date <= "2020-01-01" & 
                                    quants$date >= "2019-12-01" &
                                    quants$var %in% c("confirmations","onsets"),]) + 
  geom_line(aes(x=date,y=median,col=var)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  scale_x_date(breaks="2 days") +
  scale_y_continuous(breaks=seq(0,160,by=20)) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none", 
        panel.grid.minor=element_blank())
vp <- viewport(width=0.5,height=0.5, x=0.4,y=0.7)
png("hubei_incidence.png", height=5,width=8,res=300,units="in")
print(hubei_plot)
print(hubei_plot_inset, vp=vp)
dev.off()

factor_levels <- confirmed_data2 %>% filter(province != "Hubei") %>% 
  group_by(province) %>%
  summarise(x=sum(n,na.rm=TRUE)) %>%
  arrange(-x) %>% pull(province)

quants <- quants[!(quants$province == "Hubei" & quants$var == "infections"),]
quants1 <- quants[!(quants$province == "Hubei" & quants$date > "2020-01-23"),]
quants1 <- quants1 %>% filter(province != "Hubei")
confirmed_data2 <- confirmed_data2 %>% filter(province != "Hubei")

confirmed_data2$province <- factor(confirmed_data2$province, levels=factor_levels)
quants1$province <- factor(quants1$province, levels=factor_levels)

p <- ggplot(quants1[quants1$var %in% c("infections","observations","onsets") & quants1$date >= "2020-01-01",]) +
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=var),alpha=0.25) +
  geom_line(aes(x=date,y=median,col=var)) +
  geom_point(data=confirmed_data2[confirmed_data2$date >= "2020-01-01",],aes(x=date,y=n),size=0.5) +
  scale_x_date(breaks="7 days") +
  geom_vline(xintercept=as.Date("2020-01-23",origin="2019-11-01"),linetype="dashed") +
  theme_pubr() + 
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~province,ncol=5,scales="free_y")
p


png("all_fits.png",height=10,width=12,units="in",res=300)
p + theme(legend.position="top")
dev.off()

library(data.table)
tmp <- chain[,colnames(chain) %like% "local_r"]
tmp <- tmp[,2:ncol(tmp)]
colnames(tmp) <- unique(confirmed_data2$province)
tmp1 <- reshape2::melt(tmp)
import_probs_melt <- reshape2::melt(import_probs)
import_probs_melt <- import_probs_melt %>% group_by(Var1) %>% summarise(val=mean(value)) %>% filter(Var1 != "Hubei") %>% ungroup()
colnames(import_probs_melt) <- c("variable", "Mean import probability")
tmp1 <- tmp1 %>% left_join(import_probs_melt)
tmp_order <- tmp1 %>% group_by(variable) %>% summarise(x=mean(value)) %>% arrange(-x) %>% pull(variable)
tmp1$variable <- factor(tmp1$variable, levels=factor_levels)
local_r_plot <- ggplot(tmp1) + 
  geom_violin(aes(x=variable, y=value, fill=`Mean import probability`),
              draw_quantiles=c(0.025,0.5,0.975),scale="width") + 
  scale_fill_gradient(low="blue",high="red") +
  geom_hline(yintercept=1,linetype="dashed") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="bottom") + 
  ylab("Average number of local cases per import") + xlab("Province") + 
  scale_y_continuous(expand=c(0,0),breaks=seq(0,10,by=1),limits=c(0,10))

png("local_r_plot.png",height=5,width=7,units="in",res=300)
local_r_plot
dev.off()

r0_ests <- read_csv("data/r0_ests.csv")

r0_ests <- r0_ests %>% select(Location,Value)
colnames(r0_ests) <- c("province","r0")

tmp2 <- plyr::ddply(tmp1, ~variable, function(x) mean(x$value))
colnames(tmp2) <- c("province","local_r")
wow <- merge(r0_ests, tmp2)

inc_p <- plot_incubation_period(chain,nsamp=200,xmax=40, prior_pars=inc_period_draws) + ylab("Probability density") + xlab("Delay from infection to symptom onset")
conf_p <- plot_confirm_delay(chain,nsamp=200,xmax=40) + ylab("Probability density") + xlab("Delay from symptom onset to confirmation")
delay_plots <- inc_p / conf_p
png("delay_plots.png",height=5,width=7,res=300,units="in")
delay_plots
dev.off()
