## Fig 2: prevalence and incidence in Hubei over time
quants <- read_csv("prevalence_estimates_summary_logistic_growth.csv")

quants_prev <- quants %>% filter(var %in% c("total_prevalence","onset_prev","infection_prev"))
var_key <- c("total_prevalence"="All infections","onset_prev"="Infected, post incubation period, not yet confirmed",
             "infection_prev"="Infected, in incubation period")
quants_prev$var <- var_key[quants_prev$var]
quants_prev$var <- factor(quants_prev$var, levels=c("All infections",
                                                    "Infected, in incubation period",
                                                    "Infected, post incubation period, not yet confirmed"))
colnames(quants_prev)[3] <- "Prevalence definition"
hubei_plot <- quants_prev %>% filter(province == "Hubei" & date >= "2019-12-08") %>%
  ggplot() + 
  geom_vline(xintercept=(as.Date("2020-01-23", format="%Y-%m-%d", tz="LMT", origin="2019-11-01")),linetype="dashed") +
  geom_line(aes(x=date,y=median,col=`Prevalence definition`)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=`Prevalence definition`),alpha=0.25) +
  scale_fill_manual(values=c("#E69F00","#CC79A7","#0072B2"))+
  scale_color_manual(values=c("#E69F00","#CC79A7","#0072B2"))+
  scale_x_date(breaks="7 days") + 
  theme_bw() +
  ylab("Daily prevalence (absolute numbers)") +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        axis.title = element_text(size=8),
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.position="bottom",
        panel.grid.minor=element_blank())
hubei_plot

hubei_plot_inset <- quants_prev %>% filter(province == "Hubei" & date >= "2019-12-08" & date <= "2020-01-01") %>%
  ggplot() + 
  geom_vline(xintercept=(as.Date("2020-01-23", format="%Y-%m-%d", tz="LMT", origin="2019-11-01")),linetype="dashed") +
  geom_line(aes(x=date,y=median,col=`Prevalence definition`)) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=`Prevalence definition`),alpha=0.25) +
  scale_fill_manual(values=c("#E69F00","#CC79A7","#0072B2"))+
  scale_color_manual(values=c("#E69F00","#CC79A7","#0072B2"))+
  scale_x_date(breaks="2 days") + 
  theme_bw() +
  ylab("") +
  xlab("") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=6),
        panel.background = element_blank(),
        axis.text.y=element_text(size=6),
        legend.position="none",
        panel.grid.minor=element_blank())
vp <- viewport(width=0.4,height=0.5, x=0.3,y=0.7)

pdf("figures/fig2A.pdf",height=4,width=8)
print(hubei_plot)
print(hubei_plot_inset, vp=vp)
dev.off()
