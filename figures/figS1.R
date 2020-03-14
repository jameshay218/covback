## Fig S1 - get incidence of case confirmations by province
confirmed_data <- read_csv("confirmed_data.csv")

confirmed_case_order <- confirmed_data %>% filter(country_region=="Mainland China") %>%
  group_by(province) %>% summarise(total_n = sum(n,na.rm=TRUE)) %>% arrange(-total_n) %>% pull(province)

confirmed_data <- confirmed_data %>% filter(country_region == "Mainland China")
confirmed_data$province <- factor(confirmed_data$province, levels=confirmed_case_order)
figS1 <- confirmed_data %>% filter(province != "Tibet") %>% 
  ggplot() + 
  geom_bar(aes(x=date,y=diff),stat="identity",
           fill="grey80",col="grey30",size=0.5) + 
  facet_wrap(~province, scales="free_y",ncol=5) +
  geom_vline(xintercept=(as.POSIXct("2020-01-23", format="%Y-%m-%d")),
             linetype="dashed",col="blue") +
  scale_x_datetime(breaks="5 days")+#,limits=as.POSIXct(c("2020-01-15","2020-03-03"),format="%Y-%m-%d")) +
  theme_bw() +
  ylab("Number of new case confirmations per day") +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        strip.background = element_rect(fill="grey90"),
        strip.text=element_text(size=8),
        legend.position=c(0.2,0.2),
        panel.grid.minor=element_blank())

pdf("figures/figS1.pdf",height=7,width=8)
figS1
dev.off()

figS1B <- confirmed_data %>% filter(province != "Tibet") %>% 
  ggplot() + 
  #geom_bar(aes(x=date,y=n),stat="identity",
  #         fill="grey80",col="grey30",size=0.5) + 
  geom_line(aes(x=date,y=n)) +
  facet_wrap(~province, scales="free_y",ncol=5) +
  geom_vline(xintercept=(as.POSIXct("2020-01-23", format="%Y-%m-%d")),
             linetype="dashed",col="blue") +
  scale_x_datetime(breaks="5 days")+#,limits=as.POSIXct(c("2020-01-15","2020-03-03"),format="%Y-%m-%d")) +
  theme_bw() +
  ylab("Cumulative number of confirmed cases") +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        strip.background = element_rect(fill="grey90"),
        strip.text=element_text(size=8),
        legend.position=c(0.2,0.2),
        panel.grid.minor=element_blank())

pdf("figures/figS1B.pdf",height=7,width=8)
figS1B
dev.off()
