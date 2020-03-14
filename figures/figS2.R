## Real export probs
export_probs <- read.csv("data/export_probs.csv")
export_probs <- export_probs$prob_leaving

## Real import probs
import_probs <- read.csv("data/import_probs.csv",header = TRUE)
import_probs <- t(import_probs[,2:ncol(import_probs)])
import_probs <- rbind(rep(1, ncol(import_probs)),import_probs)
row.names(import_probs)[1] <- "Hubei"

tmin <- as.POSIXct("2019-11-01",format="%Y-%m-%d", tz="UTC")
tmax <- as.POSIXct("2020-03-03",format="%Y-%m-%d",tz="UTC")
times <- seq(tmin, tmax, by="1 day")
export_probs_plot <- tibble(y=export_probs,x=times)
export_probs_plot <- export_probs_plot %>% filter(times >= "2019-12-08")

pA <- ggplot(data=export_probs_plot) +
  geom_line(aes(x=x,y=y))+
  geom_vline(xintercept=(as.POSIXct("2020-01-23", format="%Y-%m-%d")),
             linetype="dashed",col="blue") +
  scale_x_datetime(breaks="5 days")+
  theme_bw() +
  scale_y_continuous(limits=c(0,0.015), breaks=seq(0,0.015,by=0.001)) +
  ylab("Daily probability of leaving Wuhan\n to another Chinese province") +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        strip.background = element_rect(fill="grey90"),
        strip.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position=c(0.2,0.2),
        panel.grid.minor=element_blank())

import_probs_melted <- reshape2::melt(import_probs)
import_probs_melted$Var2 <- times[import_probs_melted$Var2]
import_probs_melted <- import_probs_melted %>% 
  filter(Var2 >= "2019-12-08" & Var1 != "Hubei")
colnames(import_probs_melted)[1] <- "Province"

import_order <- import_probs_melted %>% group_by(Province) %>% summarise(total_prob=sum(value)) %>%
  arrange(-total_prob) %>% pull(Province)
import_probs_melted$Province <- factor(import_probs_melted$Province, levels=import_order)


pB <- ggplot(data=import_probs_melted) +
  geom_line(aes(x=Var2,y=value, col=Province))+
  geom_vline(xintercept=(as.POSIXct("2020-01-23", format="%Y-%m-%d")),
             linetype="dashed",col="blue") +
  scale_x_datetime(breaks="5 days")+
  scale_y_continuous(limits=c(0,0.3),breaks=seq(0,0.3,by=0.05)) +
  guides(col=guide_legend(ncol=3)) +
  theme_bw() +
  ylab("Proportion of cases exported\n from Wuhan received") +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
        axis.text.y=element_text(size=7),
        strip.background = element_rect(fill="grey90"),
        strip.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.position=c(0.8,0.5),
        legend.text = element_text(size=6),
        panel.grid.minor=element_blank())
pB

figS2 <- pA / pB


pdf("figures/figS2.pdf",height=8,width=8)
figS2
dev.off()
