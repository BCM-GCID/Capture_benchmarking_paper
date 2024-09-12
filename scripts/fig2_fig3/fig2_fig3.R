library(tidyverse)
library(readxl)
library(scales)
library(RColorBrewer)
library(patchwork)


## Data input and wrangle ##
# Genome complete if: within length range, >90% completeness & >20x coverage
# Genome low coverage complete if: within length range, >90% completeness & <20x coverage
# Genome incomplete if: below length range, <90% completeness & <20x coverage
# NoV genome length range: 7.4-7.7kb
# RSV genome length range: 15-15.3 kb

input_raw <- readxl::read_excel("./benchmarking_captureEnrichment_rawInput.xlsx")

df <- input_raw %>%
  mutate(PercVirInput = ReadCount_target * 100 / Virmap_Input,
         Virus=factor(Virus, levels = c("RSV", "HuNoV")),
         comp = case_when(Virus == "HuNoV" & Length_genome >= 7400 & Length_genome <= 7750 
                          & Completeness_genome >= 90  & Coverage_genome >= 20  ~ 'Complete',
                          Virus == "HuNoV" & Length_genome >= 7400 & Length_genome <= 7750 
                          & Completeness_genome >= 90  & Coverage_genome < 20  ~ 'Complete-lowCov',
                          Virus == "HuNoV" & Length_genome < 7400 & Completeness_genome < 90  ~ 'Incomplete',
                          Virus == "HuNoV" & Length_genome < 7400 ~ 'Incomplete',
                          Virus == "HuNoV" & Completeness_genome < 90  ~ 'Incomplete',
                          Virus == "RSV" & Length_genome >= 15000 & Length_genome <= 15350 
                          & Completeness_genome >= 90  & Coverage_genome >= 20  ~ 'Complete',
                          Virus == "RSV" & Length_genome >= 15000 & Length_genome <= 15350 
                          & Completeness_genome >= 90  & Coverage_genome < 20  ~ 'Complete-lowCov',
                          Virus == "RSV" & Length_genome < 15000 & Completeness_genome < 90  ~ 'Incomplete',
                          Virus == "RSV" & Length_genome < 15000 ~ 'Incomplete',
                          Virus == "RSV" & Completeness_genome < 90  ~ 'Incomplete'),
         Genotype3 = case_when(Genotype2 == "GI.1" ~ 'GI.1',
                               Genotype2 == "GII.4" ~ 'GII.4',
                               Genotype2 == "GII.17" ~ 'Other GII',
                               Genotype2 == "GII.6" ~ 'Other GII',
                               Genotype2 == "GII.3" ~ 'Other GII',
                               Genotype2 == "RSV-A" ~ 'RSV-A',
                               Genotype2 == "RSV-B" ~ 'RSV-B'))

df$Virus <- factor(df$Virus, levels=c('RSV','HuNoV'))

#Add CT_Type designation.
df$CT_Type <- df$CT 
df[df$CT =='< LLQ', ]$CT <- '100'
df$CT <-round(as.double(df$CT),2)
df[df$CT == 100, ]$CT_Type <- 'ND'
df[df$CT < 20, ]$CT_Type <- 'CT < 20'
df[df$CT >= 20 & df$CT <= 30,]$CT_Type <- 'CT 20 to 30'
df[df$CT > 30 & df$CT != 100,]$CT_Type <- 'CT > 30'

## Figure 2
fig2 <- df %>% 
  ggplot(aes(x=Genotype3, y=PercVirInput)) +
  geom_boxplot(aes(x=Genotype3, y=PercVirInput), outlier.shape = NA, color = "darkgrey") +
  geom_jitter(aes(color = CT_Type, shape = Type), alpha=0.8, width = 0.2, size = 3.5) +
  facet_wrap(~ Virus, nrow = 1, scales = "free_x") +
  scale_color_manual(name = "CT Range",
                     breaks = c('CT < 20', 'CT 20 to 30', 'CT > 30', 'ND'),
                     values = c("#CC0033", "#56B4E9", "#339900", "#CC33FF")) +
  scale_shape_manual(name = "Sequencing method",
                     breaks = c("pre-capture", "capture"),
                     values=c(16, 17))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_text(size=11),
    axis.text.y  = element_text(size=11),
    axis.title = element_text(size=11),
    plot.title = element_text(size = 11)) +
  labs(y="Percent viral reads (%)", x="Expected Genotype") 

ggsave(fig2, file="./fig2.pdf", width = 11.69, height = 8.27, units = c("in"))



## Figure 3 
fig3 <- df %>% 
  ggplot() +
  geom_boxplot(aes(x=Genotype3, y=Coverage_genome), outlier.shape = NA, color = "darkgrey") +
  geom_jitter(aes(x=Genotype3, y=Coverage_genome, color = CT_Type, shape = Type), 
              alpha=0.8, width = 0.2, size = 3.5) +
  geom_hline(yintercept = 20, color = "black", linetype = "dotted") +
  
  facet_wrap(~Virus+comp, scales = "free_x", nrow = 1) +
  scale_y_log10(labels = label_number_si()) +
  scale_color_manual(name = "CT Range",
                      breaks = c('CT < 20', 'CT 20 to 30', 'CT > 30', 'ND'),
                      values = c("#CC0033", "#56B4E9", "#339900", "#CC33FF")) +
  scale_shape_manual(name = "Sequencing method",
                     breaks = c("pre-capture", "capture"),
                     values=c(16, 17))+
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_text(size=11),
    axis.text.y  = element_text(size=11),
    axis.title = element_text(size=11),
    plot.title = element_text(size = 11)) +
  labs(x="Genome completness", y="Average genome coverage (log10)")

ggsave(fig3, file="./fig3.pdf", width = 11.69, height = 8.27, units = c("in"))

