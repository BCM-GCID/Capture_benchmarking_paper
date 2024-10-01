library(tidyverse)
library(readxl)
library(ggpubr)
library(scales)
library(janitor)
library(RColorBrewer)
library(egg)
library(grid)

# Colors Pre-capture: #CA6B02; Capture: #2E068C

# Uploading table with coverage, ct and other information from each sample.

filtered_cov_tab = read_csv(file = "filt_cov_tab.csv")

# Making the plot for fig 4
color_mapping <- c("pre_capture" = "#CA6B02", "capture" = "#2E068C")

# Custom function to replace ct values of 45 with "ND" in the axis labels
label_custom <- function(x) {
  ifelse(x == 45, "ND", as.character(x))  # Replace 45 with ND
}

fig4= ggplot(data = filtered_cov_tab, 
                       aes(x = ct, y = perc_20x, color = treatment, group = SampleID, label = label, shape = treatment)) +
  geom_point(size = 3) +
  geom_line(linewidth=0.1, color="darkgrey") +  # Add lines to connect pairs
  labs(y = "Breadth of coverage\n(% of bases with 20x coverage)", x = "Ct", color = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  scale_color_manual(values = color_mapping, name = "Sequencing method",
                     labels = c("Capture", "Pre-capture")) +
  scale_shape_manual(values = c(17, 16), labels = c("Capture", "Pre-capture"),
                     guide ="none") +
  facet_grid(~genotype2, space="free_x")+
  theme(panel.spacing.x = unit(c(0.15, 1, 0.15, 0.15), "cm")) +
  guides(color = guide_legend(override.aes = list(shape = c(17, 16)))) +
  scale_x_continuous(labels = label_custom) 

ggsave(filename = "plots/breadth_of_coverage/newFig4.png",plot = fig4 , width=12, dpi = 600, bg="white")
ggsave(filename = "plots/breadth_of_coverage/newFig4.pdf",plot = fig4 , width=12, dpi = 600, bg="white")
ggsave(filename = "plots/breadth_of_coverage/newFig4.svg",plot = fig4 , width=12, dpi = 600, bg="white")
ggsave(filename = "plots/breadth_of_coverage/newFig4.tiff",plot = fig4 , width=12, dpi = 600, bg="white")