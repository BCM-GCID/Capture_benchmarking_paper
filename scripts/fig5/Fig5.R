

library(ggplot2)
library(scales)
library("ggpubr")
library("reshape")



#RSV-A capture
data <- read.table ("./p1442-capture.RSV-A.allreadbam-s1.FPKM.txt", check.name = FALSE, header = T)
dat <-  as.data.frame(data)
melt.dat <-data.frame(reshape::melt(dat, id.vars = "gene"))
melt.dat$gene <- as.character(melt.dat$gene)
melt.dat$gene <- factor(melt.dat$gene,  levels = unique(melt.dat$gene))

pl = ggplot(melt.dat, aes(x = as.factor(gene), y = value)) +
  geom_boxplot(color="grey", fill="grey") +
  xlab("Genes") +
  ylab("Normalized read pair count (FPKM)" ) +
  theme(axis.title.x = element_text(size = 14, face="bold")) +
  theme(axis.title.y = element_text(size = 14, face="bold")) +
  theme(axis.text.x = element_text(hjust = 1, size = 36)) +
  theme(axis.text.y = element_text(size = 36)) +
  geom_jitter(shape=16, color="darkblue", size=2.5, position=position_jitter(0.1)) +
  scale_y_continuous(limits = c(0, 250000), labels = unit_format(unit = "k", scale = 1e-3)) +
  theme_bw()

ggsave("p1442-capture.RSV-A.allreadbam-s1.gene.FPKM.box.pdf", width = 12, height = 7)




#RSV-B capture
data <- read.table ("./p1442-capture.RSV-B.allreadbam-s1.FPKM.txt", check.name = FALSE, header = T)
dat <-  as.data.frame(data)
melt.dat <-data.frame(reshape::melt(dat, id.vars = "gene"))
melt.dat$gene <- as.character(melt.dat$gene)
melt.dat$gene <- factor(melt.dat$gene,  levels = unique(melt.dat$gene))

pl = ggplot(melt.dat, aes(x = as.factor(gene), y = value)) +
  geom_boxplot(color="grey", fill="grey") +
  xlab("Genes") +
  ylab("Normalized read pair count (FPKM)" ) +
  theme(axis.title.x = element_text(size = 14, face="bold")) +
  theme(axis.title.y = element_text(size = 14, face="bold")) +
  theme(axis.text.x = element_text(hjust = 1, size = 36)) +
  theme(axis.text.y = element_text(size = 36)) +
  geom_jitter(shape=16, color="darkblue", size=2.5, position=position_jitter(0.1)) +
  scale_y_continuous(limits = c(0, 250000), labels = unit_format(unit = "k", scale = 1e-3)) +
  theme_bw()

ggsave("p1442-capture.RSV-B.allreadbam-s1.gene.FPKM.box.pdf", width = 12, height = 7)




#RSV-A pre-capture
data <- read.table ("./p1442-precapture.RSV-A.allreadbam-s1.FPKM.txt", check.name = FALSE, header = T)
dat <-  as.data.frame(data)
melt.dat <-data.frame(reshape::melt(dat, id.vars = "gene"))
melt.dat$gene <- as.character(melt.dat$gene)
melt.dat$gene <- factor(melt.dat$gene,  levels = unique(melt.dat$gene))

pl = ggplot(melt.dat, aes(x = as.factor(gene), y = value)) +
  geom_boxplot(color="grey", fill="grey") +
  xlab("Genes") +
  ylab("Normalized read pair count (FPKM)" ) +
  theme(axis.title.x = element_text(size = 14, face="bold")) +
  theme(axis.title.y = element_text(size = 14, face="bold")) +
  theme(axis.text.x = element_text(hjust = 1, size = 36)) +
  theme(axis.text.y = element_text(size = 36)) +
  geom_jitter(shape=16, color="darkblue", size=2.5, position=position_jitter(0.1)) +
  scale_y_continuous(limits = c(0, 250000), labels = unit_format(unit = "k", scale = 1e-3)) +
  theme_bw()

ggsave("p1442-precapture.RSV-A.allreadbam-s1.gene.FPKM.box.pdf", width = 12, height = 7)



#RSV-B pre-capture RSV-B
data <- read.table ("../p1442-precapture.RSV-B.allreadbam-s1.FPKM.txt", check.name = FALSE, header = T)
dat <-  as.data.frame(data)
melt.dat <-data.frame(reshape::melt(dat, id.vars = "gene"))
melt.dat$gene <- as.character(melt.dat$gene)
melt.dat$gene <- factor(melt.dat$gene,  levels = unique(melt.dat$gene))

pl = ggplot(melt.dat, aes(x = as.factor(gene), y = value)) +
  geom_boxplot(color="grey", fill="grey") +
  xlab("Genes") +
  ylab("Normalized read pair count (FPKM)" ) +
  theme(axis.title.x = element_text(size = 14, face="bold")) +
  theme(axis.title.y = element_text(size = 14, face="bold")) +
  theme(axis.text.x = element_text(hjust = 1, size = 36)) +
  theme(axis.text.y = element_text(size = 36)) +
  geom_jitter(shape=16, color="darkblue", size=2.5, position=position_jitter(0.1)) +
  scale_y_continuous(limits = c(0, 250000), labels = unit_format(unit = "k", scale = 1e-3)) +
  theme_bw()

ggsave("p1442-precapture.RSV-B.allreadbam-s1.gene.FPKM.box.pdf", width = 12, height = 7)

