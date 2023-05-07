
###################
# what about those that show significant mthylation diffs? do they show diff fst/allele freq divergence?
###################
# get fst

## calculate fst:

library(poolfstat)
library(scales) 
setwd("~/hudsonica_genomics/analysis/RNA")


dat <- read.table("filtered_variants.txt", header=TRUE, stringsAsFactors=FALSE, nrow=2)

pops <- colnames(dat)[11:ncol(dat)]

pdat <- popsync2pooldata(sync.file="variants.sync", poolsizes = rep(100, length(pops)),
    poolnames = pops, min.maf=0, min.rc=0,
     max.cov.per.pool = 1e+07)

fst <- computePairwiseFSTmatrix(pdat, method = "Anova",
  min.cov.per.pool = -1, max.cov.per.pool = 1e+07, min.maf = -1,
  output.snp.values = TRUE)

#global.fsts <- computeFST(pdat, method = "Anova", snp.index = NA)

# will return the max of the pair. so values < 0 will be turned to 0
fst_mat <- pmax(fst$PairwiseFSTmatrix,0)

library(tidyr)
library(dplyr)
library(ggplot2)
library(plotrix)
# correlation style plot- nice but not using.

#pdf("~/tonsa_genomics/figures/fst.corr.pdf",height=15, width=15)
#par(mfrow = c(1, 1), mar=c(3, 3, 1, 1), mgp=c(3, 1, 0), las=0)
#
#corrplot(fst_mat, method="color", type="upper",cl.lim=c(0,max(fst$PairwiseFSTmatrix)),
#    col=colorRampPalette(c("blue","white","firebrick3"))(200),
#    is.corr=FALSE, addCoef.col="black",
#    tl.col="black", cl.pos ="n",
#    diag=FALSE, tl.cex = 1,tl.srt=45,
#    number.digits = 2)

#dev.off()

# turn to df
dat <- as.data.frame(as.table((fst_mat)))
colnames(dat) <- c("Samp1", "Samp2", "fst")
head(dat)

dat <- separate(dat, Samp1, into = c("Treatment_1", "Generation_1", "Rep_1"), sep = "_", remove=FALSE)
dat <- separate(dat, Samp2, into = c("Treatment_2", "Generation_2", "Rep_2"), sep = "_", remove=FALSE)

dat$trt_gen1 <- paste(dat$Treatment_1, dat$Generation_1, sep="_")
dat$trt_gen2 <- paste(dat$Treatment_2, dat$Generation_2, sep="_")

# drop self comps:
dat <- dat[dat$Samp1 != dat$Samp2,]

# find means and sd:
meanOut <- dat %>%
  group_by(trt_gen1, trt_gen2) %>%
  summarise(mean_value = mean(fst), sd_val = sd(fst))

write.table(file="~/hudsonica_genomics/analysis/RNA/mean_fst.txt", meanOut, 
            sep="\t", col.names=T, quote=F, row.names=F)

# write all out
dat2 <- dat[dat$Treatment_1 == dat$Treatment_2,]
dat3 <- dat2[dat2$Generation_1 != dat2$Generation_2,]
dat4 <- subset(dat3, !duplicated(fst)) # drop duplicate rows

write.table(file="~/hudsonica_genomics/analysis/RNA/group_fst.txt", dat4, 
            sep="\t", col.names=T, quote=F, row.names=F)
dat3 <- dat4
#plot fst across gens
factor_levels <- c("F02", "F04", "F11")
dat3$Generation_1 <- factor(dat3$Generation_1, levels = factor_levels)


p1 <- ggplot(dat3, aes(x=Generation_1, y=fst, fill=Treatment_1, shape=Treatment_1)) + 
        geom_point(color="black",position = position_dodge(width = 0.5), 
                    alpha=0.3, size = 3) +
        theme_bw() +
        scale_fill_manual(values=c('#6699CC',"#F2AD00",
                            "#00A08A", "#CC3333"),
                    labels = c("Ambient", "Acidification",
                               "Warming", "OWA"))+
        scale_shape_manual(values=c( 21,22, 23,24), guide = "none")+
        guides(fill=guide_legend(override.aes=list(
            shape=c(21,22, 23,24),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
            fill=c('#6699CC',"#F2AD00",
                            "#00A08A", "#CC3333")),order = 2)) +
        stat_summary(fun.data = "mean_cl_boot", 
                    position = position_dodge(width = 0.5),
                    size = 2,
                    show.legend=FALSE) +
        guides(color = guide_legend(
        override.aes = list(size = 2)
            )) +
        xlab("") +
        ylab("mean Fst vs. F0")

ggsave(file="~/hudsonica_genomics/analysis/RNA/pairwise_fst_vs_F0.pdf", plot=p1, h=4, w=6)

# get just the inter comp values

dat <- as.data.frame(as.table((fst_mat)))
colnames(dat) <- c("Samp1", "Samp2", "fst")
head(dat)

library(tidyr)
dat <- separate(dat, Samp1, into = c("Treatment_1", "Generation_1", "Rep_1"), sep = "_", remove=FALSE)
dat <- separate(dat, Samp2, into = c("Treatment_2", "Generation_2", "Rep_2"), sep = "_", remove=FALSE)

dat$trt_gen1 <- paste(dat$Treatment_1, dat$Generation_1, sep="_")
dat$trt_gen2 <- paste(dat$Treatment_2, dat$Generation_2, sep="_")


# drop self comps:
dat <- dat[dat$Samp1 != dat$Samp2,]
dat2 <- dat[dat$Treatment_1 != dat$Treatment_2,]
dat3 <- dat2[dat2$Generation_1 == dat2$Generation_2,]
dat4 <- subset(dat3, !duplicated(fst))

dat4$group <- paste(dat4$Treatment_2, dat4$Treatment_1, sep="_vs_")
dat4$group_gen <- paste(dat4$group, dat4$Generation_1, sep="_")

p1 <- ggplot(dat4, aes(x=group_gen, y=fst, fill=Treatment_1, shape=Treatment_1)) + 
        geom_point(color="black",position = position_dodge(width = 0.5), 
                    alpha=0.3, size = 3) +
        theme_bw() +
        scale_fill_manual(values=c("#F2AD00",
                            "#00A08A", "#CC3333"),
                    labels = c("Acidification",
                               "Warming", "OWA"))+
        scale_shape_manual(values=c(22, 23,24), guide = "none")+
        guides(fill=guide_legend(override.aes=list(
            shape=c(21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
            fill=c("#F2AD00",
                            "#00A08A", "#CC3333")),order = 2)) +
        stat_summary(fun.data = "mean_cl_boot", 
                    position = position_dodge(width = 0.5),
                    size = 2,
                    show.legend=FALSE) +
        guides(color = guide_legend(
        override.aes = list(size = 2)
            )) +
        xlab("") +
        ylab("mean Fst vs. F0")+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


ggsave(file="~/hudsonica_genomics/analysis/RNA/pairwise_fst_Group.pdf", plot=p1, h=4, w=8)


# parse down to only vs aa comparisons

dat5 <- dat4[which(dat4$Treatment_2 == "AA"),]

dat5$Treatment_1[dat5$Treatment_1 == "AH"] <- c("OA")
dat5$Treatment_1[dat5$Treatment_1 == "HA"] <- c("OW")
dat5$Treatment_1[dat5$Treatment_1 == "HH"] <- c("OWA")


write.table(file="~/hudsonica_genomics/analysis/RNA/fst_figure.txt", dat5, 
            sep="\t", col.names=T, quote=F, row.names=F)


meanOut <- dat5 %>%
  group_by(trt_gen1, trt_gen2) %>%
  summarise(mean_value = mean(fst), sd_err = std.error(fst))




p1 <- ggplot(dat5, aes(x=group_gen, y=fst, fill=Treatment_1, shape=Treatment_1)) + 
        geom_point(color="black",position = position_dodge(width = 0.5), 
                    alpha=0.3, size = 3) +
        theme_bw() +
        scale_fill_manual(values=c("#F2AD00",
                            "#00A08A", "#CC3333"),
                    labels = c("Acidification",
                               "Warming", "OWA"))+
        scale_shape_manual(values=c(22, 23,24), guide = "none")+
        guides(fill=guide_legend(override.aes=list(
            shape=c(21,22, 23),
        #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
            fill=c("#F2AD00",
                            "#00A08A", "#CC3333")),order = 2)) +
        stat_summary(fun.data = "mean_se", 
                    position = position_dodge(width = 0.5),
                    size = 1.6,
                    show.legend=FALSE) +
        guides(color = guide_legend(
        override.aes = list(size = 2)
            )) +
        xlab("") +
        ylab("mean Fst vs. F0")+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p1

ggsave(file="~/hudsonica_genomics/analysis/RNA/pairwise_fst_vs_AA.pdf", plot=p1, h=4, w=8)

# run anova
# just trying to ask, do we see differences as generations increase. 
# bc no F11 and also no F2 for JJ, makes stats weird. Just check if effect through time. and/or interacrtion
dat5$Generation_mod <-dat5$Generation_1
dat5$Generation_mod[which(dat5$Generation_mod == "F04" & 
                          dat5$Treatment_1 == "OWA")] <- c("F02")
dat5$Generation_mod[which(dat5$Generation_mod == "F11")] <- c("F04")
mod1 <- lm(fst~Treatment_1*Generation_mod, 
                            data=dat5)


m <- car::Anova(mod1, type = 3)
TukeyHSD(aov(mod1))

#



