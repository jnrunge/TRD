window <- 300 # for cleaning the allele frequencies (window = number of loci to look at at once)
percentile=0.7 # the middle (1-percentile) values are taken and then 1.5* of those upper and lower values taken

source("~/BrusselSprouts/scripts/functions.R")
library(naturalsort)
crosses = readLines("~/data/trd/mapped_reads/TRD.vcf.gz.samples")
crosses = crosses[startsWith(crosses, "YJNRC") | startsWith(crosses, "Chris")]
crosses
crosses_xlsx=readxl::read_xlsx("~/data/trd/Crosses.xlsx", sheet=2)
cc=data.frame("Cross ID"=c(paste0("ChrisC",1:8)),
                                  "Short name 1"=c("ACP","BAP","CCD","ATE","ACK","AKE","BAH","ANG"),
                                  "Short name 2"=c("BFP","CMP","CPG","SACE_YCR","ACV","BAH","CGD","CEI"), stringsAsFactors=FALSE)
colnames(cc)=str_replace_all(colnames(cc), fixed("."), " ")
crosses_xlsx <- bind_rows(crosses_xlsx[,c("Cross ID","Short name 1","Short name 2")],
                        cc)

s=as.numeric(args[1])

getAD <- function(x) {
  allele_count <- str_count(OH_vs_cross$ADcross[x], ",")
  AD3_0 <- as.numeric(strsplit(OH_vs_cross$ADcross[x], ",", fixed = TRUE)[[1]][1])
  AD3_1 <- as.numeric(strsplit(OH_vs_cross$ADcross[x], ",", fixed = TRUE)[[1]][2])
  if (allele_count > 1) {
    AD3_2 <- as.numeric(strsplit(OH_vs_cross$ADcross[x], ",", fixed = TRUE)[[1]][3])
    if (allele_count > 2) {
      AD3_3 <- as.numeric(strsplit(OH_vs_cross$ADcross[x], ",", fixed = TRUE)[[1]][4])
    } else {
      AD3_3 <- NA
    }
  } else {
    AD3_2 <- NA
    AD3_3 <- NA
  }
  return(data.frame(AD3_0 = AD3_0, AD3_1 = AD3_1, AD3_2 = AD3_2, AD3_3 = AD3_3))
}

isSNP <- function(x) {
  Alleles <- OH_vs_cross$alleles[x]
  Alleles <- strsplit(Alleles, ",", fixed = TRUE)[[1]]
  return(sum(str_length(Alleles) == 1) == length(Alleles))
}
PosMinus1Except1 <- function(x) {
  if (1 %in% x) {
    x[x == 1] <- x[x == 1] + 1
  }
  return(x - 1)
}


sample <- crosses[s]
print(sample)
hetLoci <- fread(paste0("~/data/trd/mapped_reads/", sample, ".hetLoci.gz"))
colnames(hetLoci) <- c("chr", "pos", "alleles", "ADcross")

hetLoci$DPcross <- NA
for (i in 1:nrow(hetLoci)) {
  hetLoci$DPcross[i] <- sum(as.numeric(strsplit(hetLoci$ADcross[i], ",", fixed = TRUE)[[1]]))
}

hetLoci <- subset(hetLoci, DPcross <= quantile(DPcross, 0.95))



OH1 <- fread(paste0("/home/jnrunge/data/trd/mapped_reads/", crosses_xlsx$`Short name 1`[crosses_xlsx$`Cross ID` == sample], ".homLoci.gz"))
OH2 <- fread(paste0("/home/jnrunge/data/trd/mapped_reads/", crosses_xlsx$`Short name 2`[crosses_xlsx$`Cross ID` == sample], ".homLoci.gz"))

colnames(OH1) <- c("chr", "pos", "alleles", "GT1")
colnames(OH2) <- c("chr", "pos", "alleles", "GT2")
OH <- full_join(OH1, OH2, by = c("chr", "pos", "alleles"))

OH <- subset(OH, !is.na(GT1) & !is.na(GT2))
OH <- subset(OH, substr(GT1, 1, 1) == substr(GT1, 3, 3))
OH <- subset(OH, substr(GT2, 1, 1) == substr(GT2, 3, 3))
OH <- subset(OH, substr(GT1, 1, 1) != substr(GT2, 1, 1))

OH_vs_cross <- full_join(OH, hetLoci, by = c("chr", "pos", "alleles"))

OH_vs_cross <- OH_vs_cross[naturalorder(OH_vs_cross$chr), ]

if (sum(duplicated(OH_vs_cross[, c(1, 2)])) > 0) {
  stop("Different alleles in OH and hetLoci tables for the same position")
}
OH_vs_cross <- subset(OH_vs_cross, !is.na(ADcross) & !is.na(GT1) & !is.na(GT2))
if (nrow(OH_vs_cross) == 0) {
  df_fit$countLociFit[i] <- 0
  next
}

OH_vs_cross <- OH_vs_cross[unlist(lapply(1:nrow(OH_vs_cross), isSNP)), ]
OH_vs_cross <- bind_cols(OH_vs_cross, bind_rows(lapply(1:nrow(OH_vs_cross), getAD)))
OH_vs_cross$Allele1 <- substr(OH_vs_cross$GT1, 1, 1)
OH_vs_cross$Allele2 <- substr(OH_vs_cross$GT2, 1, 1)


OH_vs_cross$AD_A1 <- NA
OH_vs_cross$AD_A2 <- NA

for (i in 1:nrow(OH_vs_cross)) {
  OH_vs_cross$AD_A1[i] <- OH_vs_cross[i, paste0("AD3_", OH_vs_cross$Allele1[i])]
  OH_vs_cross$AD_A2[i] <- OH_vs_cross[i, paste0("AD3_", OH_vs_cross$Allele2[i])]
}


chrs <- summarise(group_by(OH_vs_cross, chr), maxpos = max(pos))
chrs <- chrs[naturalorder(chrs$chr), ]

OH_vs_cross$global_pos <- OH_vs_cross$pos
for (c in 2:length(unique(OH_vs_cross$chr))) {
  chr <- unique(OH_vs_cross$chr)[c]
  OH_vs_cross$global_pos[OH_vs_cross$chr == chr] <- OH_vs_cross$pos[OH_vs_cross$chr == chr] + sum(chrs$maxpos[chrs$chr %in% unique(OH_vs_cross$chr)[1:(c - 1)]])
}


chrs$global_pos <- cumsum(chrs$maxpos)

OH_vs_cross$sumCount <- OH_vs_cross$AD_A1 + OH_vs_cross$AD_A2
OH_vs_cross <- subset(OH_vs_cross, sumCount > 0)
OH_vs_cross <- subset(OH_vs_cross, sumCount <= quantile(OH_vs_cross$sumCount, 0.95) & sumCount >= quantile(OH_vs_cross$sumCount, 0.05))

OH_vs_cross <- subset(OH_vs_cross, (sumCount / DPcross) > 0.99)

clean_data <- function(from, to, ovc) {
    ovc_sub<-slice(ovc, from:to)
  AFs <- ovc_sub %>%
    mutate(AF = AD_A1 / sumCount) %>%
    select(AF) %>%
    pull()
  q1 <- quantile(AFs, 0+(percentile/2))
  q3 <- quantile(AFs, 1-(percentile/2))
  iqr <- q3 - q1

  # Calculate the lower and upper bounds
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr

  keep<- which(AFs >= lower & AFs <= upper)

    return(ovc_sub%>%slice(keep))
    
}

clean_data_wrap=function(x,windows,ovc){
    return(clean_data(windows$from[x],windows$to[x],ovc))
}

windows=data.table()
for(chr in unique(OH_vs_cross$chr)){
    print(chr)
    which_rows=which(OH_vs_cross$chr==chr)
    print(min(which_rows))
    print(max(which_rows))
    sequence_from=round(seq(from=min(which_rows), to=max(which_rows), length.out = round(length(which_rows)/window)))
    if(length(sequence_from)==1){
        sequence_from=c(min(which_rows),max(which_rows))
    }
    print(head(sequence_from))
    sequence_to=sequence_from[-1]-1
    sequence_to[length(sequence_to)]=sequence_to[length(sequence_to)]+1
    sequence_from=sequence_from[-length(sequence_from)]
    windows=bind_rows(windows, data.table(from=sequence_from,to=sequence_to))
}


OH_vs_cross_cleaned=bind_rows(lapply(1:nrow(windows), FUN=clean_data_wrap, ovc=OH_vs_cross, windows=windows))


p<-ggplot(OH_vs_cross, aes((AD_A1 / sumCount))) +
  geom_histogram()
ggsave(paste0("/home/jnrunge/data/trd/quick_plots/", sample,"-AF-histo.png"), p)

p<-ggplot(OH_vs_cross_cleaned, aes((AD_A1 / sumCount))) +
  geom_histogram()
ggsave(paste0("/home/jnrunge/data/trd/quick_plots/", sample,"-AF-histo-cleaned.png"), p)

library(scales)

p<-ggplot(OH_vs_cross, aes(global_pos, AD_A1 / sumCount)) +
  geom_point(alpha = 0.1, color = "grey") +
  scale_color_viridis_c(option = "A", limits = c(0, 0.5)) +
  ylim(c(0, 1)) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = chrs$global_pos) +
  theme_bw(16) +
  ylab("Allele Frequency") +
  xlab("POS") +
  theme(legend.position = "none") +
  # geom_hline(yintercept = c(0.4,0.6))+
  ggtitle(sample) +
  labs(alpha = "Coverage") +
  scale_x_continuous(labels = comma)

ggsave(paste0("/home/jnrunge/data/trd/quick_plots/", sample,"-transmission-dirty.png"), p, width=10, height=3)


OH_vs_cross <- OH_vs_cross_cleaned # from now on only use this


print(paste("Fraction AD_A1==sumCount:", sum(OH_vs_cross$AD_A1 == OH_vs_cross$sumCount) / nrow(OH_vs_cross)))
print(paste("Mean AF:", mean(OH_vs_cross$AD_A1 / OH_vs_cross$sumCount)))

chr_summary <- summarise(group_by(OH_vs_cross, chr), meanAF = mean(AD_A1 / sumCount))

loessMod50 <- loess((AD_A1 / sumCount) ~ global_pos, data = OH_vs_cross, span = 0.01)
split_pred <- round(seq(from = 1, to = nrow(OH_vs_cross), length.out = round(nrow(OH_vs_cross) / 10000)))

# large prediction data (w/ SE calc) crashes predict.loess, so we do it per 10000 loci

for (k in 1:length(split_pred)) {
  if (k == length(split_pred)) {
    next
  }
  sp1 <- split_pred[k]


  if (k == (length(split_pred) - 1)) {
    sp2 <- split_pred[k + 1]
  } else {
    sp2 <- split_pred[k + 1] - 1
  }

  tmp <- predict(loessMod50, se = TRUE, newdata = OH_vs_cross[sp1:sp2, ])

  if (k == 1) {
    smoothed50 <- tmp
  } else {
    smoothed50 <- bind_rows(smoothed50, tmp)
  }
}



T <- qt(p = 0.975, df = smoothed50$df)
lwr <- smoothed50$fit - T * smoothed50$se.fit
upr <- smoothed50$fit + T * smoothed50$se.fit

OH_vs_cross$smoothed <- smoothed50$fit
OH_vs_cross$lwr <- lwr
OH_vs_cross$upr <- upr

library(scales)

OH_vs_cross <- OH_vs_cross[order(OH_vs_cross$global_pos), ]

p<-ggplot(OH_vs_cross, aes(global_pos, AD_A1 / sumCount)) +
  geom_point(alpha = 0.1, color = "grey") +
  geom_line(mapping = aes(global_pos, smoothed, color = abs(0.5 - smoothed)), inherit.aes = FALSE, linewidth = 2) +
  scale_color_viridis_c(option = "A", limits = c(0, 0.5)) +
  ylim(c(0, 1)) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = chrs$global_pos) +
  theme_bw(16) +
  ylab("Allele Frequency") +
  xlab("POS") +
  theme(legend.position = "none") +
  # geom_hline(yintercept = c(0.4,0.6))+
  ggtitle(sample) +
  labs(alpha = "Coverage") +
  scale_x_continuous(labels = comma)

ggsave(paste0("/home/jnrunge/data/trd/quick_plots/", sample,"-transmission.png"), p, width=10, height=3)

fwrite(OH_vs_cross, paste0("~/data/TRD/results/shiny/", sample, "-AF.csv.gz"))

rle_trd <- rle(x = abs(0.5 - OH_vs_cross$smoothed[-1]) > 0.1 & sign(OH_vs_cross$lwr[-1]) == sign(OH_vs_cross$upr[-1]) & OH_vs_cross$chr[-nrow(OH_vs_cross)] == OH_vs_cross$chr[-1])
rle_trd

# this is for plotting and the boundaries are not precise if distortion starts at POS=1
trd_regions <- data.frame(
  ID = 1:sum(rle_trd$values == TRUE),
  lengthSNPs = rle_trd$lengths[rle_trd$values == TRUE],
  chr_start = OH_vs_cross$chr[cumsum(rle_trd$lengths)[PosMinus1Except1(which(rle_trd$values == TRUE))] + 1],
  chr_end = OH_vs_cross$chr[cumsum(rle_trd$lengths)[(which(rle_trd$values == TRUE))]],
  global_start = OH_vs_cross$global_pos[cumsum(rle_trd$lengths)[PosMinus1Except1(which(rle_trd$values == TRUE))] + 1],
  global_end = OH_vs_cross$global_pos[cumsum(rle_trd$lengths)[(which(rle_trd$values == TRUE))]]
)
trd_regions$lengthBp <- trd_regions$global_end - trd_regions$global_start
trd_regions <- subset(trd_regions, lengthSNPs >= 100 & lengthBp >= 50000)

print(trd_regions)

if (file.exists(paste0("~/data/TRD/results/shiny/", sample, "-TRD_regions.csv.gz"))) {
  file.remove(paste0("~/data/TRD/results/shiny/", sample, "-TRD_regions.csv.gz"))
}

if (nrow(trd_regions) == 0) {
  next
}
fwrite(trd_regions, paste0("~/data/TRD/results/shiny/", sample, "-TRD_regions.csv.gz"))
for (i in 1:nrow(trd_regions)) {
  p<-ggplot(subset(OH_vs_cross, global_pos >= (trd_regions$global_start[i] - 100000) &
    global_pos <= (trd_regions$global_end[i] + 100000)), aes(global_pos, AD_A1 / sumCount)) +
    geom_point(shape = 1, color = "grey") +
    theme_bw(16) +
    ylab("Allele Frequency") +
    xlab("POS") +
    geom_line(mapping = aes(global_pos, smoothed, color = abs(0.5 - smoothed)), inherit.aes = FALSE, linewidth = 2) +
    theme(legend.position = "none") +
    scale_color_viridis_c(option = "A", limits = c(0, 0.5)) +
    geom_vline(xintercept = chrs$global_pos[chrs$global_pos >= (trd_regions$global_start[i] - 100000) &
      chrs$global_pos <= (trd_regions$global_end[i] + 100000)]) +
    ggtitle(sample) +
    labs(alpha = "Coverage") +
    geom_hline(yintercept = 0.5) +
    ylim(c(0, 1))
    
    ggsave(paste0("/home/jnrunge/data/trd/quick_plots/", sample,"-TRD-",i,".png"), p, width=10, height=3)
}
