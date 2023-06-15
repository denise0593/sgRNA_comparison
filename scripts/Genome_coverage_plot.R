if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("karyoploteR")) BiocManager::install("karyoploteR") 
library("karyoploteR")

## THE SCRIPT CREATES THE COVERAGE PLOT ALONG THE GENOME BASED ON BAM FILES

covid.genome <- toGRanges(data.frame(chr="MN908947.3", start=0, end=29903))
covid.cytobands <- toGRanges("/karyotype_covid")
kp <- plotKaryotype(genome = covid.genome, cytobands = covid.cytobands)
kpAddCytobandLabels(kp, cex=0.9, force.all = TRUE, srt=90)

# path to bam files
bam1 <- "6264_1_T0_NT_non_trattate.bam"
bam2 <- "6264_2_T0_NT_non_trattate.bam"
bam3 <- "6265_1_24h_NT.bam"
bam4 <- "6265_2_24h_NT.bam"
bam5 <- "6266_1_48h_NT.bam"
bam6 <- "6266_2_48h_NT.bam"
bam7 <- "6267_1_72h_NT.bam"
bam8 <- "6267_2_72h_NT.bam"
bam9 <- "6268_1_96h_NT.bam"
bam10 <- "6268_2_96h_NT.bam"

kp <- kpPlotBAMCoverage(kp, data=bam1, col="gold2", r0=0, r1=0.07, ymax=34500)
kpAddLabels(kp, "6264_1_0h", label.margin = 0.05, r0=0, r1=0.1)
kpAxis(kp, r0=0, r1=0.1, ymax=34500, numticks = 2)

kp <- kpPlotBAMCoverage(kp, data=bam2, col="gold2", r0=0.13, r1=0.2, ymax=34500)
kpAddLabels(kp, "6264_2_0h", label.margin = 0.05, r0=0.11, r1=0.2)
kpAxis(kp, r0=0.13, r1=0.2, ymax=34500, numticks = 2)

kp <- kpPlotBAMCoverage(kp, data=bam3, col="deepskyblue3", r0=0.23, r1=0.3, ymax=34500)
kpAddLabels(kp, "6265_1_24h", label.margin = 0.05, r0=0.21, r1=0.3)
kpAxis(kp, r0=0.23, r1=0.3, ymax=34500, numticks = 2)

kp <- kpPlotBAMCoverage(kp, data=bam4, col="deepskyblue3", r0=0.33, r1=0.4, ymax=34500)
kpAddLabels(kp, "6265_2_24h", label.margin = 0.05, r0=0.31, r1=0.4)
kpAxis(kp, r0=0.33, r1=0.4, ymax=34500, numticks = 2)

kp <- kpPlotBAMCoverage(kp, data=bam5, col="paleturquoise3",  r0=0.43, r1=0.5, ymax=34500)
kpAddLabels(kp, "6266_1_48h", label.margin = 0.05, r0=0.41, r1=0.5)
kpAxis(kp, r0=0.43, r1=0.5, ymax=34500, numticks = 2)

kp <- kpPlotBAMCoverage(kp, data=bam6, col="paleturquoise3", r0=0.53, r1=0.6, ymax=34500)
kpAddLabels(kp, "6266_2_48h", label.margin = 0.05, r0=0.51, r1=0.6)
kpAxis(kp, r0=0.53, r1=0.6, ymax=34500, numticks = 2)

kp <- kpPlotBAMCoverage(kp, data=bam7, col="indianred3", r0=0.63, r1=0.7, ymax=34500)
kpAddLabels(kp, "6267_1_72h", label.margin = 0.05, r0=0.61, r1=0.7)
kpAxis(kp, r0=0.63, r1=0.7, ymax=34500, numticks = 2)

kp <- kpPlotBAMCoverage(kp, data=bam8, col="indianred3", r0=0.73, r1=0.8, ymax=34500)
kpAddLabels(kp, "6267_2_72h", label.margin = 0.05, r0=0.71, r1=0.8)
kpAxis(kp, r0=0.73, r1=0.8, ymax=34500, numticks = 2)

kp <- kpPlotBAMCoverage(kp, data=bam9, col="springgreen4", r0=0.83, r1=0.9, ymax=34500)
kpAddLabels(kp, "6268_1_96h", label.margin = 0.05, r0=0.81, r1=0.9)
kpAxis(kp, r0=0.83, r1=0.9, ymax=34500, numticks = 2)

kp <- kpPlotBAMCoverage(kp, data=bam10, col="springgreen4", r0=0.93, ymax=34500)
kpAddLabels(kp, "6268_2_96h", label.margin = 0.05, r0=0.9)
kpAxis(kp, r0=0.93, ymax=34500, numticks = 2)