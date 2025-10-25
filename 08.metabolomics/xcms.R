############### (1) LC-MS data preprocessing and analysis with xcms (6)
############### (2) LC-MS/MS data analysis with xcms (560)
############### (3) Grouping FTICR-MS data with xcms
############### (4) LC-MS feature grouping

rm(list = lm())  # Clean up
############### (1) LC-MS data preprocessing and analysis with xcms

## ----biocstyle, echo = FALSE, results = "asis"--------------------------------
BiocStyle::markdown()

## ----init, message = FALSE, echo = FALSE, results = "hide"--------------------
## Silently loading all packages
library(BiocStyle)
library(xcms)
library(faahKO)
library(pander)
register(SerialParam())


## ----load-libs-pheno, message = FALSE-----------------------------------------
library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(pheatmap)
library(MsExperiment)

## Get the full path to the CDF files
cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
            recursive = TRUE)[c(1, 2, 5, 6, 7, 8, 11, 12)]
## Create a phenodata data.frame
pd <- data.frame(sample_name = sub(basename(cdfs), pattern = ".CDF",
                                   replacement = "", fixed = TRUE),
                 sample_group = c(rep("KO", 4), rep("WT", 4)),
                 stringsAsFactors = FALSE)

## -----------------------------------------------------------------------------
faahko <- readMsExperiment(spectraFiles = cdfs, sampleData = pd)
faahko

## -----------------------------------------------------------------------------
spectra(faahko)

## -----------------------------------------------------------------------------
table(fromFile(faahko))

## -----------------------------------------------------------------------------
sampleData(faahko)

## -----------------------------------------------------------------------------
faahko_3 <- faahko[3]
spectra(faahko_3)
sampleData(faahko_3)

## ----data-inspection-bpc, message = FALSE, fig.align = "center", fig.width = 12, fig.height = 6, fig.cap = "Base peak chromatogram."----
## Get the base peak chromatograms. This reads data from the files.
bpis <- chromatogram(faahko, aggregationFun = "max")
## Define colors for the two groups
group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
names(group_colors) <- c("KO", "WT")

## Plot all chromatograms.
plot(bpis, col = group_colors[sampleData(faahko)$sample_group])

## ----data-inspection-chromatogram, message = FALSE----------------------------
bpi_1 <- bpis[1, 1]
rtime(bpi_1) |> head()
intensity(bpi_1) |> head()

## -----------------------------------------------------------------------------
faahko <- filterRt(faahko, rt = c(2550, 4250))
## creating the BPC on the subsetted data
bpis <- chromatogram(faahko, aggregationFun = "max")

## ----data-inspection-tic-boxplot, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 4, fig.cap = "Distribution of total ion currents per file."----
## Get the total ion current by file
tc <- spectra(faahko) |>
    tic() |>
    split(f = fromFile(faahko))
boxplot(tc, col = group_colors[sampleData(faahko)$sample_group],
        ylab = "intensity", main = "Total ion current")

## ----data-inspection-bpc-heatmap, message = FALSE, fig.align = "center", fig.width = 7, fig.height = 6, fig.cap = "Grouping of samples based on similarity of their base peak chromatogram."----
## Bin the BPC
bpis_bin <- bin(bpis, binSize = 2)

## Calculate correlation on the log2 transformed base peak intensities
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- bpis_bin$sample_name

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = bpis_bin$sample_group)
rownames(ann) <- bpis_bin$sample_name

## Perform the cluster analysis
pheatmap(cormat, annotation = ann,
         annotation_color = list(group = group_colors))

## ----peak-detection-plot-eic, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 5, fig.cap = "Extracted ion chromatogram for one peak."----
## Define the rt and m/z range of the peak area
rtr <- c(2700, 2900)
mzr <- c(334.9, 335.1)
## extract the chromatogram
chr_raw <- chromatogram(faahko, mz = mzr, rt = rtr)
plot(chr_raw, col = group_colors[chr_raw$sample_group])

## ----peak-detection-plot-ms-data, message = FALSE, warning = FALSE, fig.aligh = "center", fig.width = 14, fig.height = 14, fig.cap = "Visualization of the raw MS data for one peak. For each plot: upper panel: chromatogram plotting the intensity values against the retention time, lower panel m/z against retention time plot. The individual data points are colored according to the intensity."----
faahko |>
    filterRt(rt = rtr) |>
    filterMz(mz = mzr) |>
    plot(type = "XIC")

## ----peak-detection-eic, message = FALSE--------------------------------------
xchr <- findChromPeaks(chr_raw, param = CentWaveParam(snthresh = 2))

## ----peak-detection-eic-chromPeaks--------------------------------------------
chromPeaks(xchr)

## ----peak-detection-chromatogram-chromPeakData--------------------------------
chromPeakData(xchr)

## ----peak-detection-eic-plot, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Signal for an example peak. Red and blue colors represent KO and wild type samples, respectively. Peak area of identified chromatographic peaks are highlighted in the sample group color."----
## Define a color for each sample
sample_colors <- group_colors[xchr$sample_group]
## Define the background color for each chromatographic peak
bg <- sample_colors[chromPeaks(xchr)[, "column"]]
## Parameter `col` defines the color of each sample/line, `peakBg` of each
## chromatographic peak.
plot(xchr, col = sample_colors, peakBg = bg)


## ----peak-detection-centwave, message = FALSE, results = "hide"---------------
cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000,
                     prefilter = c(6, 5000))
faahko <- findChromPeaks(faahko, param = cwp)

## ----peak-detection-chromPeaks, message = FALSE-------------------------------
chromPeaks(faahko) |>
    head()

## ----peak-detection-chromPeakData---------------------------------------------
chromPeakData(faahko)

## ----peak-postprocessing, message = FALSE-------------------------------------
mpp <- MergeNeighboringPeaksParam(expandRt = 4)
faahko_pp <- refineChromPeaks(faahko, mpp)

## ----peak-postprocessing-merged, fig.widht = 10, fig.height = 5, fig.cap = "Result from the peak refinement step. Left: data before processing, right: after refinement. The splitted peak was merged into one."----
mzr_1 <- 305.1 + c(-0.01, 0.01)
chr_1 <- chromatogram(faahko[1], mz = mzr_1)
chr_2 <- chromatogram(faahko_pp[1], mz = mzr_1)
par(mfrow = c(1, 2))
plot(chr_1)
plot(chr_2)

## ----peak-postprocessing-not-merged, fig.widht = 10, fig.height = 5, fig.cap = "Result from the peak refinement step. Left: data before processing, right: after refinement. The peaks were not merged."----
mzr_1 <- 496.2 + c(-0.01, 0.01)
chr_1 <- chromatogram(faahko[1], mz = mzr_1)
chr_2 <- chromatogram(faahko_pp[1], mz = mzr_1)
par(mfrow = c(1, 2))
plot(chr_1)
plot(chr_2)

## ----peak-postprocessing-chr, fig.width = 5, fig.height = 5-------------------
#' Too low minProp could cause merging of isomers!
res <- refineChromPeaks(chr_1, MergeNeighboringPeaksParam(minProp = 0.05))
chromPeaks(res)
plot(res)

## -----------------------------------------------------------------------------
faahko <- faahko_pp

## ----peak-detection-peaks-per-sample, message = FALSE, results = "asis"-------
summary_fun <- function(z)
    c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))

T <- chromPeaks(faahko) |>
    split.data.frame(f = chromPeaks(faahko)[, "sample"]) |>
    lapply(FUN = summary_fun) |>
    do.call(what = rbind)
rownames(T) <- basename(fileNames(faahko))
pandoc.table(
    T,
    caption = paste0("Summary statistics on identified chromatographic",
                     " peaks. Shown are number of identified peaks per",
                     " sample and widths/duration of chromatographic ",
                     "peaks."))

## ----peak-detection-chrom-peak-table-selected, message = FALSE----------------
chromPeaks(faahko, mz = c(334.9, 335.1), rt = c(2700, 2900))

## ----peak-detection-chrom-peaks-plot, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 8, fig.cap = "Identified chromatographic peaks in the m/z by retention time space for one sample."----
plotChromPeaks(faahko, file = 3)

## ----peak-detection-chrom-peak-image, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Frequency of identified chromatographic peaks along the retention time axis. The frequency is color coded with higher frequency being represented by yellow-white. Each line shows the peak frequency for one file."----
plotChromPeakImage(faahko, binSize = 10)

## ----peak-detection-eic-example-peak, message = FALSE-------------------------
chr_ex <- chromatogram(faahko, mz = mzr, rt = rtr)
chromPeaks(chr_ex)

## ----peak-detection-highlight-chrom-peaks-plot-polygon, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Signal for an example peak. Red and blue colors represent KO and wild type samples, respectively. The signal area of identified chromatographic peaks are filled with a color."----
sample_colors <- group_colors[chr_ex$sample_group]
plot(chr_ex, col = group_colors[chr_raw$sample_group], lwd = 2,
     peakBg = sample_colors[chromPeaks(chr_ex)[, "sample"]])

## ----peak-detection-highlight-chrom-peaks-plot, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Signal for an example peak. Red and blue colors represent KO and wild type samples, respectively. The rectangles indicate the identified chromatographic peaks per sample."----
plot(chr_ex, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_ex)[, "sample"]],
     peakBg = NA)

## ----peak-detection-chrom-peak-intensity-boxplot, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Peak intensity distribution per sample."----
## Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(faahko)[, "into"]),
              f = chromPeaks(faahko)[, "sample"])
boxplot(ints, varwidth = TRUE, col = sample_colors,
        ylab = expression(log[2]~intensity), main = "Peak intensities")
grid(nx = NA, ny = NULL)

## ----alignment-obiwarp, message = FALSE, results = "hide"---------------------
faahko <- adjustRtime(faahko, param = ObiwarpParam(binSize = 0.6))

## ----alignment-rtime, message = FALSE-----------------------------------------
## Extract adjusted retention times
adjustedRtime(faahko) |> head()

## Or simply use the rtime method
rtime(faahko) |> head()

## Get raw (unadjusted) retention times
rtime(faahko, adjusted = FALSE) |> head()

## ----alignment-obiwarp-plot, message = FALSE, fig.align = "center", fig.width = 12, fig.height = 8, fig.cap = "Obiwarp aligned data. Base peak chromatogram before (top) and after alignment (middle) and difference between adjusted and raw retention times along the retention time axis (bottom)."----
## Get the base peak chromatograms.
bpis_adj <- chromatogram(faahko, aggregationFun = "max", chromPeaks = "none")
par(mfrow = c(3, 1), mar = c(4.5, 4.2, 1, 0.5))
plot(bpis, col = sample_colors)
grid()
plot(bpis_adj, col = sample_colors)
grid()
## Plot also the difference of adjusted to raw retention time.
plotAdjustedRtime(faahko, col = sample_colors)
grid()

## ----alignment-peak-groups-example-peak, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 10, fig.cap = "Example extracted ion chromatogram before (top) and after alignment (bottom)."----
par(mfrow = c(2, 1))
## Plot the raw data
plot(chr_raw, col = sample_colors)
grid()
## Extract the chromatogram from the adjusted object
chr_adj <- chromatogram(faahko, rt = rtr, mz = mzr)
plot(chr_adj, col = sample_colors, peakType = "none")
grid()

## ----subset-define, message = FALSE, warning = FALSE--------------------------
faahko <- dropAdjustedRtime(faahko)

## Define the experimental layout
sampleData(faahko)$sample_type <- "study"
sampleData(faahko)$sample_type[c(1, 4, 7)] <- "QC"

## ----alignment-subset, message = FALSE, warning = FALSE-----------------------
## Initial peak grouping. Use sample_type as grouping variable
pdp_subs <- PeakDensityParam(sampleGroups = sampleData(faahko)$sample_type,
                             minFraction = 0.9)
faahko <- groupChromPeaks(faahko, param = pdp_subs)

## Define subset-alignment options and perform the alignment
pgp_subs <- PeakGroupsParam(
    minFraction = 0.85,
    subset = which(sampleData(faahko)$sample_type == "QC"),
    subsetAdjust = "average", span = 0.4)
faahko <- adjustRtime(faahko, param = pgp_subs)

## ----alignment-subset-plot-2, message = FALSE, warning = FALSE, fig.align = "center", fig.width = 10, fig.height = 10, fig.cap = "Subset-alignment results with option average. Difference between adjusted and raw retention times along the retention time axis. Samples on which the alignment models were estimated are shown in green, study samples in grey."----
clrs <- rep("#00000040", 8)
clrs[sampleData(faahko)$sample_type == "QC"] <- c("#00ce0080")
par(mfrow = c(2, 1), mar = c(4, 4.5, 1, 0.5))
plot(chromatogram(faahko, aggregationFun = "max", chromPeaks = "none"),
     col = clrs)
grid()
plotAdjustedRtime(faahko, col = clrs, peakGroupsPch = 1,
                  peakGroupsCol = "#00ce0040")
grid()

## ----correspondence-example-slice, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 10, fig.cap = "Example for peak density correspondence. Upper panel: chromatogram for an mz slice with multiple chromatographic peaks. lower panel: identified chromatographic peaks at their retention time (x-axis) and index within samples of the experiments (y-axis) for different values of the bw parameter. The black line represents the peak density estimate. Grouping of peaks (based on the provided settings) is indicated by grey rectangles."----
## Define the mz slice.
mzr <- c(305.05, 305.15)

## Extract and plot the chromatograms
chr_mzr <- chromatogram(faahko, mz = mzr)
## Define the parameters for the peak density method
pdp <- PeakDensityParam(sampleGroups = sampleData(faahko)$sample_group,
                        minFraction = 0.4, bw = 30)
plotChromPeakDensity(chr_mzr, col = sample_colors, param = pdp,
                     peakBg = sample_colors[chromPeaks(chr_mzr)[, "sample"]],
                     peakCol = sample_colors[chromPeaks(chr_mzr)[, "sample"]],
                     peakPch = 16)


## ----correspondence, message = FALSE------------------------------------------
## Perform the correspondence using fixed m/z bin sizes.
pdp <- PeakDensityParam(sampleGroups = sampleData(faahko)$sample_group,
                        minFraction = 0.4, bw = 30)
faahko <- groupChromPeaks(faahko, param = pdp)

## -----------------------------------------------------------------------------
## Drop feature definitions and re-perform the correspondence
## using m/z-relative bin sizes.
faahko_ppm <- groupChromPeaks(
    dropFeatureDefinitions(faahko),
    PeakDensityParam(sampleGroups = sampleData(faahko)$sample_group,
                     minFraction = 0.4, bw = 30, ppm = 10))

## ----fig.cap = "Relationship between a feature's m/z and the m/z width (max - min m/z) of the feature. Red points represent the results with the fixed m/z bin size, blue with the m/z-relative bin size."----
## Calculate m/z width of features
mzw <- featureDefinitions(faahko)$mzmax - featureDefinitions(faahko)$mzmin
mzw_ppm <- featureDefinitions(faahko_ppm)$mzmax -
    featureDefinitions(faahko_ppm)$mzmin
plot(featureDefinitions(faahko_ppm)$mzmed, mzw_ppm,
     xlab = "m/z", ylab = "m/z width", pch = 21,
     col = "#0000ff20", bg = "#0000ff10")
points(featureDefinitions(faahko)$mzmed, mzw, pch = 21,
       col = "#ff000020", bg = "#ff000010")


## -----------------------------------------------------------------------------
featureDefinitions(faahko) |> head()

## ----correspondence-feature-values, message = FALSE---------------------------
featureValues(faahko, value = "into") |> head()

## ----featureChromatograms, message = FALSE------------------------------------
feature_chroms <- featureChromatograms(faahko, features = 1:4)

feature_chroms

## ----feature-eic, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 8, fig.cap = "Extracted ion chromatograms for features 1 to 4."----
plot(feature_chroms, col = sample_colors,
     peakBg = sample_colors[chromPeaks(feature_chroms)[, "sample"]])


## -----------------------------------------------------------------------------
eic_2 <- feature_chroms[2, ]
chromPeaks(eic_2)

## ----fill-chrom-peaks, message = FALSE----------------------------------------
faahko <- fillChromPeaks(faahko, param = ChromPeakAreaParam())

featureValues(faahko, value = "into") |> head()

## ----export-result, eval = FALSE, echo = FALSE--------------------------------
#  save(faahko, file = "faahko.RData")

## -----------------------------------------------------------------------------
library(SummarizedExperiment)
res <- quantify(faahko, value = "into", method = "sum")
res

## -----------------------------------------------------------------------------
rowData(res)

## -----------------------------------------------------------------------------
colData(res)

## -----------------------------------------------------------------------------
assayNames(res)

## -----------------------------------------------------------------------------
assay(res) |> head()

## -----------------------------------------------------------------------------
assays(res)$raw_nofill <- featureValues(faahko, filled = FALSE, method = "sum")

## -----------------------------------------------------------------------------
assayNames(res)

## -----------------------------------------------------------------------------
assay(res, "raw_nofill") |> head()

## -----------------------------------------------------------------------------
metadata(res)

## -----------------------------------------------------------------------------
processHistory(faahko)[[1]]

## -----------------------------------------------------------------------------
processHistory(faahko)[[1]] |> processParam()

## ----final-pca, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 8, fig.cap = "PCA for the faahKO data set, un-normalized intensities."----
## Extract the features and log2 transform them
ft_ints <- log2(assay(res, "raw"))

## Perform the PCA omitting all features with an NA in any of the
## samples. Also, the intensities are mean centered.
pc <- prcomp(t(na.omit(ft_ints)), center = TRUE)

## Plot the PCA
pcSummary <- summary(pc)
plot(pc$x[, 1], pc$x[,2], pch = 21, main = "",
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                   digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                   digits = 3), " % variance"),
     col = "darkgrey", bg = sample_colors, cex = 2)
grid()
text(pc$x[, 1], pc$x[,2], labels = res$sample_name, col = "darkgrey",
     pos = 3, cex = 2)


## -----------------------------------------------------------------------------
# Set up parameters for RsdFilter
rsd_filter <- RsdFilter(threshold = 0.3,
                        qcIndex = sampleData(faahko)$sample_type == "QC")

# Apply the filter to faakho object
filtered_faahko <- filterFeatures(object = faahko, filter = rsd_filter)

# Now apply the same strategy to the res object
rsd_filter <- RsdFilter(threshold = 0.3, qcIndex = res$sample_type == "QC")
filtered_res <- filterFeatures(object = res, filter = rsd_filter, assay = "raw")

## -----------------------------------------------------------------------------
# Set up parameters for DratioFilter
dratio_filter <- DratioFilter(
    threshold = 0.5,
    qcIndex = sampleData(filtered_faahko)$sample_type == "QC",
    studyIndex = sampleData(filtered_faahko)$sample_type == "study")

# Apply the filter to faahko object
filtered_faakho <- filterFeatures(object = filtered_faahko,
                                  filter = dratio_filter)

# Now same but for the res object
dratio_filter <- DratioFilter(
    threshold = 0.5,
    qcIndex = filtered_res$sample_type == "QC",
    studyIndex = filtered_res$sample_type == "study")

filtered_res <- filterFeatures(object = filtered_res,
                               filter = dratio_filter)

## -----------------------------------------------------------------------------
# To set up parameter `f` to filter only based on QC samples
f <- sampleData(filtered_faakho)$sample_type
f[f != "QC"] <- NA

# To set up parameter `f` to filter per sample type excluding QC samples
f <- sampleData(filtered_faakho)$sample_type
f[f == "QC"] <- NA

missing_filter <- PercentMissingFilter(threshold = 30,
                                       f = f)

# Apply the filter to faakho object
filtered_faakho <- filterFeatures(object = filtered_faakho,
                                  filter = missing_filter)

# Apply the filter to res object
missing_filter <- PercentMissingFilter(threshold = 30,
                                       f = f)
filtered_res <- filterFeatures(object = filtered_res,
                               filter = missing_filter)

## -----------------------------------------------------------------------------
# Retrieve documentation for the main function and the specific filter.
?filterFeatures
?BlankFlag

## -----------------------------------------------------------------------------
ref <- loadXcmsData("xmse")
tst <- loadXcmsData("faahko_sub2")

## -----------------------------------------------------------------------------
f <- sampleData(ref)$sample_type
f[f != "QC"] <- NA
ref <- filterFeatures(ref, PercentMissingFilter(threshold = 0, f = f))
ref_mz_rt <- featureDefinitions(ref)[, c("mzmed","rtmed")]
head(ref_mz_rt)
nrow(ref_mz_rt)

## -----------------------------------------------------------------------------
param <- LamaParama(lamas = ref_mz_rt, method = "loess", span = 0.5,
                    outlierTolerance = 3, zeroWeight = 10, ppm = 20,
                    tolerance = 0, toleranceRt = 20, bs = "tp")

#' input into `adjustRtime()`
tst_adjusted <- adjustRtime(tst, param = param)
tst_adjusted <- applyAdjustedRtime(tst_adjusted)

## ----fig.height=12, fig.width=5-----------------------------------------------
#' evaluate the results with BPC
bpc <- chromatogram(ref, chromPeaks = "none")
bpc_tst_raw <- chromatogram(tst, chromPeaks = "none")
bpc_tst_adj <- chromatogram(tst_adjusted, chromPeaks = "none")

## ----fig.height=4, fig.width=10-----------------------------------------------
#' BPC of a sample
par(mfrow = c(1, 2),  mar = c(4, 2.5, 1, 0.5))
plot(bpc[1, 1], col = "#00000080", main = "Before Alignment")
points(rtime(bpc_tst_raw[1, 1]), intensity(bpc_tst_raw[1, 1]), type = "l",
       col = "#ff000080")
grid()

plot(bpc[1, 1], col = "#00000080", main = "After Alignment")
points(rtime(bpc_tst_adj[1, 1]), intensity(bpc_tst_adj[1, 1]), type = "l",
       col = "#0000ff80")
grid()

## -----------------------------------------------------------------------------
param <- matchLamasChromPeaks(tst, param = param)
mtch <- matchedRtimes(param)

#' BPC of the first sample with matches to lamas overlay
par(mfrow = c(1, 1))
plot(bpc[1, 1], col = "#00000080", main = "Distribution CP matched to Lamas")
points(rtime(bpc_tst_adj[1, 1]), intensity(bpc_tst_adj[1, 1]), type = "l",
       col = "#0000ff80")
grid()
abline(v = mtch[[1]]$obs)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 1))
plot(bpc[1, 2], col = "#00000080", main = "Distribution CP matched to Lamas")
points(rtime(bpc_tst_adj[1, 2]), intensity(bpc_tst_adj[1, 2]), type = "l",
       col = "#0000ff80")
grid()
abline(v = mtch[[2]]$obs)

## -----------------------------------------------------------------------------
#' access summary of matches and model information
summary <- summarizeLamaMatch(param)
summary

#' coverage for each file
summary$Matched_peaks / summary$Total_peaks * 100

#' access the information on the model of for the first file
summary$Model_summary[[1]]

#' Plot obs vs. ref with fitting line
plot(param, index = 1L, main = "ChromPeaks versus Lamas for the first file",
     colPoint = "red")
abline(0, 1, lty = 3, col = "grey")
grid()

## ----multicore, message = FALSE, eval = FALSE---------------------------------
#  register(bpstart(MulticoreParam(2)))

## ----snow, message = FALSE, eval = FALSE--------------------------------------
#  register(bpstart(SnowParam(2)))

## ----si-----------------------------------------------------------------------
sessionInfo()

############### (2) LC-MS/MS data analysis with xcms
############### (3) Grouping FTICR-MS data with xcms
############### (4) LC-MS feature grouping
