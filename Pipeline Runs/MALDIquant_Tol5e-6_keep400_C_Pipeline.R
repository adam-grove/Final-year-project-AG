library("MALDIquant")
library("MALDIquantForeign")

df <- importMzMl("C:\\Users\\adamg\\OneDrive - The University of Liverpool\\Project Work\\BenchMarkingData\\Final Combined Dataset")

file_names <- list.files("C:\\Users\\adamg\\OneDrive - The University of Liverpool\\Project Work\\BenchMarkingData\\Final Combined Dataset")
file_names

# Looking for empties
any(sapply(df, isEmpty))

# Number of spectra 
table(sapply(df, length))

plot(df[[90]])
# Averaging the spectra 
samples <- factor(sapply(df,function(x)metaData(x)$file))
avgSpectra <- df
# avgSpectra <- averageMassSpectra(df, labels=samples,method="mean")
plot(avgSpectra[[1]])
# Estimating noise
noise <- estimateNoise(avgSpectra[[1]])
noise
plot(avgSpectra[[1]])
lines(noise, col="red")

# Detecting peaks 
peaks <- detectPeaks(avgSpectra, method="MAD",halfWindowSize=20, SNR=3.5)
plot(avgSpectra[[1]])
points(peaks[[1]], col="red", pch=4)
lines(noise[,1], noise[, 2]*3.5, col="blue")

peaks <- binPeaks(peaks, tolerance=2e-6)

featureMatrix <- intensityMatrix(peaks,avgSpectra)

library(data.table)

file_names <- sub("Combined_(.*)\\.mzML", "\\1", file_names)

keep <- NULL

#for(k in 1:ncol(featureMatrix)){
#  x <- length(which(featureMatrix[,k]>0))
#  if(x>400){
#    keep <- c(keep,k)
#  }
#}

n_chars_to_remove <- 3
library(stringr)
file_names <- str_sub(file_names, 1, nchar(file_names) - n_chars_to_remove)
rownames(featureMatrix) <- file_names


head(featureMatrix)

write.csv(featureMatrix,"MALDIquant_Tol5e-6_keep_C.csv") 

