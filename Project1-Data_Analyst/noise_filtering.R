setwd("/projectnb/bf528/users/hedgehog/project_1/noise_filtering/")

# take in args instead of script defining table
args = commandArgs(trailingOnly=TRUE)

# assign args[1] to example file so I can be lazy
if (length(args) == 0){
  # file <- "example_intensity_data.csv"
  file <- "hedgehog1_data.csv"
} else {
  file <- args[1]
}

# split for writing filters out again
filesplit <- strsplit(file, ".", fixed = T)

# not actually a csv >.>, format is genes as rows, samples as cols
# samples <- read.table(file, header = T)
samples <- read.csv(file, header = T)

print(paste("Sample rows: ", dim(samples)[1], 
            "sample cols: ", dim(samples)[2]))
## filter1 ##
# the first filter, returns T/F if 20% or more of samples in x are greater
# than log2(15)
filter1 <- function(x, thresh) {
  return(sum(x > log(15, base = 2)) >= thresh)
}
# determine 20% of samples based on row width
min_samples <- dim(samples)[2] * .2
# apply to each row using apply(MARGIN = 1)
filt1_bool <- apply(samples[,-1], 1, function(x) filter1(x, min_samples))
print(paste("Filter 1 removed", dim(samples)[1] - sum(filt1_bool), "rows."))
filt1_applied <- samples[filt1_bool,]
f1_file <- paste0(filesplit[[1]][1], "-filter1.csv")
print(paste("Saving filter 1 to", f1_file, 
            "with dimensions:", paste(dim(filt1_applied), collapse = " x ")))
write.csv(filt1_applied, f1_file)


## filter 2 ##
# pretty sure probe sets is the list of genes, so determine median variance
# across all rows

probe_med_var <- median(apply(samples[,-1], 1, var))
probe_med_var2 <- median(apply(filt1_applied[,-1], 1, var))
# degrees of freedom <- N-1
dof <- dim(samples[,-1])[2] - 1
# calculate test statistics
test_stat <- function(x){
  return(dof*((sd(x)/probe_med_var)**2))
}
test_stat2 <- function(x){
  return(dof*(sd(x))/probe_med_var2)
}
# apply test_statistic formula row-wise, using apply()
 
filt2_bool <- apply(samples[,-1], 1, test_stat2) > qchisq(0.01, df = dof, lower.tail = F)
print(paste("Filter 2 removed", dim(samples)[1] - sum(filt2_bool), "rows."))
filt2_applied <- samples[filt2_bool,]
f2_file <- paste0(filesplit[[1]][1], "-filter2.csv")
print(paste("Saving filter 2 to", f2_file, 
            "with dimensions:", paste(dim(filt2_applied), collapse = " x ")))
write.csv(filt2_applied, f2_file)


## filter 3 ##
cv <- function(vector) {
  # determine coefficient of variation of a vector
  return(sd(vector)/mean(unlist(vector)))
}
# apply cv() to rows, return T/F of which are greater than 0.186
filt3_bool <- apply(samples[,-1], 1, cv) > 0.186
print(paste("Filter 3 removed", dim(samples[,-1])[1] - sum(filt3_bool), "rows."))
filt3_applied <- samples[filt3_bool,]
f3_file <- paste0(filesplit[[1]][1], "-filter3.csv")
print(paste("Saving filter 3 to", f3_file, 
            "with dimensions:", paste(dim(filt3_applied), collapse = " x ")))
write.csv(filt3_applied, f3_file)

## combined filters ##
# We can just use two AND gates to combine the three boolean vectors, and return
# a vector where TRUE means all three filters say TRUE
filtall_bool <- filt1_bool & filt2_bool & filt3_bool
print(paste("All filters removed", 
            dim(samples[,-1])[1] - sum(filtall_bool), 
            "rows."))
filtall_applied <- samples[filtall_bool,]
fall_file <- paste0(filesplit[[1]][1], "-filter_all.csv")
print(paste("Saving all filters combined to", fall_file, 
            "with dimensions:", paste(dim(filtall_applied), collapse = " x ")))
write.csv(filtall_applied, fall_file)

# determine how many of the original probes match our generated set
probesTru <- read.csv("probe_sets.csv")

count = 0
for (i in c(1, length(probesTru))) {
  if (probesTru[i,] %in% filtall_applied$X) {
    count <- count + 1
  }
}

recordProbes <- probesTru[apply(probesTru, MARGIN = 1, 
                                function(x){x %in% filtall_applied$X}),]
print(paste0(length(recordProbes), "/1459 probes match the original paper data set."))
write.csv(recordProbes, "probes_matching_paper.csv")
