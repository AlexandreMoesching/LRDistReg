## Code to prepare the weight for age dataset

# Load data
# devtools::install_github("warnes/SASxport")
growthdata <- SASxport::read.xport("NHANES III (1988-1994).xpt")

# Rename variables to lower-case
names(growthdata) <- tolower(names(growthdata))
names(growthdata)

# Rename some variables
names(growthdata)[1] <- "number"
names(growthdata)[5] <- "stat.wt"
names(growthdata)[6] <- "age.years"
names(growthdata)[7] <- "age.months"
names(growthdata)[9] <- "std.height"
names(growthdata)[10] <- "weight"
names(growthdata)[13] <- "bwt.r"
names(growthdata)[14] <- "l.bwt.flg.r"
names(growthdata)[15] <- "head.circ"
names(growthdata)[16] <- "rec.height"
names(growthdata)[17] <- "bwt.c"
names(growthdata)[18] <- "vl.bwt.flg.c"
names(growthdata)[19] <- "l.bwt.flg.c"
names(growthdata)[20] <- "dif.flg"
names(growthdata)[21] <- "vl.bwt.flg.r"
names(growthdata)

# Swap some columns
seq <- c(8, 6:7, 9, 16, 20, 10:11, 13:14, 21, 17, 19, 18, 15, 1:5, 12)
names(growthdata)[seq]
growthdata <- growthdata[,seq]
names(growthdata)

# Create dataset in the "data/" folder
usethis::use_data(growthdata, overwrite = TRUE)
