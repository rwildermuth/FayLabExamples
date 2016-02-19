# example change for git repository branch

# Fay lab meeting to go over R applications in stats and plotting
# Date created: 8/11/2015

# Read in data files
# Can also use the window to right to Environment -> Import Dataset
goa01to05 <- read.csv("/Users/rwildermuth/DropBox/PhD UMass/R Stats Class/goa2001_2005/goa2001_2005.csv", 
                      stringsAsFactors=FALSE)

goa07to13 <- read.csv("/Users/rwildermuth/DropBox/PhD UMass/R Stats Class/goa2007_2013/goa2007_2013.csv", 
                      stringsAsFactors=FALSE)

# merge all data
allData <- rbind(goa01to05, goa07to13)
dim(allData) # 114678 x 17
names(allData)

# Try subsetting for 2007 data
data2007 <- allData[allData$YEAR == 2007, ]
nrow(data2007)
# Or create index vector:
pick <- which(allData$YEAR == 2007)
length(pick) # !!! One less

head(allData$DATETIME)

data2007 <- data2007[data2007$COMMON == "Pacific halibut", ]
dim(data2007)
summary(data2007)
str(data2007)

# Plot halibut data for 2007
par(mfrow=c(2,2), mar=c(2,2,2,2), oma=c(2,2,0,0))
plot(data2007$LATITUDE, data2007$WTCPUE)
plot(data2007$STRATUM, data2007$WTCPUE)
plot(data2007$BOT_TEMP, data2007$WTCPUE)
plot(data2007$BOT_DEPTH, data2007$WTCPUE)
hist(log(data2007$WTCPUE), breaks = 100)
class(data2007$BOT_DEPTH)
data2007$BOT_DEPTH <- as.numeric(data2007$BOT_DEPTH)
summary(data2007$BOT_DEPTH) # One NA created

# Start looking at full dataset
# Boxplot of halibut CPUE over time
boxplot(log(WTCPUE)~YEAR, data=allData[allData$COMMON == "Pacific halibut", ])
mtext("Year", side = 1, line = 2) # adds margin text
mtext("log(CPUE)", side = 2, line = 2)

# extract halibut data
datHalibut <- allData[allData$COMMON == "Pacific halibut", ]
summary(log(datHalibut$WTCPUE))
# Find for all years
aggregate(datHalibut$WTCPUE,
          by = list(Year = datHalibut$YEAR), 
          summary)

# Number of rows for each year value
table(datHalibut$YEAR)

# test out apply() statements
xx <- matrix(rnorm(100), ncol=10)
row_means <- apply(xx, MARGIN=1, FUN=mean)
rowMeans(xx)

# test linear regression for change in CPUE over time
test1 <- lm(formula = log(WTCPUE) ~ YEAR, data = datHalibut)
summary(test1)
plot(test1)

pairs(datHalibut[, c("LATITUDE", "LONGITUDE", "YEAR", "WTCPUE", "BOT_TEMP")])
?pairs
