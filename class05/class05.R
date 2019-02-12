#' ---
#' title: "class05"
#' author: "julien"
#' date: "20190124"
#' output: github_document
#' ---


#Class 05 R graphics intro

#My frirst boxplot

#?rnorm or help(rnorm) to learn about normal distribution
x <- rnorm(1000,0)
boxplot(x)
#'I have generated **bold** x **italics** and it has `r code(x)`

# summary(x)
# hist(x)
# ?boxplot
# boxplot(x, horizontal = TRUE)
# ?read.table

weight <- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)
plot(weight, type="o", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", main="Weight Chart from Text File")

features <- read.table("bimm143_05_rstats/feature_counts.txt", sep = "\t", header = TRUE)
features
barplot(features$Count, names.arg = features$Feature, xlab = "Feature", horiz = TRUE, las=1, mar=c(1,1,1,1), main = "Feature Counts from Text File")
#margin parameters need to be done by guess and check
par(mar=c(3.1, 12.1, 4.1, 2.1))
barplot(features$Count, names.arg = features$Feature, xlab = "Feature", horiz = TRUE, las=1, main = "Feature Counts from Text File", xlim = c(0,80000))

#2C: Histograms
hist_values <- c(rnorm(1000), rnorm(1000)+4)
par(mar=c(5,4,4,2)+0.1)
hist(hist_values, breaks=50, col = rainbow(50))

#3B Coloring by Value
updown <- read.table("bimm143_05_rstats/up_down_expression.txt", sep = "\t", header = TRUE)
# updown
table(updown$State)
par(mar=c(5,4,4,2)+0.1)
plot(updown$Condition1, updown$Condition2, col = updown$State, xlab = "Expression Condition 1", ylab = "Expression Condition 2", main = "R day eRRYday")
levels(updown$State)
palette(c("blue", "gray", "red"))


