library(chipseqDBData)
tf.data <- NFYAData()
tf.data <- head(tf.data, -1) # skip the input.
bam.files <- tf.data$Path

cell.type <- sub("NF-YA ([^ ]+) .*", "\\1", tf.data$Description)
design <- model.matrix(~factor(cell.type))
colnames(design) <- c("intercept", "cell.type")

library(csaw)
param <- readParam(minq=20)
data <- windowCounts(bam.files, ext=110, width=10, param=param)

binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
keep <- filterWindowsGlobal(data, binned)$filter > log2(5)
data <- data[keep,]

data <- normFactors(binned, se.out=data)

library(edgeR)
y <- asDGEList(data)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)

merged <- mergeResults(data, results$table, tol=1000L)
