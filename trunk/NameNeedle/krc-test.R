sf2Names  <- read.table("sf2Names.tsv",
                        sep="\t", header=TRUE, row.names=1, as.is=TRUE)
sf2Names <- sf2Names$x
rppaNames <- read.table("rppaNames.tsv",
                        sep="\t", header=TRUE, row.names=1, as.is=TRUE)
rppaNames <- rppaNames$x
illuNames <- read.table("IlluminaNames.tsv",
                        sep="\t", header=TRUE, row.names=1, as.is=TRUE)
illuType <- factor(illuNames$Type)
illuNames <- illuNames$illuNames
save(sf2Names, rppaNames, illuNames, illuType, file="cellLineNames.Rda")

# load the fast implementation of Needleman-Wunsch
source("R/needles.R")

needles("hcc123", "hcc-123")

myParams <- mainparams
myParams$MISMATCH <- -2
myParams$MATCH <- 2

# match RPPA names to SF2 names
matchscore <- matchcode <- rep(NA, length(sf2Names))
for (j in 1:length(sf2Names)) {
  scores <- sapply(rppaNames, function(x, y) {
    needles(x, y,myParams)$score
  }, y=sf2Names[j])
  matchcode[j] <- paste(which(scores==max(scores)), collapse=',')
  matchscore[j] <- max(scores)
}

rppaMatch <- sapply(matchcode, function(x) {
  y <- as.numeric(strsplit(x, ',')[[1]])
  paste(rppaNames[y], collapse="; ")
})

# match Illumina names to SF2 names
imatchscore <- imatchcode <- rep(NA, length(sf2Names))
for (j in 1:length(sf2Names)) {
  scores <- sapply(illuNames, function(x, y) {
    needles(x, y, myParams)$score
  }, y=sf2Names[j])
  imatchcode[j] <- paste(which(scores==max(scores)), collapse=',')
  imatchscore[j] <- max(scores)
}


illuMatch <- sapply(imatchcode, function(x) {
  y <- as.numeric(strsplit(x, ',')[[1]])
  paste(illuNames[y], collapse="; ")
})



matcher <- data.frame(rppaMatch=rppaMatch, rppaScore=matchscore,
                      illuMatch=illuMatch, illuScore=imatchscore,
                      combined)

write.table(matcher, file="namesMatched.tsv", sep="\t", quote=FALSE, col.names=NA)

