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

save(sf2Names, rppaNames, illuNames, illuType, file="../../data/cellLineNames.Rda")
