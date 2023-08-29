library("irr")
data <- read.csv("~/Concordanza_prova.csv", header=FALSE, sep=";" , row.names = 1)
icc(data, model = "twoway", type = "agreement", unit = "single")


data1 <- read.csv("Concordance_LETRS_tot.csv", header=FALSE, sep=";" , row.names = 1)
icc(data1, model = "twoway", type = "agreement", unit = "single")

data2 <- read.csv("Concordance_sgDE-tector_tot.csv", header=FALSE, sep=";" , row.names = 1)
icc(data2, model = "twoway", type = "agreement", unit = "single")

data3 <- read.csv("Concordance_Periscope_canonic.csv", header=FALSE, sep=";" , row.names = 1)
icc(data3, model = "twoway", type = "agreement", unit = "single")

data4 <- read.csv("Concordance_LETRS_canonic.csv", header=FALSE, sep=";" , row.names = 1)
icc(data4, model = "twoway", type = "agreement", unit = "single")

data5 <- read.csv("Concordance_sgDE-tector_canonic.csv", header=FALSE, sep=";" , row.names = 1)
icc(data5, model = "twoway", type = "agreement", unit = "single")

data6 <- read.csv("Concordance_Periscope_non_canonic.csv", header=FALSE, sep=";" , row.names = 1)
icc(data6, model = "twoway", type = "agreement", unit = "single")

data7 <- read.csv("Concordance_LETRS_non_canonic.csv", header=FALSE, sep=";" , row.names = 1)
icc(data7, model = "twoway", type = "agreement", unit = "single")

data8 <- read.csv("Concordance_sgDE-tector_non_canonic.csv", header=FALSE, sep=";" , row.names = 1)
icc(data8, model = "twoway", type = "agreement", unit = "single")
