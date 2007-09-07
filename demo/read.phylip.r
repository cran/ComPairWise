oldwd <- getwd()
setwd(file.path(.Library, "ComPairWise", "examples"))
phy.data <- read.phylip("sample.phy")
setwd(oldwd)
rm(oldwd)
phy.data

