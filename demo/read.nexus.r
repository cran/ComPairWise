oldwd <- getwd()
setwd(file.path(.Library, "ComPairWise", "examples"))
nex.data <- read.nexus("sample.nex")
setwd(oldwd)
rm(oldwd)
nex.data
