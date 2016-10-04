readData <- function(filename, header = TRUE) {
    data <- read.csv(paste0("Data/",filename),header)
}