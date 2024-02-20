devtools::install_github("jlmelville/coil20")
library(coil20)
library(R.matlab)
coil20 <- download_coil20(verbose = TRUE)

# object: 1 to 20 pose: 0-71
show_object(coil20, object = 1, pose = 0)

for (object in 1:20) {
  for (pose in 0:71) {
    r <- as.matrix(coil20[paste(object, pose, sep = "_"), 1:(ncol(coil20) - 1)])
    dim(r) <- c(128, 128)
    img <- png::readPNG(png::writePNG(r))
    #graphics::plot(c(100, 250), c(300, 450), type = "n", xlab = "", ylab = "",axes = FALSE)
    #graphics::rasterImage(img, 100,300, 250, 450, interpolate = FALSE)
    name <- paste("coil20_",paste(object, pose, sep = "_"),".mat",sep="")
    writeMat(name, coil20=img)
  }
}

# object: 1 to 100 pose: 0-71
options(download.file.method = "wininet")
coil100 <- download_coil100(verbose = TRUE)

for (object in 1:100) {
  for (pose in 0:71) {
    pose <- 5*pose
    r <- as.matrix(coil100[paste(object, pose, sep = "_"), 1:(ncol(coil100) - 1)])
    dim(r) <- c(128, 128, 3)
    img <- png::readPNG(png::writePNG(r))
    name <- paste("coil100_",paste(object, pose, sep = "_"),".mat",sep="")
    writeMat(name, coil100=img)
  }
}