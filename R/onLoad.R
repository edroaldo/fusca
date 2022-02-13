.onLoad <- function(libname, pkgname) {
  rJava::.jpackage(pkgname, lib.loc = libname)
  # Remove when building the package.
  rJava::.jaddClassPath(dir(file.path(getwd(), "inst/java"), full.names = TRUE))
}
