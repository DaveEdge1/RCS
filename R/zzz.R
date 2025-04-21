.onLoad = function (libname, pkgname) {
  ns = topenv()
  ns$TN_rwl = read.rwl(system.file("extdata", "TreeNobAllLumped10-7.csv", package = "RCS"))
  ns$TN_po = read.csv(system.file("extdata", "TN_POLumped_Oct_7_2020.csv", package = "RCS"))
}
