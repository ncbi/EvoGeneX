library("devtools")
pkg_path = '/home/pals2/Work/tissue-specificity/repou'

devtools::build(pkg_path)
devtools::install(pkg_path)
