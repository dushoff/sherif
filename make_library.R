
library(Rcpp)

# Input the name of the package here
mypackage <- "sherif"

cfiles <- c(
	grep(".h$", input_files, value=TRUE)
	, grep(".cpp$", input_files, value=TRUE)
)
cfolder <- "./"
cpp.path <- paste0(cfolder, cfiles)

# This generates all the necessary files 
# when creating an R package from scrath 
# with the view of interfacing with C++ (using Rcpp)
#?Rcpp::Rcpp.package.skeleton
Rcpp.package.skeleton(name =  mypackage,example_code = FALSE,
	author = "David Champredon",
	cpp_files = cpp.path, 
	force = TRUE
)
