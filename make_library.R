library(Rcpp)

# Input the name of the package here
mypackage <- "sherif"

# List of C++ files to be included in the package:
rootnames <- c("EventNumberContainer",
               "dcTools",
               "RV",
               "globalVar",
               "individual",
               "mc",
               "dcDataFrame",
               "simulator",
               "dcMatrix")
nf <- length(rootnames)
fullnames <- paste(rootnames,c(rep("h",nf),rep("cpp",nf)),sep=".")
cppfiles <- c(fullnames, "Rwrap_sherif.cpp")
cpp.folder <- "./"
cpp.path <- paste0(cpp.folder,cppfiles)

# This generates all the necessary files 
# when creating an R package from scrath 
# with the view of interfacing with C++ (using Rcpp)
#?Rcpp::Rcpp.package.skeleton
Rcpp.package.skeleton(name =  mypackage,
	author = "David Champredon",
	cpp_files = cpp.path, 
	force = TRUE
)
