.datatable.aware <- TRUE

# In addition, it's a good practice to define global variables
# to avoid 'no visible binding for global variable' notes from R CMD check.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", ".I", ".N", ".SD", ".SDcols"))
}