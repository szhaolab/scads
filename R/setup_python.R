#' Set up Python Environment for scads Package
#'
#' This function creates and configures a Python environment with the required packages.
#' @export
#' 
setup_python_env <- function(python_path=NULL) {
  
  # Check if python_path is null and Python 3 is available
  if(is.null(python_path)){
    python_path <- Sys.which("python3")
    if (python_path == "") {
      stop("Python3 is not installed or not found in system PATH.")
    }
  }
  # assignInNamespace("is_conda_python", function(x){return(FALSE)}, ns="reticulate") 
  # Create a virtual environment
  reticulate::virtualenv_create(envname = "scads_env", python = python_path)
  # Install required packages
  # reticulate::py_install(packages = c("setuptools"), envname = "scads_env", 
  #                        python_version = "/dartfs/rc/lab/S/Szhao/liyang/conda/base_clone/bin/python3")
  reticulate::virtualenv_install("scads_env", packages = c("pandas", "numpy", "scipy", "tqdm",
                                                           "pandas_plink"))
  # Inform the user
  message("Python environment 'scads_env' has been set up with required packages.")
  
}

