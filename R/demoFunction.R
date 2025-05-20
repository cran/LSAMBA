#' Get monolix demo project path.
#'
#' @return path to the monolix demo from 'Monolix' software
#' @export
#'
#' @examples
#' \dontrun{
#' print(getMLXdir())
#' }
getMLXdir <- function(){
  mlx.initializeLixoftConnectors()
  r <- mlx.getDemoPath()
  return(paste0(r,"/5.models_for_individual_parameters/5.2.covariate_model/warfarin_covariate1_project.mlxtran"))
}
