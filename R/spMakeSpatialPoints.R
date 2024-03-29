#' Spatial - Generates an S4 SpatialPoints object from X/Y coordinates.
#' 
#' Generates an S4 SpatialPoints object with defined projection from a data
#' table or matrix including X and Y coordinates, with option to export as an
#' ArcGIS shapefile (*.shp).
#' 
#' 
#' @param xyplt Data frame object or String. Name of layer with xy coordinates
#' and unique identifier. Can be layer with xy_dsn, full pathname, including
#' extension, or file name (with extension) in xy_dsn folder.
#' @param xyplt_dsn String. Name of database or folder were xyplt is. The dsn
#' varies by driver. See gdal OGR vector formats
#' (https://www.gdal.org/ogr_formats.html).
#' @param xy.uniqueid String. Unique identifier of xyplt rows.
#' @param xvar String. Name of variable in xyplt defining x coordinate.
#' @param yvar String. Name of variable in xyplt defining y coordinate.
#' @param xy.crs PROJ.4 String or CRS object or Integer EPSG code defining
#' Coordinate Reference System. (e.g., EPSG:4269-Geodetic coordinate system for
#' North America, NAD83).
#' @param addxy Logical. If TRUE, adds x and y variables to spatial sf object.
#' @param exportsp Logical. If TRUE, exports spatial object.
#' @param savedata_opts List. See help(savedata_options()) for a list
#' of options. Only used when exportsp = TRUE.  
#' @return \item{spplt}{ sf obect with spatial points and defined CRS. }
#' 
#' If exportsp = TRUE, the sf object is written to specified output.
#' @note If exportsp=TRUE and a shp output format is specified:\cr The ESRI
#' shapefile driver truncates variable names to 10 characters or less.
#' Variable names are changed before export using an internal function
#' (trunc10shp). Name changes are output to the outfolder,
#' 'outshpnm'_newnames.csv.  The returned Spatial object will have original
#' names, before truncating.
#' 
#' If Spatial object has more than 1 record, it cannot be exported.
#' 
#' @author Tracey S. Frescino
#' @keywords data
#' @examples 
#' # Generate an `sf` points object with `spMakeSpatialPoints` for Wyoming plot
#' # data, stored in `FIESTA`
#' spMakeSpatialPoints(xyplt = WYplt, 
#'                     xy.uniqueid = "CN", 
#'                     xvar = "LON_PUBLIC", 
#'                     yvar = "LAT_PUBLIC", 
#'                     xy.crs = 4269)
#' @export spMakeSpatialPoints
spMakeSpatialPoints <- function(xyplt, 
                                xyplt_dsn = NULL, 
                                xy.uniqueid = NULL, 
                                xvar = NULL, 
                                yvar = NULL, 
                                xy.crs = 4269, 
                                addxy = FALSE, 
                                exportsp = FALSE, 
                                savedata_opts = NULL){
  ##############################################################################
  ## DESCRIPTION:
  ## Generates sf points object with defined projection.
  ## Arguments
  ## crs - can be EPSG or PROJ.4 string
  ## ... - Arguments passed to spExportSpatial
  ##############################################################################
 
  ## IF NO ARGUMENTS SPECIFIED, ASSUME GUI=TRUE
  gui <- ifelse(nargs() == 0, TRUE, FALSE)


  ## Check input parameters
  input.params <- names(as.list(match.call()))[-1]

  formallst <- names(formals(spMakeSpatialPoints))
  if (!all(input.params %in% formallst)) {
    miss <- input.params[!input.params %in% formallst]
    stop("invalid parameter: ", toString(miss))
  }

  ## Check parameter lists
  pcheck.params(input.params, savedata_opts=savedata_opts)
  
  ## Set savedata defaults
  savedata_defaults_list <- formals(savedata_options)[-length(formals(savedata_options))]
  
  for (i in 1:length(savedata_defaults_list)) {
    assign(names(savedata_defaults_list)[[i]], savedata_defaults_list[[i]])
  }
  
  ## Set user-supplied savedata values
  if (length(savedata_opts) > 0) {
    for (i in 1:length(savedata_opts)) {
      if (names(savedata_opts)[[i]] %in% names(savedata_defaults_list)) {
        assign(names(savedata_opts)[[i]], savedata_opts[[i]])
      } else {
        stop(paste("Invalid parameter: ", names(savedata_opts)[[i]]))
      }
    }
  }

  ##################################################################
  ## CHECK INPUT PARAMETERS
  ##################################################################
  ## check xyplt
  xypltx <- pcheck.table(xyplt, xyplt_dsn, tabnm="xyplt", 
		caption="XY data table", stopifnull=TRUE, returnDT=FALSE)

  ## check xy.uniqueid
  xypltnmlst <- names(xypltx)
  xy.uniqueid <- pcheck.varchar(var2check=xy.uniqueid, varnm="xy.uniqueid", 
		checklst=xypltnmlst, caption="UniqueID variable - xyplt", 
		warn="xy.uniqueid not in xyplt", gui=gui, stopifnull=TRUE)

  ## check for NA or duplicate values
  if (sum(is.na(xypltx[[xy.uniqueid]])) > 0) stop("NA values in ", xy.uniqueid)
  if (length(unique(xypltx[[xy.uniqueid]])) < nrow(xypltx)) 
    warning("plt records are not unique")

  ## check xvar and yvar
  x <- pcheck.varchar(xvar, varnm="xvar", checklst=xypltnmlst, 
		caption="X variable", gui=gui, stopifnull=TRUE)
  y <- pcheck.varchar(yvar, varnm="yvar", checklst=xypltnmlst, 
		caption="Y variable", gui=gui, stopifnull=TRUE)

  ## check if x = y
  if (x == y) stop("x and y are the same value")

  ## set x and y variables to numeric
  if (!is.numeric(xypltx[[x]])) xypltx[[x]] <- as.numeric(xypltx[[x]])
  if (!is.numeric(xypltx[[y]])) xypltx[[y]] <- as.numeric(xypltx[[y]])

  ## check for NA values
  if (any(is.na(xypltx[[x]])) || any(is.na(xypltx[[y]]))) {
    xypltx <- setDT(xypltx)
    missCN <- xypltx[list(NA,NA), on=c(x,y), xy.uniqueid, nomatch=0, with = FALSE][[1]]
    xypltx <- as.data.frame(xypltx)
    rowp <- ifelse(length(missCN) == 1, "row", "rows")
    message(paste("removing", length(missCN), rowp, "with", xy.uniqueid, "values = NA:", 
		paste(missCN, collapse=", ")))
    setDT(xypltx)
    xypltx <- xypltx[!list(NA,NA), on=c(x,y)]
    xypltx <- as.data.frame(xypltx)
  }

  ## check xy.crs   
  if (is.null(xy.crs)) {
    stop("must include xy.crx") 
  }  
  ### check exportsp
  exportsp <- pcheck.logical(exportsp, varnm="exportsp", 
		title="Export spatial layer?", first="NO", gui=gui)

  ### check addxy
  addxy <- pcheck.logical(addxy, varnm="addxy", 
		title="Add xy variables?", first="NO", gui=gui)

  ##################################################################
  ## DO WORK
  ##################################################################

  ## Make xy.uniqueid a character
  xypltx[[xy.uniqueid]] <- as.character(xypltx[[xy.uniqueid]])
 
  ## Generate sf layer  
  spplt <- sf::st_as_sf(xypltx, coords=c(x,y), crs=xy.crs, 
		stringsAsFactors=FALSE, agr="identity")
  if (is.na(sf::st_crs(spplt))) stop("invalid crs: ", xy.crs) 
  

  ## Add coordinates
  if (addxy) {
    xy.coords <- data.frame(sf::st_coordinates(spplt))
    names(xy.coords) <- c(x,y)
    spplt <- sf::st_sf(data.frame(spplt, xy.coords)) 
  }

  if (exportsp) {
    spExportSpatial(spplt, savedata_opts=savedata_opts) 
  }   
  return(spplt)
}


