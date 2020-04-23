#' Read FASTA file into CGR array
#'
#' \code{fasta_to_array} Returns an array of CGRs with dimension (\code{sqrt(k**4), sqrt(k**4), n })
#'
#' This is a function that wraps the python cgr-view function \code{from_fasta()} and
#' reloads and returns the saved CGR .npy file as an R array.
#'
#' Requires python cgr-view package and Jellyfish on path.
#'
#' @param fasta_file FASTA file to load
#' @param outfile outfile to save
#' @param as_single If TRUE treats all entries as single sequence and return one CGR. If FALSE, treats all entries individually and returns many CGR
#' @param k length of kmer to use
#' @return array of CGRs
#'
fasta_to_array <- function(fasta_file, outfile = "my_cgrs", as_single = FALSE, k = 7) {

  if ( grepl(".*\\.npy$", outfile) ){
    outfile <- gsub("\\.npy", "", outfile)
  }
  k = as.character(k)
  cgr <- reticulate::import("cgr")
  cgr$from_fasta(fasta_file, outfile, as_single,k)
  cgr$load_npy(paste0(outfile, ".npy"))
}

#' Read numpy npy file to R array
#'
#' \code{npy_to_array} Returns a native R array of a Python numpy .npy file in sample, height, width, channel order/ Works for one channnel arrays only
#'
#' @param npy_file: numpy .npy file path
#' @return array
#'
npy_to_array <- function(npy_file){
  cgr <- reticulate::import("cgr")
  a <- cgr$load_npy(npy_file)
  a <- aperm(a, c(3,1,2))
  dim(a) <- c(dim(a), 1)
  return(a)
}

#' Save a single 2D CGR matrix to a PNG image file
#'
#' \code{save_png} saves a single 2D CGR matrix to a PNG file
#'
#' @param cgr_matrix a CGR matrix
#' @param outfile the filename to write
#' @param w width of image in pixels
#' @param h height of image in pixels
#' @param option ViridisLite pallete to use (one of A, B, C, D, E (default))
#' @return NULL
#'
save_png <- function(cgr_matrix, outfile = "my_cgr.png", option = "E", w = 600, h = 600){

  png(outfile, width = w, height = h)
  plot(raster::raster(cgr_matrix),
       col=viridisLite::viridis(n = max(cgr_matrix), begin = 0, end = 1, direction = 1,
                                option = option) )

  dev.off()
}

#' Draw a single 2D CGR matrix to the screen
#'
#' \code{plot_cgr} draws a single 2D CGR matrix to screen
#'
#' @param cgr_matrix a CGR matrix
#' @param option ViridisLite pallete to use (one of A, B, C, D, E (default))
#' @return NULL
#'
plot_cgr <- function(cgr_matrix, option = "E"){
  plot(raster::raster(cgr_matrix),
       col=viridisLite::viridis(n = max(cgr_matrix), begin = 0, end = 1, direction = 1,
                                option = option) )

}
