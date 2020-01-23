run_gemma <- function(gemma_path = NULL,
                      wd = ".",
                      geno = NULL,
                      fixed = NULL,
                      kin = NULL,
                      d = NULL,
                      u = NULL,
                      lmm = 2,
                      miss = 0.05,
                      maf = 0.01,
                      r2 = 0.9999,
                      trait_col = 1,
                      out_name = NULL){
  if (!is.null(kin)) {
    d <- NULL
    u <- NULL}
  if (is.null(fixed)) {
    if (!is.null(d) &
        !is.null(u)) {
      system(command = paste(gemma_path, 
                             "--bfile",
                             paste0(wd,"/", geno),
                             "-lmm", lmm,
                             "-miss", miss,
                             "-maf",maf,
                             "-r2", r2,
                             "-n", paste(trait_col, collapse = " "),
                             "-d", paste0(wd,"/", d),
                             "-u", paste0(wd,"/", u),
                             "-o", out_name))
    } else {
      system(command = paste(gemma_path, 
                             "--bfile",
                             paste0(wd,"/", geno),
                             "-lmm", lmm,
                             "-miss", miss,
                             "-maf",maf,
                             "-r2", r2,
                             "-n", paste(trait_col, collapse = " "),
                             "-k", paste0(wd,"/", kin),
                             "-o", out_name))
    }
  } else {
    if (!is.null(d) &
        !is.null(u)) {
      system(command = paste(gemma_path, 
                             "--bfile",
                             paste0(wd,"/", geno),
                             "-lmm", lmm,
                             "-miss", miss,
                             "-maf",maf,
                             "-r2", r2,
                             "-c",fixed,
                             "-n", paste(trait_col, collapse = " "),
                             "-d", paste0(wd,"/", d),
                             "-u", paste0(wd,"/", u),
                             "-o", out_name))
    } else {
      system(command = paste(gemma_path, 
                             "--bfile",
                             paste0(wd,"/", geno),
                             "-lmm", lmm,
                             "-miss", miss,
                             "-maf",maf,
                             "-r2", r2,
                             "-c",fixed,
                             "-n", paste(trait_col, collapse = " "),
                             "-k", paste0(wd,"/", kin),
                             "-o", out_name))
    }
  }
}
