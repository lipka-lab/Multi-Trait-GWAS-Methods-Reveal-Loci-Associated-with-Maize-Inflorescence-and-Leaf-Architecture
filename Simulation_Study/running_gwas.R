setwd("C:/Users/samuelf/Box Sync/brians_paper")
source("run_gemma.R")
library(parallel)
files <- dir("./SimulatedPhenotypes/Setting_1_Null_TRAIT/gemma/")
files <- files[grep("Simulated_Data_", files)]
files <- split(files, 1:length(files))
dir.create("./GWAS")
setwd("./GWAS")
for (i in 1:10){
  file.copy(paste0("../marker_data/k",i,".eigenD.txt"), "./")
  file.copy(paste0("../marker_data/k",i,".eigenU.txt"), "./")
}
file.copy("../marker_data/PCA.txt", "./")
mclapply(files, function(x) {
  for (i in 1:10){
    genoname <- paste0(gsub(".fam", "", x),"_chr", i)
    file.copy(paste0("../SimulatedPhenotypes/Setting_1_Null_TRAIT/gemma/",x), "./")
    file.copy(paste0("../marker_data/chr",i,".bed"), "./")
    file.copy(paste0("../marker_data/chr",i,".bim"), "./")
    file.rename(paste0("chr",i,".bed"), paste0(genoname, ".bed"))
    file.rename(paste0("chr",i,".bim"), paste0(genoname, ".bim"))
    file.rename(x, paste0(genoname,".fam"))
  run_gemma(gemma_path = "~/gemma",
            geno = genoname, 
            fixed = "PCA.txt",
            d = paste0("k",i,".eigenD.txt"), 
            u = paste0("k",i,".eigenU.txt"),
            trait_col = c(1, 2, 3), 
            lmm = 2,
            miss = 0.99,
            maf = 0.001,
            r2 = 0.99999,
            out_name = genoname)
  }
  print("Done! Starting next file...")
},
mc.preschedule = FALSE,
mc.cores = detectCores())

