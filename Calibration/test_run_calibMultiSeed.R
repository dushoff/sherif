n <- as.integer(read.table("nseed.csv"))
for(i in 1:n){
  cmd <- paste0("Rscript test_calibration_multiseed.R 1 ",i," ")
  system(command = cmd, intern = F)
}