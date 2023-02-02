######################################################################################
#### Procedure to download and decompress the bioclimatic data from WorldClim 2.1  ###
#### and from CHELSA-TraCE21K from the contemporary period until 22,000 years ago, ###
#### in steps of 500 years                                                         ### 
######################################################################################

# Download WorldClim 2.1 contemporary data

download.file("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_bio.zip", 
              "/home/shared/pedro.braga/chapter-BiogHistPhyloRelatedSpatScales/data/rasters/wc2.1_2.5m_bio.zip")

# Download CHELSA Trace (using WorldClim 2.1 as baseline climate) data

chelsa_trace_prefix_url <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/"
export_dir <- "~/PaleoClimPoolCommDivPhyloStruct/R1/data/rasters/CHELSA_TraCE21k/bioclim"

bio_vars <- c(
  "bio01", "bio04",
  "bio12", "bio15"
)

time_ID <- seq(from = -200,
               to = 20, 
               by = 5)

### Iterative procedure that downloads one dataset at a time
## See preferred alternative below this chunk

filename <- paste("CHELSA_TraCE21k", 
                  bio_vars, 
                  time_ID, sep = "_")

for(bio_vars.i in bio_vars){
  for(time_ID.k in time_ID){
    download.file(paste(chelsa_trace_prefix_url, filename, "_V1.0.tif", sep = ""),
                  paste(export_dir, "/", filename, "_V1.0.tif", sep = ""))
  }
}

### Iterative procedure that downloads multiple sets at once
# requires 'purr'

# We will separately specify periodds

# for CHELSA

walk(bio_vars, function(x, y) {
  purrr::map(x, 
             ~sprintf("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_%s_%s_V1.0.tif", ., time_ID)) %>% 
    flatten_chr() -> urls
  #print(paste(export_dir, "/", basename(urls), sep = ""))
  download.file(urls, paste(export_dir, "/", basename(urls), sep = ""), method = "libcurl")
}) 

# test not run 
# download.file(
#   "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_bio01_20_V1.0.tif",
#   paste(export_dir, "/", "CHELSA_TraCE21k_bio01_20_V1.0.tif", sep = ""))
 
# download.file(
#   "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_bio12_20_V1.0.tif",
#   paste(export_dir, "/", "CHELSA_TraCE21k_bio12_20_V1.0.tif", sep = ""))
# 

renaming.filenames <- data.frame(time_ID = seq(from = -200,
                                               to = 20, 
                                               by = 5),
                                 time_ID_pad = c(sprintf("%04d", time_ID[1:40]),
                                                 sprintf("%03d", time_ID[41:45]))
)

renaming.filenames <- rbind(renaming.filenames, 
                            renaming.filenames,
                            renaming.filenames)

renaming.filenames[1,1]

file.rename(from = rasterFilenames, 
            to = sub(pattern = "CHELSA_TraCE21k_ ", 
                     replacement= "CHELSA_TraCE21k_", 
                     rasterFilenames))

for(i in 1:nrow(renaming.filenames)){
  file.rename(from = rasterFilenames, 
              to = sub(pattern = paste0(renaming.filenames[i, 1]), 
                       replacement= paste0(renaming.filenames[i, 2]), 
                       rasterFilenames))
}

