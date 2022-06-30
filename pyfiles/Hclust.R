library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(arrow)      # for loading dist matrix

df <- arrow::read_feather('/home/liam/anthrax/pangraph/processed_data/anthrax.arrow')
d <- dist(df, method = "euclidean")

hc_a <- agnes(df, method="ward")
sub_grp <- cutree(as.hclust(hc_a),k=4)
mutate(df, cluster = sub_grp)
arrow::write_feather(df,'anthraxcluster.arrow')
