```R
install.packages("adegenet")
library(adegenet)
install.packages("rgl")
library(rgl)
library(gridExtra)
library(htmltools)
library(plotly)

# Load the dissimilarity index CSV file
dist_matrix <- as.matrix(read.csv("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_mperc.csv", row.names = 1))
dist_matrix <- as.matrix(read.csv("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_genic.csv", row.names = 1))

# Verify the class and structure of the matrix
class(dist_matrix)
str(dist_matrix)

# Perform classical multidimensional scaling (MDS)
mds_result <- cmdscale(dist_matrix, k = 3)

# Plot the heatmap
heatmap(dist_matrix)

# Perform PCA on the dissimilarity matrix
pca <- prcomp(mds_result)
pca <- prcomp(dist_matrix)
summary(pca)

# Plot the PCA
plot(pca$x[, 1], pca$x[, 2], pch = 16, xlab = "PC1", ylab = "PC2", main = "PCA Plot")

# Plot the top 3 principal components in 3D
plot3d(pca$x[, 1], pca$x[, 2], pca$x[, 3], col = "blue", xlab = "PC1", ylab = "PC2", zlab = "PC3")
text3d(pca$x[, 1], pca$x[, 2], pca$x[, 3], texts = rownames(pca$x), adj = c(1.5,1.5), cex = 0.8)
rgl::rglwidget()

#Get sample info
sample_info <- read.csv("//jic-hpc-data/Group-Scratch/Saskia-Hogenhout/tom_heaven/Aphididae/snp_calling/Myzus/persicae/biello/PCA_file_host.csv", row.names = 1)
row.names(sample_info) <- sample_info[, 1]
sample_info <- sample_info[, -1]
point_labels <- rownames(pca$x)
sample_info <- sample_info[point_labels, ]
sample_info <- sample_info[, -c(1, 2, 3, 4, 6, 11)]
sample_info <- sample_info[, c(5, 3, 2, 4, 1)]

# Color the points based on 'JIC' or 'Bass' information
plot3d(pca$x[, 1], pca$x[, 2], pca$x[, 3], col = "black", xlab = "PC1", ylab = "PC2", zlab = "PC3")
text3d(pca$x[, 1], pca$x[, 2], pca$x[, 3], texts = rownames(pca$x), adj = c(1.5, 1.5), cex = 0.8)
label_colors <- c("red", "blue") 
for (label in rownames(pca$x)) {
  if (sample_info[label, 1] == "JIC") {
    color <- label_colors[1]  # 'JIC' color
  } else if (sample_info[label, 1] == "Bass") {
    color <- label_colors[2]  # 'Bass' color
  } else {
    color <- "black"  # Default color for other labels
  }
  points3d(pca$x[as.character(label), 1], pca$x[as.character(label), 2], pca$x[as.character(label), 3], col = color, size = 10)
}
rgl::rglwidget()

# Color the points based on Continent information
plot3d(pca$x[, 1], pca$x[, 2], pca$x[, 3], col = "black", xlab = "PC1", ylab = "PC2", zlab = "PC3")
text3d(pca$x[, 1], pca$x[, 2], pca$x[, 3], texts = rownames(pca$x), adj = c(1.5, 1.5), cex = 0.8)
label_colors <- c("red", "blue", "green", "yellow", "pink") 
for (label in rownames(pca$x)) {
  if (sample_info[label, 2] == "Europe") {
    color <- label_colors[1]  
  } else if (sample_info[label, 2] == "Australia") {
    color <- label_colors[2]  
  } else if (sample_info[label, 2] == "America") {
    color <- label_colors[3] 
  } else if (sample_info[label, 2] == "Asia") {
    color <- label_colors[4] 
  } else if (sample_info[label, 2] == "Africa") {
    color <- label_colors[5] 
  } else {
    color <- "black"
  }
  points3d(pca$x[as.character(label), 1], pca$x[as.character(label), 2], pca$x[as.character(label), 3], col = color, size = 10)
}
rgl::rglwidget()

# Color the points based on simple host information
plot3d(pca$x[, 1], pca$x[, 2], pca$x[, 3], col = "black", xlab = "PC1", ylab = "PC2", zlab = "PC3")
text3d(pca$x[, 1], pca$x[, 2], pca$x[, 3], texts = rownames(pca$x), adj = c(1.5, 1.5), cex = 0.8)
label_colors <- c("red", "blue", "green", "yellow", "pink", "orange") 
for (label in rownames(pca$x)) {
  if (sample_info[label, 4] == "Brassica") {
    color <- label_colors[1]  
  } else if (sample_info[label, 4] == "Sugar beet") {
    color <- label_colors[2]  
  } else if (sample_info[label, 4] == "Lab") {
    color <- label_colors[3] 
  } else if (sample_info[label, 4] == "Solanacae") {
    color <- label_colors[4] 
  } else if (sample_info[label, 4] == "Prunus") {
    color <- label_colors[5] 
  } else if (sample_info[label, 4] == "Tobacco") {
    color <- label_colors[6] 
  } else {
    color <- "black"
  }
  points3d(pca$x[as.character(label), 1], pca$x[as.character(label), 2], pca$x[as.character(label), 3], col = color, size = 10)
}
rgl::rglwidget()

```
With python:
```python
from sklearn.manifold import MDS
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import plotly.express as px
import plotly.graph_objects as go


#Collect input
df = pd.read_csv('snp_calling/Myzus/persicae/biello/gatk/p_distance/p_dis_193_mperc.csv', index_col=0)
acc =  df.index.tolist()

#MDS projection dimensionality
n_components = 3                                                                 
embedding = MDS(n_components=3,random_state=20348, dissimilarity='precomputed')  #dissimilarity parameter shows D is already known
nm_scores = embedding.fit_transform(df)

print('MDS projected coordinates matrix (nm_matrix) has dimensions: ', np.shape(nm_scores))

nm_x_proj = nm_scores[:,0]
nm_y_proj = nm_scores[:,1]
nm_z_proj = nm_scores[:,2]

color_discrete_mapD = {'1': 'grey', 
                      '2': 'orange', 
                      '3': 'red',
                      '4': 'green',
                      '5': 'blue',
                      '6': 'purple',
                      '7': 'cyan',
                      'U': 'black',
                      '1-4': 'pink',
}

#Colour for core diversity set:
cores =  ['A105', 'Crespys', 'FR15', 'G006', 'I1', 'K40', 'K66', 'ligustri', 'MISC67', 'NIC_410G', 'NIC_410R', 'NL_IRS', 'S2', 'S6', 'S8', 'S12', 'S15', 'S24', 'S28', 'S34', 'S38', 'S39', 'S42', 'S50', 'S55', 'S58', 'S62', 'S66', 'S67', 'S68', 'S71', 'S73', 'S75', 'S76', 'S86', 'S89', 'S90', 'S91', 'S92', 'S93', 'S99', 'S113']
colour = [ '2' if x in cores else  '1'  for x in acc]

df1 = pd.DataFrame({'acc': acc, 'colour': colour, 'PC1': nm_x_proj, 'PC2': nm_y_proj, 'PC3': nm_z_proj, },
columns=['acc', 'colour', 'PC1', 'PC2', 'PC3'],
index = acc,
)

fig = px.scatter_3d(df1, x='PC1', y='PC2', z='PC3', hover_data = ['acc', 'colour'], color="colour", 
   color_discrete_map = color_discrete_mapD, title = "p_dis_mperc")
fig.update_traces(marker_size = 4)
fig.write_html("cores-python-MDS.html")

#Collect sample info:
sampler = {}
geography = {}
host = {}
head = ['#', 'sample.id', 'EV1', 'EV2', 'EV3', 'EV4', 'host', 'Collection', 'Country', 'Continent', 'simple_host', 'Collected_by', 'dummy']
with open('snp_calling/Myzus/persicae/biello/PCA_file_host.csv') as inp:
  next(inp)
  for line in inp:
      A =  line.strip().split(',')
      d = dict(zip(head, A))
      sampler[d['sample.id']] = d['Collected_by']
      geography[d['sample.id']] = d['Continent']
      host[d['sample.id']] = d['simple_host']

#Colour for sampler
colour_settings = {
 'JIC': 'JIC', 
 'Bass': 'Bass'
}

color_discrete_mapD = {'JIC': 'grey', 
                      'Bass': 'orange', 
                      'U': 'black',
}

colour = [colour_settings[sampler[x]]  for x in acc]
df3 = pd.DataFrame({'acc': acc, 'colour': colour, 'PC1': nm_x_proj, 'PC2': nm_y_proj, 'PC3': nm_z_proj, },
columns=['acc', 'colour', 'PC1', 'PC2', 'PC3'],
index = acc,
)

fig = px.scatter_3d(df3, x='PC1', y='PC2', z='PC3', hover_data = ['acc', 'colour'], color="colour", 
   color_discrete_map = color_discrete_mapD, title = "Sampler")
fig.update_traces(marker_size = 4)
fig.write_html("Sampler-python-MDS.html")

#Colour for continent
colour_settings = {
 'Australia': 'Australia', 
 'Europe': 'Europe', 
 'Africa': 'Africa', 
 'America': 'America', 
 'Asia': 'Asia'
}

color_discrete_mapD = {'Australia': 'grey', 
                      'Europe': 'orange', 
                      'Africa': 'red',
                      'America': 'green',
                      'Asia': 'blue',
                      'U': 'black',
}

colour = [colour_settings[geography[x]]  for x in acc]
df3 = pd.DataFrame({'acc': acc, 'colour': colour, 'PC1': nm_x_proj, 'PC2': nm_y_proj, 'PC3': nm_z_proj, },
columns=['acc', 'colour', 'PC1', 'PC2', 'PC3'],
index = acc,
)

fig = px.scatter_3d(df3, x='PC1', y='PC2', z='PC3', hover_data = ['acc', 'colour'], color="colour", 
   color_discrete_map = color_discrete_mapD, title = "Geography")
fig.update_traces(marker_size = 4)
fig.write_html("Geography-python-MDS.html")

#Colour for host
colour_settings = {
 'Brassica': 'Brassica', 
 'Lab': 'Lab',
 'Solanacae': 'Solanacae',
 'Prunus': 'Prunus',
 'Sugar beet': 'Sugar beet',
 'Tobacco': 'Tobacco',
 'Others': 'Others'
}

color_discrete_mapD = {'Brassica': 'grey', 
                      'Lab': 'orange', 
                      'Solanacae': 'red',
                      'Prunus': 'green',
                      'Sugar beet': 'blue',
                      'Tobacco': 'purple',
                      'Others': 'cyan',
                      'U': 'black',
}

colour = [colour_settings[host[x]]  for x in acc]
df3 = pd.DataFrame({'acc': acc, 'colour': colour, 'PC1': nm_x_proj, 'PC2': nm_y_proj, 'PC3': nm_z_proj, },
columns=['acc', 'colour', 'PC1', 'PC2', 'PC3'],
index = acc,
)

fig = px.scatter_3d(df3, x='PC1', y='PC2', z='PC3', hover_data = ['acc', 'colour'], color="colour", 
   color_discrete_map = color_discrete_mapD, title = "Host")
fig.update_traces(marker_size = 4)
fig.write_html("Host-python-MDS.html")

#Outputs were moved to the R folder with other figures
```