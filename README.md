# Trajectory-inference_Slingshot

**1. Data Acquisition and Preparation**
   
    ~ Creating a Data Directory:
      Start by creating a directory named "data" to store the files.
      
    ~ Downloading the Data:
      Download a compressed data file from the NCBI GEO database using the "download.file" function and specify the 
      destination path.
      
    ~ Decompressing the Data:
      After downloading, you decompress the file using gunzip.
      
    ~ Loading and Converting the Data:
      load the data into R and convert it to a sparse matrix, which is then saved as an "RDS file" for efficient storage

**2. Seurat Workflow**
   
    ~ Loading Data:
      The data is loaded from the RDS file and subsetted for further analysis.
    ~ Seurat Analysis Pipeline:
      Perform standard preprocessing and clustering analysis using Seurat.
    ~ Creating Seurat Object
    ~ Normalizing data
    ~ Finding Variable Features
    ~ Scaling data
    ~ Principal Component Analysis (PCA)
    ~ Finding Neighbours and Clustering
    ~ UMAP Dimentionality Reduction
    ~ Plotting Clusters:
      Visualize the clusters using UMAP.

  ![Umap_clusters](https://github.com/Divya090597/Trajectory-inference_Slingshot/assets/156469276/aa594455-2359-484f-ad75-c76681d23c3e)





