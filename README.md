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

**3. Trajectory Inference with Slingshot**

From now on, we will start using that clustering and data reduction techniques for trajectory inference. The whole process can be done using a single function named slingshot, which is simply a wrapper for the 2 main steps for trajectory inference. The first step of the process is to define the lineages and then fit a curve through the data that defines a trajectory

    ~ Preparing Data for Slingshot:
      Extract UMAP embeddings and clustering results from Seurat.
   
    ~ Running Slingshot:
      Define cell lineages and plot them.
   
    ~ Plotting Lineages:
      Visualize the inferred trajectories.

![Plot_lineages](https://github.com/Divya090597/Trajectory-inference_Slingshot/assets/156469276/eddb09b3-88a1-46c0-b8cc-48147d94c3e8)

**Defining Principal Curves**

Once the clusters are connected, Slingshot allows you to transform them to a smooth trajectory using principal curves. This is an algorithm that iteratively changes an initial curve to better match the data points. It was developed for linear data. To apply it to single-cell data, slingshot adds two enhancements:

It will run principal curves for each ‘lineage’, which is a set of clusters that go from a defined start cluster to some end cluster

Lineages with a same set of clusters will be constrained so that their principal curves remain bundled around the overlapping clusters


**4. Differential Gene Expression with tradeSeq**

    ~ Filtering and Preparing Counts:
      Filter genes to speed up computations.
  
    ~ Fitting Generalized Additive Models (GAMs):
      Fit GAMs to the filtered counts using the inferred trajectories.
  
    ~ Plotting Differential Expression:
      Define a function to plot differentially expressed genes along pseudotime.
  
    ~ Finding Differentially Expressed Genes:
      Perform various tests to identify genes that change with pseudotime, between start and end points, and between lineages.
  
    ~ Plotting Results:
      Plot the most significant differentially expressed genes.

**Conclusion**

This detailed workflow demonstrates the power of integrating multiple tools like Seurat, Slingshot, and tradeSeq for comprehensive scRNA-seq data analysis.








