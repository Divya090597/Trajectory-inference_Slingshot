# Trajectory-inference_Slingshot

Trajectory inference (TI) in single-cell data analysis involves the reconstruction of developmental trajectories and the identification of gene expression changes over time within a cell population. This concept allows researchers to study dynamic changes in gene expression and uncover the underlying biological processes within a single-cell dataset.

When studying dynamic cellular processes like cell differentiation or cellular response to certain stimulus, cells transition from one functional state to another. So, when these cells move between states there are certain set of genes that are activated during earlier states and certain set of genes activating during later states.

**Trajectory inference is primarily used to:**

Reconstruct Developmental Trajectories: This involves ordering cells along a trajectory based on their gene expression profiles, which can 
represent biological processes such as cell differentiation.

Identify Gene Expression Changes: TI helps in identifying genes that are associated with specific lineages in the trajectory or exhibit 
differential expression between lineages, providing insights into biological processes.

**Purpose of Pseudo time**

Ordering Cells: Pseudotime assigns a numerical value to each cell, ordering them along a trajectory based on their gene expression profiles. This ordering is used to represent the progression of cells through a biological process, such as differentiation or a developmental pathway.

Capturing Developmental Dynamics: Pseudotime provides a continuous representation of the underlying biological process, allowing researchers to study and visualize the dynamic changes in gene expression and cell state over the course of development.

**Slingshot consists of two main stages:** 
1) the inference of the global lineage structure and
2) the inference of pseudotime variables for cells along each lineage

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

![Plot_Curves](https://github.com/Divya090597/Trajectory-inference_Slingshot/assets/156469276/c2f7c282-b8a4-466f-93c3-e986aeefc934)

**4. Differential Gene Expression with tradeSeq**

The main way to interpret a trajectory is to find genes that change along the trajectory. There are many ways to define differential expression along a trajectory:

    ~ Expression changes along a particular path (i.e. change with pseudotime)
    ~ Expression differences between branches
    ~ Expression changes at branch points
    ~ Expression changes somewhere along the trajectory

tradeSeq is a algorithm to find trajectory differentially expressed genes. It works by smoothing the gene expression along the trajectory by fitting a smoother using generalized additive models (GAMs), and testing whether certain coefficients are statstically different between points in the trajectory.

    ~ Filtering and Preparing Counts:
      Filter genes to speed up computations.
  
    ~ Fitting Generalized Additive Models (GAMs):
      Fit GAMs to the filtered counts using the inferred trajectories.
  
    ~ Plotting Differential Expression:
      Define a function to plot differentially expressed genes along pseudotime.

![Plot_Clusters](https://github.com/Divya090597/Trajectory-inference_Slingshot/assets/156469276/555d53cf-8721-4bbc-8224-28c22d6cc4ff)

  
    ~ Finding Differentially Expressed Genes:
      Perform various tests to identify genes that change with pseudotime, between start and end points, and between lineages.
  
    ~ Plotting Results:
      Plot the most significant differentially expressed genes.

**Genes that change with pseudo time**

![Genes@change$pseudotime](https://github.com/Divya090597/Trajectory-inference_Slingshot/assets/156469276/8b85eb7f-ecbb-4adf-83be-4c70f34bcc24)

**Genes that change between two pseudotime points**

![bet_2 pseudotime points](https://github.com/Divya090597/Trajectory-inference_Slingshot/assets/156469276/46468a17-1635-4682-ab60-ac0e5e0f5289)

**Genes that are different between lineages**

More interesting are genes that are different between two branches. There are several ways to define “different between branches”, and each have their own functions:

     ~ Different at the end points, using "diffEndTest"
     ~ Different at the branching point, using "earlyDETest"
     ~ Different somewhere in pseudotime the branching point, using "patternTest"
     ~ Note that the last function requires that the pseudotimes between two lineages are aligned.

Using "diffEndTest"

![diff bet 2 lineages](https://github.com/Divya090597/Trajectory-inference_Slingshot/assets/156469276/466263c5-8ca8-42ff-9a15-1c2fc68505ba)

Using "earlyDETest"

![Branch_point_association](https://github.com/Divya090597/Trajectory-inference_Slingshot/assets/156469276/425b081b-704e-4b12-9691-ca5188e4b166)


**Conclusion**

This detailed workflow demonstrates the power of integrating multiple tools like Seurat, Slingshot, and tradeSeq for comprehensive scRNA-seq data analysis.








