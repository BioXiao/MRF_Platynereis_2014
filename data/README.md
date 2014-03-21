##Platynereis dumerilii single cells gene expression dataset
#In this folder are located all the files required to run the EM clustering algorithm 
- `binary_86_genes.tab` contains the binarised dataset for the 32,203 cubes in the brain of Platynereis for the 86 considered genes. The first columns is a unique identifier for each cube. Each following column represents one gene in the same order as described in `genes_names_list.tab`.
- `neighbouring_graph.nei` contains the encoded spatial relationship between cubes. For each line (corresponding to a cube) the column is the number of neighbours. The following columns are the ids of the neighbouring cells.

#Additionally the following files are available
- `genes_names_list.tab` is a list of the 86 genes considered
- `3D_coordinates.csv` is a CSV files that contains the 3D coordinates for each of the 32,203 cubes in the dataset.
