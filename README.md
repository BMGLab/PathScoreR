# PathScoreR

This is Rcpp extension that has written by [Musa-Sina-Ertugrul](https://github.com/Musa-Sina-Ertugrul) for BMG Lab. To install extension, install .cpp file from this github page. Then choose a location for file such as:
```
path_to_file/bmg_lab.cpp
```
After these steps open R then instal and import Rcpp (If you already installed Rcpp, you can pass installation step):
```r
install.packages("Rcpp")
```
Then import it:
```r
library(Rcpp)
```
Finally you can import bmg_lab.cpp:
```r
Rcpp::sourceCpp("path_to_file/bmg_lab.cpp")
```
This extension has only one main function named as run() other functions are test functions, to run run():
```r
# B_data is dataframe that gene names are rownames
# list = B_Data_list is list format of B_Data
# row_names = rownames(B_Data) 
# groups = groups_2 is group information as
          # NumericVector such as c( 1,1,1,1,2,2,2,2) that data has 2 group
# edges_1 = as.character(Edges[["X.node1"]]) one of columns in edges dataframe this must be
          # CharecterVector Chracters must be names of genes and
          # they must be same with row_names parameter genes
# edges_2 = as.character(Edges[["node2"]]) same of edges_1
# alpha threshold of p_value
# D threshold for length between two nodes that are in network

final<-run( list = B_Data_list, groups = groups_2, row_names = rownames(B_Data),
 edges_1 = as.character(Edges[["X.node1"]]), edges_2 = as.character(Edges[["node2"]]),
 alpha = 0.05, D = 1000)

```
There is one more important function that named as S_i(). That function calculates S score nodes as topological:
```r
# B_data is dataframe that gene names are rownames
# rowNames = rownames(B_Data)
# edges_1 = as.character(Edges[["X.node1"]]) one of columns in edges dataframe this must be
          # CharecterVector Chracters must be names of genes and
          # they must be same with row_names parameter genes
# edges_2 = as.character(Edges[["node2"]]) same of edges_1
# alpha threshold of p_value
# D threshold for length between two nodes that are in network

S_i( rowNames = as.character(rownames(B_Data)),edges_1 = as.vector(as.character(Edges[["X.node1"]])),
edges_2 = as.vector(as.character(Edges[["node2"]])), alpha = 0.05, D = 100)

```
This work based on [NetworkHub article](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3444-7). Thanks for visiting.
