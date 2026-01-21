# Community-Size Biases in Statistical Inference of Communities in Temporal Networks

This code implements our comparison of four different methods for community detection in temporal networks.

# Acknowledgements

We thank [Christopher R. Anderson](http://www.math.ucla.edu/~anderson) for writing much of the MathVector and MathMatrix classes, as part of his Math 280 course.

# Usage

The code saves and outputs the community assignments for each of the runs of each community-detection method. 
The saved group assignments for a given layer are outputted to the files ```layer[layer]finalGroups[community 1 size][method]sigma[SBM parameter].txt```. Here [layer] is the index of the given layer, [community 1 size] is the size of community 1 without any offsets $\tau$, [method] is the method being considered (e.g. `bazzi` for the method of [Bazzi et al.](), and [SBM parameter] is the value used for the diagonal elements of the parameter matrices for the SBMs used to generate the network structure for each layer. 

Each row of this text file corresponds to a node, and each column corresponds to a different run (i.e., the entry in the 5th row and 3rd column of the output will be the community assignment of node 5 in run 3). 

