# Granule
A MATLAB script to identify granules in fluorescence image data. This scripts works in 2020b and should work in newer versions as well.

This script is described in: Jeremy C. Hunn, Katherine M. Hutchinson, Joshua B. Kelley, Daniel Reines  AN ALGORITHM TO QUANTIFY INDUCIBLE PROTEIN CONDENSATES IN EUKARYOTIC CELLS. bioRxiv DOI:...

This script takes two input images, one of which will be used to make individual masks for each cell (here the nucleus based on histone signal specifically, but whole cell masks could be used) and the other of which is the protein of interest (Nab) that will be assessed for granule formation. The default values below work for our Nab3 and histone images, but may need to be adjusted for other systems or conditions.
The script will record the granules detected over a range of size cutoffs for granules and save that data as an excel file to help assess the appropriate parameters for scoring a granule.  

We have provided examples images for yeast expressing histone-RFP to define the nucleus, and Nab3-GFP to assess for granule formation.  Nab3 will form nuclear granules in response to glucose starvation.   

The live script (.mlx) is the version we anticipate being used, but if you need the normal .m file, it is available as "granulinator_dead.m". Live scripts can be more intuitive to use, but will take up more resources. 
