# A generalized single step method that includes any number of hierarchical genomic matrices

Miguel Pérez-Enciso

The Single Step algorithm allows combining information from genotyped and un-genotyped individuals, provided they are connected by a pedigree. Current single step theory is limited to a single genotyping array. 

We present a generalized single step (GSS) method that can accommodate any number of hierarchical molecular datasets (e.g., sequence, high and low density arrays) and pedigree, avoiding imputation. We proof that a similar efficient inversion algorithm exists. The method is recursive, starting with the highest marker density scenario. We illustrate the method with simulation and show that GSS can increase predictive accuracy compared to standard single step. R code is provided so that custom scenarios can be easily compared, either with simulated or real data.

The method developed generalizes extant single step theory to any number of hierarchical molecular relationship matrices, broadening the scenarios where single step can be applied. A topic of particular interest can be ecology field data or human populations where pedigree is not available, but where samples sequenced and genotyped at different densities can exist. GSS can also be a useful tool to optimize allocation of genotyping and / or sequencing resources.
