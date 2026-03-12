# A generalized single step method that includes any number of hierarchical genomic matrices

**Miguel Pérez-Enciso**

mperezenciso@gmail.com

The Single Step algorithm allows combining information from genotyped and un-genotyped individuals, provided they are connected by a pedigree. Current single step theory is limited to a single genotyping array. 

We present a generalized single step (GSS) method that can accommodate any number of hierarchical molecular datasets (e.g., sequence, high and low density arrays) and pedigree, avoiding imputation. We proof that a similar efficient inversion algorithm exists. The method is recursive, starting with the highest marker density scenario. We illustrate the method with simulation and show that GSS can increase predictive accuracy compared to standard single step. R code is provided so that custom scenarios can be easily compared, either with simulated or real data.

The method developed generalizes extant single step theory to any number of hierarchical molecular relationship matrices, broadening the scenarios where single step can be applied. A topic of particular interest can be ecology field data or human populations where pedigree is not available, but where samples sequenced and genotyped at different densities can exist. GSS can also be a useful tool to optimize allocation of genotyping and / or sequencing resources.

## Theory
In the following, we assume the usual linear mixed model:

y = Xb + Wu + e,

where y contains the phenotypes, b is the fixed effects vector, random genetic values are in u, residuals in e, and X and W are incidence matrices. SS theory considers a subset in u, u1 comprising individuals without molecular data and subset u2 which have been genotyped. Both u1 and u2 are related through a common pedigree (P) that results in numerator relationship matrix A. When marker information on m SNPs is available, a genomic relationship matrix (GRM) G can be computed as

$$\mathbf{G} = \frac{(\mathbf{M} - 2\boldsymbol{p}) (\mathbf{M} - 2\boldsymbol{p})'}{2 \boldsymbol{p}'(1 - \boldsymbol{p})}$$

where M is a n individual by m marker matrix containing genotypes (coded as 0,1,2) and p, a m-dimension vector with corresponding allele frequencies. Legarra et al. (2009) and Christensen and Lund (2010) showed that the variance of breeding values can be approximated by

$$\mathbf{H} = Var \left( \begin{matrix} \boldsymbol{u}_1 \\ \boldsymbol{u}_2 \end{matrix} \middle| \mathbf{P}, \mathbf{M}, \dots \right) = 
\begin{pmatrix} 
\mathbf{A}_{11} - \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21} + \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}_{22}\mathbf{A}_{22}^{-1}\mathbf{A}_{21} & \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G}_{22} \\ 
\mathbf{A}_{21}\mathbf{A}_{22}^{-1}\mathbf{G}_{22} & \mathbf{G}_{22} 
\end{pmatrix} \sigma_u^2$$

We generalize SS current theory to an arbitrary number of hierarchized molecular information datasets. Initially, without loss of generalization, let us consider three subsets: no molecular information (only P available), SNP array (M) and sequence (MS) information. We consider that sequence is available only for a subset u3 of those individuals that are also genotyped with the array.  We further assume that all markers in M are also present in MS, and that individuals are connected through a pedigree P. Therefore, the goal is to obtain

$$Var \left( \begin{matrix} \boldsymbol{u}_1 \\ \boldsymbol{u}_2 \\ \boldsymbol{u}_3 \end{matrix} \middle| \mathbf{P}, \mathbf{M}, \mathbf{M}_S, \dots \right) = 
\begin{pmatrix} 
\mathbf{H}_{11} & \mathbf{H}_{12} & \mathbf{H}_{13} \\ 
\mathbf{H}_{21} & \mathbf{H}_{22} & \mathbf{H}_{23} \\ 
\mathbf{H}_{31} & \mathbf{H}_{32} & \mathbf{H}_{33} 
\end{pmatrix} \sigma_u^2 = \mathbf{H} \sigma_u^2,$$

First, we note that

$$\mathbf{H}_{33} = Var(\boldsymbol{u}_3 | \mathbf{P}, \mathbf{M}, \mathbf{M}_S) \sim Var(\boldsymbol{u}_3 | \mathbf{M}_S) = \mathbf{S} \sigma_u^2$$

where S is the genomic relationship matrix computed applying eq. (1) to sequence data. Equation above means that sequence information makes it ignorable the contribution from the array (and that of the pedigree), i.e., the same assumption as in standard SS. Next, following the same reasoning as in single step theory (eq. 2), except that G is used instead of A, and S instead of G:

$$\mathbf{H}_{22} = Var(\boldsymbol{u}_2 | \mathbf{P}, \mathbf{M}, \mathbf{M}_S) \sim Var(\boldsymbol{u}_2 | \mathbf{M}, \mathbf{M}_S) = (\mathbf{G}_{22} - \mathbf{G}_{23} \mathbf{G}_{33}^{-1} \mathbf{G}_{32} + \mathbf{G}_{23} \mathbf{G}_{33}^{-1} \mathbf{S} \mathbf{G}_{33}^{-1} \mathbf{G}_{32}) \sigma_u^2,$$

and by the same token

$$\mathbf{H}_{23} = Cov(\boldsymbol{u}_2, \boldsymbol{u}_3 | \mathbf{P}, \mathbf{M}, \mathbf{M}_S) \sim Cov(\boldsymbol{u}_2, \boldsymbol{u}_3 | \mathbf{M}, \mathbf{M}_S) = \mathbf{G}_{23} \mathbf{G}_{33}^{-1} \mathbf{S} \sigma_u^2.$$

Next, note

$$
\begin{aligned}
f(\mathbf{u}_1 | \mathbf{u}_2, \mathbf{u}_3) &= N \left[ \begin{pmatrix} \mathbf{A}_{12} & \mathbf{A}_{13} \end{pmatrix} 
\begin{pmatrix} \mathbf{A}_{22} & \mathbf{A}_{23} \\ \mathbf{A}_{32} & \mathbf{A}_{33} \end{pmatrix}^{-1} 
\begin{pmatrix} \mathbf{u}_2 \\ \mathbf{u}_3 \end{pmatrix}, \right. \\
& \quad \left. \{ \mathbf{A}_{11} - \begin{pmatrix} \mathbf{A}_{12} & \mathbf{A}_{13} \end{pmatrix} 
\begin{pmatrix} \mathbf{A}_{22} & \mathbf{A}_{23} \\ \mathbf{A}_{32} & \mathbf{A}_{33} \end{pmatrix}^{-1} 
\begin{pmatrix} \mathbf{A}_{21} \\ \mathbf{A}_{31} \end{pmatrix} \} \sigma_u^2 \right]
\end{aligned}
$$

Now let us take the submatrix of H that contains the joint variance of u2 and u3, i.e., those individuals with molecular data, either sequence or SNP array:

$$Var \left( \begin{matrix} \mathbf{u}_2 \\ \mathbf{u}_3 \end{matrix} \middle| \mathbf{P}, \mathbf{M}, \mathbf{M}_S, \dots \right) = \begin{pmatrix} \mathbf{H}_{22} & \mathbf{H}_{23} \\ \mathbf{H}_{32} & \mathbf{H}_{33} \end{pmatrix} \sigma_u^2$$

From above, integrating out (u2, u3) in (3) using (4), it follows that 

$$
\begin{aligned}
\mathbf{H}_{11} &= Var(u_1 | \mathbf{P}, \mathbf{M}, \mathbf{M}_S) \\\\
&= \left[ \mathbf{A}_{11} - \begin{pmatrix} \mathbf{A}_{12} & \mathbf{A}_{13} \end{pmatrix} 
\begin{pmatrix} \mathbf{A}_{22} & \mathbf{A}_{23} \\\\ \mathbf{A}_{32} & \mathbf{A}_{33} \end{pmatrix}^{-1} 
\begin{pmatrix} \mathbf{A}_{21} \\\\ \mathbf{A}_{31} \end{pmatrix} \right] \sigma_u^2 \\\\
&+ \left[ \begin{pmatrix} \mathbf{A}_{12} & \mathbf{A}_{13} \end{pmatrix} 
\begin{pmatrix} \mathbf{A}_{22} & \mathbf{A}_{23} \\\\ \mathbf{A}_{32} & \mathbf{A}_{33} \end{pmatrix}^{-1} 
\begin{pmatrix} \mathbf{H}_{22} & \mathbf{H}_{23} \\\\ \mathbf{H}_{32} & \mathbf{H}_{33} \end{pmatrix} 
\begin{pmatrix} \mathbf{A}_{22} & \mathbf{A}_{23} \\\\ \mathbf{A}_{32} & \mathbf{A}_{33} \end{pmatrix}^{-1} 
\begin{pmatrix} \mathbf{A}_{21} \\\\ \mathbf{A}_{31} \end{pmatrix} \right] \sigma_u^2
\end{aligned}
$$

and 

$$
\begin{aligned}
\mathbf{H}_{1, j>1} &= \begin{pmatrix} \mathbf{H}_{12} & \mathbf{H}_{13} \end{pmatrix} \\\\
&= Cov(\mathbf{u}_1, [\mathbf{u}_2, \mathbf{u}_3]' | \mathbf{P}, \mathbf{M}, \mathbf{M}_S) \\\\
&= \begin{pmatrix} \mathbf{A}_{12} & \mathbf{A}_{13} \end{pmatrix} 
\begin{pmatrix} \mathbf{A}_{22} & \mathbf{A}_{23} \\\\ \mathbf{A}_{32} & \mathbf{A}_{33} \end{pmatrix}^{-1} 
\begin{pmatrix} \mathbf{H}_{22} & \mathbf{H}_{23} \\\\ \mathbf{H}_{32} & \mathbf{H}_{33} \end{pmatrix} \sigma_u^2
\end{aligned}
$$

Because the generalized H matrix has the same structure as in the standard SS theory, its inverse is also given by

$$
\begin{aligned}
\mathbf{H}^{-1} &= \mathbf{A}^{-1} + 
\begin{pmatrix} 
0 & \mathbf{0} \\\\ 
\mathbf{0} & \begin{bmatrix} \mathbf{H}_{22} & \mathbf{H}_{23} \\\\ \mathbf{H}_{32} & \mathbf{H}_{33} \end{bmatrix}^{-1} - \begin{bmatrix} \mathbf{A}_{22} & \mathbf{A}_{23} \\\\ \mathbf{A}_{32} & \mathbf{A}_{33} \end{bmatrix}^{-1} 
\end{pmatrix}
\end{aligned}
$$

which requires the inverse of matrices with maximum dimension the number of individuals with molecular data. 

Importantly, the theory above lends itself to generalize to any number of 
marker datasets, as long as they are hierarchical. Suppose we have a list 
of $b$ marker datasets (blocks) with covariance matrices 
$\mathbf{G}_{[1]}, \mathbf{G}_{[2]}, \dots, \mathbf{G}_{[b]}$, such that 
block 1 corresponds to ungenotyped individuals and block $b$ to those 
with highest genotyped density; in the case above $\mathbf{G}_{[1]} \equiv \mathbf{A}$, 
and $\mathbf{G}_{[3]} \equiv \mathbf{S}$. Then for any $i$-th block:


$$
\begin{aligned}
\mathbf{H}_{ii} &= \left[ \mathbf{G}_{[i]ii} - \mathbf{G}_{[i]i,j>i} \mathbf{G}_{[i]j>i,j>i}^{-1} \mathbf{G}_{[i]j>i,i} \right. \\\\
&\quad + \left. \mathbf{G}_{[i]i,j>i} \mathbf{G}_{[i]j>i,j>i}^{-1} \mathbf{H}_{j>i,j>i} \mathbf{G}_{[i]j>i,j>i}^{-1} \mathbf{G}_{[i]j>i,i} \right] \sigma_u^2
\end{aligned}
$$

and

$$
\begin{aligned}
\mathbf{H}_{i, j>i} &= \left[ \mathbf{G}_{[i]i, j>i} \mathbf{G}_{[i]j>i, j>i}^{-1} \mathbf{H}_{j>i, j>i} \right] \sigma_u^2
\end{aligned}
$$

The set of equations (6) and (7) is computed recursively, starting with 
$\mathbf{H}_{bb} = \mathbf{G}_{[b]bb} \sigma_u^2$, the individuals with 
highest marker density, and can be applied to any number of hierarchical 
marker sets. However, note that each step requires the inverse of submatrix 
$\mathbf{G}_{j>i, j>i}$ (eq. 6), with order of magnitude the sum of individuals 
in each of the blocks $i+1$ to $b$. Therefore, having many blocks of large 
sizes can be computationally demanding. Further, the inverse of $\mathbf{H}$ 
can also be obtained from:

$$
\begin{aligned}
\mathbf{H}^{-1} &= \mathbf{G}_{[1]}^{-1} + 
\begin{pmatrix} 
0 & \mathbf{0} \\\\ 
\mathbf{0} & \left( \mathbf{H}_{i>1, i>1}^{-1} - \mathbf{G}_{[1]i>1, j>1}^{-1} \right)
\end{pmatrix}^{-1}
\end{aligned}
$$

where subscript 'j>i' refers to matrix block with indices i+1 to b. As in classical SS, GSS is optimal computationally when the first block, made-up of ungenotyped individuals, (G1) is the largest and equal to A, since its inverse can be computed in linear time and is sparse.

