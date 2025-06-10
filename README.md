# RPSSFC-NAFW
A Robust Possibilistic Semi-Supervised Fuzzy Clustering Algorithm with Neighborhood-Aware Feature Weighting


The Semi-Supervised Fuzzy C-Means (SSFCM) method integrates class distribution information with fuzzy logic to overcome the challenges of semi-supervised clustering methods. While the inclusion of label information in the objective function improves the quality of the clustering method, semi-supervised fuzzy techniques still encounter important limitations, including 1) sensitivity to noise and outliers, 2) uniform feature importance, 3) neglecting the influences of neighborhood in the clustering process. In this paper, an improved semi-supervised clustering algorithm is presented to address these challenges. First, the algorithm reduces the sensitivity to noise and outliers by integrating the possibilistic fuzzy C-means algorithm into the SSFCM method. Second, a dynamic feature weighting method assigns different weights to the features in each cluster, which improves the performance of the algorithm in imbalanced datasets. Third, the proposed algorithm introduces a neighborhood mechanism that incorporates the neighbor's trade-off weighting and feature weighting strategy considering a strong metric. Finally, a robust kernel metric is used to further improve the performance on complex and nonlinear datasets. Extensive experiments are conducted on several benchmark datasets to evaluate the performance of the proposed method. The results show that the proposed method outperforms the current state-of-the-art techniques. 

# Overview of the RPSSFC-NAFW:
![image](https://github.com/user-attachments/assets/01e162f9-93a0-4481-a78e-14864951b092)



# Comment:
The repository file includes the MATLAB implementation of the RPSSFC-NAFW algorithm.

Comments are written for all steps of the algorithm for better understanding the code. Also, a demo is implemented for ease of running, which runs by importing the data and other necessary algorithm parameters.

To evaluate the proposed method, the UCI benchmark datasets have been used. 





