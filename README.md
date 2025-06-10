# RPSSFC-NAFW
A Robust Possibilistic Semi-Supervised Fuzzy Clustering Algorithm with Neighborhood-Aware Feature Weighting


The Semi-Supervised Fuzzy C-Means (SSFCM) method integrates class distribution information with fuzzy logic to overcome the challenges of semi-supervised clustering methods. While the inclusion of label information in the objective function improves the quality of the clustering method, semi-supervised fuzzy techniques still encounter important limitations, including 1) sensitivity to noise and outliers, 2) uniform feature importance, 3) neglecting the influences of neighborhood in the clustering process. In this paper, an improved semi-supervised clustering algorithm is presented to address these challenges. First, the algorithm reduces the sensitivity to noise and outliers by integrating the possibilistic fuzzy C-means algorithm into the SSFCM method. Second, a dynamic feature weighting method assigns different weights to the features in each cluster, which improves the performance of the algorithm in imbalanced datasets. Third, the proposed algorithm introduces a neighborhood mechanism that incorporates the neighbor's trade-off weighting and feature weighting strategy considering a strong metric. Finally, a robust kernel metric is used to further improve the performance on complex and nonlinear datasets. Extensive experiments are conducted on several benchmark datasets to evaluate the performance of the proposed method. The results show that the proposed method outperforms the current state-of-the-art techniques. 

# Overview of the RPSSFC-NAFW:
![image](https://github.com/user-attachments/assets/01e162f9-93a0-4481-a78e-14864951b092)



# Comment:
The repository file includes the MATLAB implementation of the RPSSFC-NAFW algorithm.

Comments are written for all steps of the algorithm for better understanding the code. Also, a demo is implemented for ease of running, which runs by importing the data and other necessary algorithm parameters.

To evaluate the proposed method, the UCI benchmark datasets have been used. 

## Condition and terms to use any sources of this project (Codes, Datasets, etc.):

1) Please cite the following papers:

[1] M. Hashemzadeh, A. Golzari Oskouei, and N. Farajzadeh, "New fuzzy C-means clustering method based on feature-weight and cluster-weight learning," Applied Soft Computing, vol. 78, pp. 324-345, 2019/05/01/ 2019, doi: https://doi.org/10.1016/j.asoc.2019.02.038.

[2] A. Golzari Oskouei, M. Hashemzadeh, B. Asheghi, and M. A. Balafar, "CGFFCM: Cluster-weight and Group-local Feature-weight learning in Fuzzy C-Means clustering algorithm for color image segmentation," Applied Soft Computing, vol. 113, p. 108005, 2021/12/01/ 2021, doi: https://doi.org/10.1016/j.asoc.2021.108005.

[3] A. Golzari Oskouei and M. Hashemzadeh, "CGFFCM: A color image segmentation method based on cluster-weight and feature-weight learning," Software Impacts, vol. 11, p. 100228, 2022/02/01/ 2022, doi: https://doi.org/10.1016/j.simpa.2022.100228.

2) Please do not distribute the database or source codes to others without the authorization from Dr. Mahdi Hashemzadeh.

Authors’ Emails: arezounajafi@azaruniv.ac.ir (A.N. Moghadam), nasseraghazadeh[at]iyte.edu.tr (N. Aghazadeh), hashemzadeh[at]azaruniv.ac.ir (M. Hashemzadeh), and a.golzari[at]tabrizu.ac.ir (A. G. Oskouei)




