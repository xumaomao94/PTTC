# Probabilistic Tensor Train Completion based on Gaussian-product-Gamma Prior, with Automatic Rank Estimation
## Functions

- demo_rank_accuracy_test.m

Run this demo to test the accuracy of the rank estimation. It randomly generate noise 10 times and tests VITTC for different settings.

- demo_image_completion.m

Run this demo to test VITTC on image completion. The provided image in "TestImages/missing_rate80_1.mat" is the 'jellybeans' image.

- f_tensorfolding/

Includes functions for ket-folding and the improved ket-folding in the reference.

- f_vittc/

Includes functions realizing the variational inference algorithms proposed for the probabilistic tensor train model.

- f_datageneration/

Includes functions that generate synthetic data.

- f_perfevaluate/

Includes functions that evaluate the performance of the recovered tensor.

- TestImages/

A 'jellybeans' image, and a mask with 80% entries missing

- experiment results/

A folder used for storing results.

### Reference

[1] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2020). Learning tensor train representation with automatic rank determination from incomplete noisy data. arXiv preprint arXiv:2010.06564.

[2] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2021, December). Overfitting Avoidance in Tensor Train Factorization and Completion: Prior Analysis and Inference. In 2021 IEEE International Conference on Data Mining (ICDM) (pp. 1439-1444). IEEE.

[3] Xu, L., Cheng, L., Wong, N., & Wu, Y. C. (2023). Tensor train factorization under noisy and incomplete data with automatic rank estimation. Pattern Recognition, 141, 109650.
