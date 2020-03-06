# FIF
Fast Iterative Filtering for the decompostion of non-stationary signals

Please refer to "Example_v7.m" and "Example_real_life_v5.m" for examples of how to use the code.

It is based on FFT, which makes FIF to be really fast.
This implies that it is required a periodical extension at the boundaries.

To overcome this limitation we can preextend the signal under investigation.
We do it thanks to the function "Extend_sig_v1_1.m". 
See "Example_real_life_v5.m" for an example of application.


Please cite our works: 

A. Cicone, J. Liu, H. Zhou. "Adaptive Local Iterative Filtering 
for Signal Decomposition and Instantaneous Frequency analysis". 
Applied and Computational Harmonic Analysis, Volume 41, Issue 2, 
September 2016, Pages 384-411. doi:10.1016/j.acha.2016.03.001 
Arxiv http://arxiv.org/abs/1411.6051

 A. Cicone, H. Zhou. "Numerical Analysis for Iterative Filtering with
 New Efficient Implementations Based on FFT"
 ArXiv http://arxiv.org/abs/1802.01359

 A. Cicone. "Iterative Filtering as a direct method for the decomposition 
 of nonstationary signals". Numerical Algorithms, Volume 373, 2020,  112248. 
 doi: 10.1007/s11075-019-00838-z
 ArXiv http://arxiv.org/abs/1811.03536
 
 
