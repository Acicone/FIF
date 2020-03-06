% Test Example
%
% Example 7 page 25 - Length of the day dataset
%
%  A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for
%  Signal Decomposition and Instantaneous Frequency analysis'. Applied and
%  Computational Harmonic Analysis, Volume 41, Issue 2, September 2016,
%  Pages 384-411. doi:10.1016/j.acha.2016.03.001
%  ArXiv http://arxiv.org/abs/1411.6051
%
%  Please cite: 
%
%  A. Cicone, H. Zhou. 'Numerical Analysis for Iterative Filtering with
%  New Efficient Implementations Based on FFT'
%  ArXiv http://arxiv.org/abs/1802.01359
%
%  A. Cicone. 'Iterative Filtering as a direct method for the decomposition 
%  of nonstationary signals'. Numerical Algorithms, Volume 373, 2020,  112248. 
%  doi: 10.1007/s11075-019-00838-z
%  ArXiv http://arxiv.org/abs/1811.03536%
%
% dataset obtained from http://hpiers.obspm.fr/eoppc/eop/eopc04/eopc04.62-now
%
clear all
clc
load LengthOftheDay_LOD_ALIF_paper

plot(x)
title('Length of the Day')

%% Extension outside the boundaries to reduce the end effects in the decomposition

L=2*length(x);
x_ext = Extend_sig_v1_1(x,'asymw',L);

%%
opts=Settings_FIF_v3('delta',10^-2,'Xi',3,'NIMFs',4,'alpha','Almost_min');
clc
tic
[IMF_1,logM1] = FIF_v2_8(x,opts);%y,opts);
toc
tic
[IMF_2,logM2] = FIF_v2_8(x_ext,opts);%y,opts);
toc
%
plot_imf_v11(IMF_2(:,L+1:end-L),IMF_1,[],5,[],[],[],[],'pre-extended signal IMFs','original signal IMFs (end effects present)');





