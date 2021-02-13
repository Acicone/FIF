% Test Example
%
%  Example 1 page 16
%
%  Ref: A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for 
%  Signal Decomposition and Instantaneous Frequency analysis'. Applied and 
%  Computational Harmonic Analysis, Volume 41, Issue 2, September 2016, 
%  Pages 384-411. doi:10.1016/j.acha.2016.03.001
%  ArXiv http://arxiv.org/abs/1411.6051
%
%  Please cite: 
%
%  A. Cicone, H. Zhou. "Numerical Analysis for Iterative Filtering with 
%  New Efficient Implementations Based on FFT". Numerische Mathematik, 147 (1), pages 1-28, 2021. 
%  doi: 10.1007/s00211-020-01165-5
%  ArXiv http://arxiv.org/abs/1802.01359
%
%  A. Cicone. 'Iterative Filtering as a direct method for the decomposition 
%  of nonstationary signals'. Numerical Algorithms, Volume 373, 2020,  112248. 
%  doi: 10.1007/s11075-019-00838-z
%  ArXiv http://arxiv.org/abs/1811.03536
%

dt=0.001;

t=0:0.001:1;

x=(2*(t-0.5).^2+0.2).*sin(20*pi*t+0.2*cos(40*pi*t));

y=4*(t-0.5).^2;

z=x+y+1;

plot_imf_v10([x;y+ones(1,length(t))],t,2)
title('Ground Truth')

figure
plot(t,z,'b','linewidth',2)
set(gca,'fontsize', 20);
title('Signal')

%% FIF with standard settings

[IMF0,logM] = FIF_v2_12(z);

plot_imf_v10(IMF0,[],4)

%% FIF with some tuning of the parameters

opts=Settings_FIF_v3('Xi',2);

[IMF,logM] = FIF_v2_12(z,opts);

plot_imf_v10(IMF,[],2)

%%
plot_imf_v11(IMF,[x;y+ones(1,length(t))],t,2,[],[],[],[],'IMFs','Ground truth')


