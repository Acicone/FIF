function [IMF,stats] = FIF_v2_12(f,options,M)


%
%  function IMF = FIF_v2_12(f,options,M)
%
% It generates the decomposition of the signal f :
%
%  f = IMF(1,:) + IMF(2,:) + ... + IMF(K, :)
%
% where the last row in the matrix IMF is the trend and the other rows
% are actual IMFs
%
%                                Inputs
%
%   f         Signal to be decomposed
%
%   options    Structure, generated using function Settings_FIF_v3, containing
%              all the parameters needed in the various algorithms
%
%   M         Mask length values for each Inner Loop
%
%                               Output
%
%   IMF       Matrices containg in row i the i-th IMF. The last row
%              contains the remainder-trend.
%
%   stats     Statistics regarding the IMFs
%               logM      Mask length values used for each IMF
%               posF      position of the first minimum in the filter DFT which is forced to become zero
%               valF      filter DFT first minimum value before the downward shift
%
%   See also SETTINGS_FIF_V3, GETMASK_V1, PLOT_IMF_V10.
%
%  Please cite:
%
%  A. Cicone, H. Zhou. "Numerical Analysis for Iterative Filtering with 
%  New Efficient Implementations Based on FFT". Numerische Mathematik, 2020. 
%  doi: 10.1007/s00211-020-01165-5
%  ArXiv http://arxiv.org/abs/1802.01359
%
%  A. Cicone. 'Iterative Filtering as a direct method for the decomposition
%  of nonstationary signals'. Numerical Algorithms, Volume 373, 2020,  112248.
%  doi: 10.1007/s11075-019-00838-z
%  ArXiv http://arxiv.org/abs/1811.03536


%% we deal with the input

tol=10^-12;

if nargin < 1,  help FIF_v2_12; return; end
if nargin < 2, options = Settings_FIF_v3; end
if nargin < 3, M = []; end

FigCol = 'ckmygr'; % Plot Colors
N = length(f);
if size(f,1)>size(f,2)
    f = f.';
end
if size(f,1)>1
    disp('Wrong dataset, the signal must be a single row vector')
    disp('If you have a multivariate signal you can try the MFIF code')
    disp('If, instead, is a multidimensional data set, plese try the FIF2 code')
    IMF=[];
    stats=[];
    return
end
IMF =zeros(options.NIMFs,N);

nameFile=sprintf('%1.0d',sum(round(clock*1000)));

Norm1f=norm(f,1); % to avoid dealing with way too small values
f=f/Norm1f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Iterative Filtering                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('prefixed_double_filter','MM');
%load('TriangularFilter','MM');
%load('New_doubly_convolved_filter','MM');


% k = length(maxmins);
% diffMaxmins=diff(maxmins);
%% Create a signal without zero regions and compute the number of extrema
f_pp=f;
f_pp(abs(f)<=tol)=[];
if isempty(f_pp)
    disp('Signal too small')
    IMF=[];
    stats=[];
    return
end

N0=length(f_pp);
pos=1:N0+10;% to include potential extrema at the boundaries we extend periodicaly of 10 points the signal
Mins=pos(islocalmin([f_pp f_pp(1:10)], 'FlatSelection', 'center'));
Mins(Mins>N0)=[];
Maxs=pos(islocalmax([f_pp f_pp(1:10)], 'FlatSelection', 'center'));
Maxs(Maxs>N0)=[];
maxmins_pp=sort([Mins Maxs]);
diffMaxmins_pp=diff(maxmins_pp);
N_pp=length(f_pp);
k_pp = length(maxmins_pp);

countIMFs=0;

while countIMFs < options.NIMFs && k_pp>=options.ExtPoints
    countIMFs=countIMFs+1;
    
    SD=1;
    h=f;
    
    if isempty(M) || length(M)<countIMFs
        
        if isa(options.alpha,'char')
            if strcmp(options.alpha,'ave') % Using an average mask length
                m = 2*round(N_pp/k_pp*options.Xi);
            elseif strcmp(options.alpha,'Almost_min') % Using an almost min mask length
                if 2*round(options.Xi*prctile(diffMaxmins_pp,30))<2*round(N_pp/k_pp*options.Xi)
                    m = 2*round(options.Xi*prctile(diffMaxmins_pp,30));
                else
                    m = 2*round(N_pp/k_pp*options.Xi);
                end
                elseif strcmp(options.alpha,'Median') % Using a median mask length                
                    m = 2*round(options.Xi*median(diffMaxmins_pp));   
            else
                disp(' Value of alpha not recognized')
                return
            end
        else % using a fixed value alpha
            m = 2*round(options.Xi*prctile(diffMaxmins_pp,options.alpha));
            %             if not(prctile(diffMaxmins_pp,options.alpha)==200)
            %             plot([0 10 20 30 40 50 60 70 80 90 100], prctile(diffMaxmins_pp,[0 10 20 30 40 50 60 70 80 90 100]))
            %             end
        end
        if countIMFs>1
            if m<=stats(countIMFs-1).logM
                if options.verbose>0
                    fprintf('Warning mask length is decreasing at step %1d. ',countIMFs)
                end
                if options.MonotoneMaskLength==true
                    m=ceil(stats(countIMFs-1).logM*1.1);
                    if options.verbose>0
                        fprintf('The old mask length is %1d whereas the new one is forced to be %1d.\n',stats(countIMFs-1).logM,ceil(stats(countIMFs-1).logM*1.1))
                    end
                else
                    if options.verbose>0
                        fprintf('The old mask length is %1d whereas the new one is %1d.\n',stats(countIMFs-1).logM,m)
                    end
                end
            end
            
            
            
        end
    else
        m=M(countIMFs);
    end
    
    inStepN=0;
    if options.verbose>0
        fprintf('\n IMF # %1.0d   -   # Extreme points %5.0d\n',countIMFs,k_pp)
        fprintf('\n  step #            SD             Mask length \n\n')
    end
    stats(countIMFs).logM=m;
    a = get_mask_v1_1(MM,m,options.verbose,tol);
    ExtendSig=1==0;
    if N < length(a) % we need to extend the signal
        ExtendSig=1==1;
        Nxs=ceil(length(a)/N);
        N_old=N;
        if rem(Nxs,2)==0
            Nxs=Nxs+1;
        end
        h_n=[];
        for ii=1:Nxs
            h_n=[h_n h];
        end
        h=h_n;
        N=Nxs*N;
    end
    
    Nza=N-length(a);
    if rem(Nza,2)==0
        a = [zeros(1,Nza/2) a zeros(1,Nza/2)];
        ifftA=real(fft([a((length(a)-1)/2+1:end) a(1:(length(a)-1)/2)]));
        % figure,plot(circshift(a,(length(a)-1)/2+1)-ifft(real(fft(circshift(a,(length(a)-1)/2+1)))),'r')
    else
        a = [zeros(1,(Nza-1)/2) a zeros(1,(Nza-1)/2+1)];
        %csA=circshift(a,(length(a))/2+1);
        ifftA=real(fft([a((length(a))/2:end) a(1:(length(a))/2-1)]));
        % figure,plot(circshift(a,(length(a))/2+1)-ifft(real(fft(circshift(a,(length(a))/2+1)))),'r')
    end
    if options.plots>0 %&& rem(inStepN,5)==0
        if gcf > 30
            close all
        end
        figN=figure;
        set(figN,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    fftH=fft(h);
    fft_h_new=fftH;
    
    %% we compensate for the aliasing effect in the DFT of the filter
    % we look for the first minimum in ifftA
    
    posF=find((diff(ifftA)>0)~=0,1,'first');
    
    stats(countIMFs).posF = posF;
    stats(countIMFs).valF = ifftA(posF);
    
    ifftA=ifftA-ifftA(posF);
    ifftA(ifftA<0)=0;
    
    %plot([ifftA(end/2+2:end) ifftA(1:end/2+1)])
    %plot([a((length(a)-1)/2+1:end) a(1:(length(a)-1)/2)]-ifft(ifftA))
    % plot(ifft(ifftA))
    if options.plots>=1
        figMask=figure;
        figRem=figure;
        set(figMask,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    while SD>options.delta && inStepN < options.MaxInner
        inStepN=inStepN+options.NumSteps;
        
        fft_h_old=(1-ifftA).^(inStepN-1).*fftH;
        fft_h_new=(1-ifftA).^inStepN.*fftH;
        
        %%%%%%%%%%%%%%%% Updating stopping criterium %%%%%%%%%%%%%%%%%
        
        SD=norm(fft_h_new-fft_h_old)^2/norm(fft_h_old)^2;
        
        %%%%%%%%%%%%%%%%%% generating f_n %%%%%%%%%%%%%%%%%%
        
        if options.verbose>0
            fprintf('    %2.0d      %1.40f          %2.0d\n',inStepN,SD,m)
        end
        
        if options.plots>=1  && rem(inStepN,2)==0
            figure(figMask)
            title(['IMF ' num2str(countIMFs) ' step # ' num2str(inStepN) ])
            plot(ifft(fft_h_new)*Norm1f,'linewidth',2)
            figure(figRem)
            title(['Remainder after IMF ' num2str(countIMFs) ' step # ' num2str(inStepN) ])
            plot((f-ifft(fft_h_new))*Norm1f,'linewidth',2)
            
            pause(0.01)
        end
        
    end
    
    h=ifft(fft_h_new);
    
    if ExtendSig % we reduce the signal
        N=N_old;
        h=h(N*(Nxs-1)/2+1:N*((Nxs-1)/2+1));
    end
    if inStepN >= options.MaxInner
        disp('Max # of inner steps reached')
        %return
    end
    stats(countIMFs).inStepN=inStepN;
    IMF(countIMFs,:) = h;
    
    if options.plots>0
        figure(100)
        periodogram(f)
        title(['Signal before removing IMF = ' num2str(countIMFs)])
    end
    f=f-h;
    if options.plots>0
        figure(101)
        periodogram(f)
        title(['Signal after removing IMF = ' num2str(countIMFs)])
        pause
    end
    
    %% Create a signal without zero regions and compute the number of extrema
    f_pp=f;
    f_pp(abs(f)<=tol)=[];
    if isempty(f_pp)
        break
    end
    N0=length(f_pp);
    pos=1:N0+10;% to include potential extrema at the boundaries we extend periodicaly of 10 points the signal
    Mins=pos(islocalmin([f_pp f_pp(1:10)], 'FlatSelection', 'center'));
    Mins(Mins>N0)=[];
    Maxs=pos(islocalmax([f_pp f_pp(1:10)], 'FlatSelection', 'center'));
    Maxs(Maxs>N0)=[];
    maxmins_pp=sort([Mins Maxs]);
    if isempty(maxmins_pp)
        break
    end
    diffMaxmins_pp=diff(maxmins_pp);
    N_pp=length(f_pp);
    k_pp = length(maxmins_pp);
    
    if options.saveInter==1
        save([nameFile '_intermediate_FIF_v2_12.mat'],'IMF','f','stats','-v7.3');
    end
end %end of while

IMF = [IMF(1:countIMFs,:); f];

IMF=IMF*Norm1f; % we scale back to the original values

if options.plots>=1
    if gcf > 30
        close all
    end
    figN=plot_imf_v10(IMF,1:N);
    for i=1:length(figN)
        set(figN(i),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        if options.saveplots>0
            saveas(figN(i),[nameFile '_IMFs'], 'fig')
            saveas(figN(i),[nameFile '_IMFs'], 'epsc')
            saveas(figN(i),[nameFile '_IMFs'], 'png')
        end
        
    end
end

if options.saveEnd == 1
    save([ 'Final_' nameFile '_FIF_v2_12.mat'],'IMF','stats','-v7.3');
end

end


%% Auxiliar functions

function a=get_mask_v1_1(y,k,verbose,tol)
%
% Rescale the mask y so that its length becomes 2*k+1.
% k could be an integer or not an integer.
% y is the area under the curve for each bar

n=length(y);
m=(n-1)/2;

if k<=m % The prefixed filter contains enough points
    
    if mod(k,1)==0     % if the mask_length is an integer
        
        a=zeros(1,2*k+1);
        
        for i=1:2*k+1
            s=(i-1)*(2*m+1)/(2*k+1)+1;
            t=i*(2*m+1)/(2*k+1);
            
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            %t2=ceil(t)-t;
            
            if floor(t)<1
                disp('Ops')
            end
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        
    else   % if the mask length is not an integer
        new_k=floor(k);
        extra = k-new_k;
        c=(2*m+1)/(2*new_k+1+2*extra);
        
        a=zeros(1,2*new_k+3);
        
        t=extra*c+1;
        t1=t-floor(t);
        %t2=ceil(t)-t;
        if k<0
            disp('Ops')
            a=[];
            return
        end
        a(1)=sum(y(1:floor(t)))+t1*y(floor(t));
        
        for i=2:2*new_k+2
            s=extra*c+(i-2)*c+1;
            t=extra*c+(i-1)*c;
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            
            
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        t2=ceil(t)-t;
        
        a(2*new_k+3)=sum(y(ceil(t):n))+t2*y(ceil(t));
    end
else % We need a filter with more points than MM, we use interpolation
    dx=0.01;
    % we assume that MM has a dx = 0.01, if m = 6200 it correspond to a
    % filter of length 62*2 in the physical space
    f=y/dx; % function we need to interpolate
    dy=m*dx/k;
    b=interp1(0:m,f(m+1:2*m+1),0:m/k:m);
    if size(b,1)>size(b,2)
        b=b.';
    end
    if size(b,1)>1
        fprintf('\n\nError!')
        disp('The provided mask is not a vector!!')
        a=[];
        return
    end
    a=[fliplr(b(2:end)) b]*dy;
    if abs(norm(a,1)-1)>tol
        if verbose>0
            fprintf('\n\n Warning!\n\n')
            fprintf(' Area under the mask equals %2.20f\n',norm(a,1))
            fprintf(' it should be equal to 1\n We rescale it using its norm 1\n\n')
        end
        a=a/norm(a,1);
    end
end

end

