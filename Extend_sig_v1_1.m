function s_ext = Extend_sig_v1_1(s,ext_type,L)

%
% It makes use of wextend and then make it periodical outside the boundaries. 
%
%   The valid extension modes (MODE) are:
%     'zpd' zero extension.
%     'sym' (or 'symh') symmetric extension (half-point).
%     'symw' symmetric extension (whole-point).
%     'asym' (or 'asymh') antisymmetric extension (half-point).
%     'asymw' antisymmetric extension (whole-point).
%     'ppd' periodized extension (1).
%     'per' periodized extension (2):
%        If the signal length is odd, WEXTEND adds an extra-sample
%        equal to the last value on the right and performs extension
%        using the 'ppd' mode. Otherwise, 'per' reduces to 'ppd'.
%        The same kind of rule stands for images.
%
%
% Please refer to the wextend help for more details
% 
% Based on an original idea by Emanuele Papini


if nargin < 2, help Extend_sig_v1_1, end
if nargin < 2, ext_type='asymw'; end
if nargin < 3, L = length(s); end

figOn=1==1;

if size(s,1)>size(s,2)
    s=s';
end
if size(s,1)>1
    size_s=2;
else
    size_s=1;
end


s_ext_o = wextend(size_s,ext_type,s,L);

C=@(x,L) 1/2*(cos(x*pi/L-pi)+1);

ave=mean(s_ext_o);

if size_s==1
    if strcmp(ext_type,'per')
        s_ext=(s_ext_o-ave).*[C(1:L,L) ones(1,length(s)+1) C(-L:-1,L)]+ave;
    else
    s_ext=(s_ext_o-ave).*[C(1:L,L) ones(1,length(s)) C(-L:-1,L)]+ave;
    end
    if figOn
    fig=figure;
    plot(s_ext_o,'r')
    hold on
    plot(s_ext)
    plot(L+1:L+length(s),s,'g')
    plot(L+1,s_ext_o(L+1),'bx')
    plot(length(s_ext_o)-L-1,s_ext_o(end-L-1),'bx')
    hold off
    lg=legend(['Plain ''' ext_type ''' extension'],['Smart ''' ext_type ''' extension'],'Original signal');
    set(lg,'Interpreter','latex')
    set(gca,'fontsize', 32); 
    axis tight
    set(fig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %pause(0.5)
    end
else
    disp('Code missing')
    s_ext=[];
end

end