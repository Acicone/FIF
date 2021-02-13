function s_ext = Extend_sig_v2(s,ext_type,L,figOn)

% s_ext = Extend_sig_v2(s,ext_type,L,figOn)
%
% It extends a signal left and right by means of wextend
% and then make it periodical outside the boundaries.
%
%  Inputs
%  s        = signal to be extended
%  ext_type = cell containing either one entry or two depending is we want
%             to have the same extension both left and right or if we prefer to have
%             different extensions. In this second case we insert first the left and
%             then that right extension option.
%             The valid extension modes are:
%             'zpd' zero extension.
%             'sym' (or 'symh') symmetric extension (half-point).
%             'symw' symmetric extension (whole-point).
%             'asym' (or 'asymh') antisymmetric extension (half-point).
%             WARNING: this option is going to do point reflection as if
%             the signal had an intermidiate point (i.e. the half point)
%             immediately outside the baundaries of value 0!!!
%             This can create big jumps at the bouandaries! Be careful
%             'asymw' antisymmetric extension (whole-point). This option is
%             going to do point reflection w.r.t. the actual bouandry value
%             of your signal.
%             'ppd' periodized extension (1).
%             'per' periodized extension (2):
%                   If the signal length is odd, WEXTEND adds an extra-sample
%                   equal to the last value on the right and performs extension
%                   using the 'ppd' mode. Otherwise, 'per' reduces to 'ppd'.
%                   The same kind of rule stands for images.
%              Examples: if we want to extend both sides symmetrically we
%              input ext_type with {'symw'}. If, instead, we want to have
%              left extension antisymmetrical and the right one periodical
%              we input {'asymw','per'}.
%              Default value is {'asymw'}
%  L         = extension length. Default value is the length of s.
%  figOn     = boolean which allows to switch on or off the plot showing the
%              extension produced by the function. Default value true
%
% Please refer to the wextend help for more details
%
% Based on an original idea by Emanuele Papini


if nargin < 2, help Extend_sig_v2, end
if nargin < 2, ext_type={'asymw'}; end
if nargin < 3, L = length(s); end
if nargin < 4, figOn=true; end

if size(s,1)>size(s,2)
    s=s';
end
if size(s,1)>1
    size_s=2;
else
    size_s=1;
end

L=round(L); % we want to have natural numbers
Ls=size(s,2);

if size_s>1
    disp('Code missing')
    s_ext=[];
    return
end
if size(ext_type,2)==1
    s_ext_o = wextend(size_s,ext_type{1},s,L);
else
    
    N_ext=floor(L/Ls);
    temp_s=s;
    for ii=1:N_ext
        if rem(ii,2)==1
            s_ext_o_1 = wextend(size_s,ext_type{1},temp_s,Ls);
            s_ext_o_2 = wextend(size_s,ext_type{2},temp_s,Ls);
        else
            s_ext_o_1 = wextend(size_s,ext_type{2},temp_s,Ls);
            s_ext_o_2 = wextend(size_s,ext_type{1},temp_s,Ls);
        end
        temp_s = [s_ext_o_1(1:end-Ls) s_ext_o_2(end-Ls+1:end)];
    end
    if rem(L,Ls)>0
        if rem(N_ext,2)==0
            s_ext_o_1 = wextend(size_s,ext_type{1},temp_s,rem(L,Ls));
            s_ext_o_2 = wextend(size_s,ext_type{2},temp_s,rem(L,Ls));
        else
            s_ext_o_1 = wextend(size_s,ext_type{2},temp_s,rem(L,Ls));
            s_ext_o_2 = wextend(size_s,ext_type{1},temp_s,rem(L,Ls));
        end
        temp_s = [s_ext_o_1(1:end-rem(L,Ls)) s_ext_o_2(end-rem(L,Ls)+1:end)];
    end
    s_ext_o = temp_s;
end
    C=@(x,L) 1/2*(cos(x*pi/L-pi)+1);
    
    ave=mean(s_ext_o);
    
    s_ext=(s_ext_o-ave).*[C(1:L,L) ones(1,length(s_ext_o)-2*L) C(-L:-1,L)]+ave;
    if figOn
        fig=figure;
        plot(s_ext_o,'r');
        hold on
        plot(s_ext,'b');
        plot(L+1:L+Ls,s,'g');
        plot(L+1,s_ext_o(L+1),'bx')
        plot(length(s_ext_o)-L,s_ext_o(end-L),'bx')
        hold off
        lg=legend('Plain extension','Smart extension','Original signal');
        set(lg,'Interpreter','latex')
        set(gca,'fontsize', 32);
        axis tight
        set(fig,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        %pause(0.5)
    end
    
    
end