function[h]=plot_imf_v11(imf1,imf2,T,NsubPltos,titLeg,h,TextSize,date,legText1,legText2)

% Input
% imf1 = First set of IMFs to be plotted
% imf2 = Second set of IMFs to be plotted
% T   = time axis
% titLeg = Title to be displayed on top to each figure
% h = handles to the figures
% TextSize = size of the axis and title
% 
% legText1 = legend text for the first IMF series
% legText2 = legend text for the second IMF series


[m,n]=size(imf1);

if nargin<3, T=1:n; end

if isempty(T), T=1:n; end

if nargin <4, NsubPltos=5; end

if isempty(NsubPltos), NsubPltos=5; end

if nargin <6, for j=1:ceil(m/NsubPltos), h(j)=figure; end, end

if nargin <7, TextSize=24; end

if isempty(TextSize), TextSize=24; end

if nargin <8, date=-1; end

LeftCorner=0.05;%0.2;
BottomCorner=0.06;
widthImage=0.9;%0.75;
DeltaYImages=0.01;
heightImage=(1-BottomCorner-BottomCorner/2)/NsubPltos-DeltaYImages;


if length(h)<ceil(m/NsubPltos)
    for i= length(h)+1:ceil(m/NsubPltos)
        h(i)=figure;
    end
end

for j=1:ceil(m/NsubPltos)
    figure(h(j));
    if j<=floor(m/NsubPltos)
        for i=1:NsubPltos
            subplot('Position',[LeftCorner BottomCorner+(NsubPltos-i)*(heightImage+DeltaYImages) widthImage heightImage]);
            I2=plot(T,imf2(i+(j-1)*NsubPltos,:),'r','LineWidth',2);
            hold on
            I1=plot(T,imf1(i+(j-1)*NsubPltos,:),'k','LineWidth',2);
            axis([T(1) T(end) -Inf Inf])
            set(gca,'fontsize', TextSize);
            if i==NsubPltos
                if date>=0
                    datetickzoom('x',date)
                end
            else
                set(gca,'XTick',[]);
            end            
            %pause(0.1)
        end
    else
        for i=1:rem(m,NsubPltos)
            subplot('Position',[LeftCorner BottomCorner+(NsubPltos-i)*(heightImage+DeltaYImages) widthImage heightImage]);
            I2=plot(T,imf2(i+(j-1)*NsubPltos,:),'r','LineWidth',2);
            hold on
            I1=plot(T,imf1(i+(j-1)*NsubPltos,:),'k','LineWidth',2); 
            axis([T(1) T(end) -Inf Inf])
            set(gca,'fontsize', TextSize);
            if i==rem(m,NsubPltos)
                if date>=0
                    datetickzoom('x',date)
                end
            else
                set(gca,'XTick',[]);
            end            
            %pause(0.1)
        end
    end
    if nargin > 8 && not(isempty(legText1)) && not(isempty(legText2))
        hh=legend([I1,I2],{legText1,legText2},'Location','none');
        set(hh,'Interpreter','latex')
    end
    set(h(j),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
end

if nargin >= 5  && not(isempty(titLeg))
    for j=1:ceil(m/NsubPltos)
        figure(h(j));
        set(gcf,'NextPlot','add');
        axes;
        if iscell(titLeg)
            if length(titLeg)>=j
                tt= title(titLeg{j},'FontSize',TextSize);
            else
                tt= title(titLeg{end},'FontSize',TextSize);
            end
        else
            tt= title(titLeg,'FontSize',TextSize);
        end
        set(tt,'Interpreter','latex')
        set(gca,'Visible','off');
        set(tt,'Visible','on');
    end
end


end

function [YoN,checkAve]=CheckForIMF(f)

maxmins = Maxmins(f,'No');
if length(maxmins)<2
    YoN=1==0;
    checkAve=[];
else
    if f(maxmins(1))>f(maxmins(2))
        checkAve1=sum(f(maxmins(1:2:end))>0)/length(maxmins(1:2:end));
        checkAve2=sum(f(maxmins(2:2:end))<0)/length(maxmins(1:2:end));
    else
        checkAve1=sum(f(maxmins(2:2:end))>0)/length(maxmins(2:2:end));
        checkAve2=sum(f(maxmins(1:2:end))<0)/length(maxmins(1:2:end));
    end
    
    if checkAve1 >= 0.5 && checkAve2 >= 0.5
        YoN=1==1;
    else
        YoN=1==0;
    end
    checkAve=[checkAve1 checkAve2];
end

end

function maxmins = Maxmins(f,extensionType)

if nargin == 1, extensionType = 'p'; end

N = length(f);
maxmins=zeros(1,N);
df = diff(f);


h = 1;
cIn=0;
if strcmp(extensionType,'p') && df(1) == 0 && df(end) == 0
    while df(h)==0
        cIn=cIn+1;
        h=h+1;
    end
end

c = 0;
cmaxmins=0;
for i=h:N-2
    if   df(i)*df(i+1) <= 0
        if df(i+1) == 0
            if c == 0
                posc = i;
            end
            c = c + 1;
        else
            if c > 0
                cmaxmins=cmaxmins+1;
                maxmins(cmaxmins)=posc+floor((c-1)/2)+1;
                c = 0;
            else
                cmaxmins=cmaxmins+1;
                maxmins(cmaxmins)=i+1;
            end
        end
    end
end
if c > 0
    cmaxmins=cmaxmins+1;
    maxmins(cmaxmins)=mod(posc+floor((c+cIn-1)/2)+1,N);
    if maxmins(cmaxmins)==0
        maxmins(cmaxmins)=N;
    end
end

maxmins=maxmins(1:cmaxmins);

if strcmp(extensionType,'p') % we deal with a periodic signal
    if isempty(maxmins)
        maxmins = 1;
    else
        if maxmins(1)~=1 && maxmins(end)~=N
            if (f(maxmins(end)) > f(maxmins(end)+1) && f(maxmins(1)) > f(maxmins(1)-1)) || (f(maxmins(end)) < f(maxmins(end)+1) && f(maxmins(1)) < f(maxmins(1)-1))
                maxmins=[1 maxmins];
            end
        end
    end
elseif strcmp(extensionType,'c')
    if isempty(maxmins)
        maxmins = [1, N];
    else
        if maxmins(1) ~= f(1) && maxmins(end) ~= f(end)
            maxmins = [1, maxmins, N];
        elseif f(maxmins(1)) ~= f(1)
            maxmins = [1, maxmins];
        elseif  f(maxmins(end)) ~= f(end)
            maxmins = [maxmins, N];
        end
    end
elseif strcmp(extensionType,'r')
    if isempty(maxmins)
        maxmins = [1, N];
    else
        if maxmins(1) ~= f(1) && maxmins(end) ~= f(end)
            maxmins = [1, maxmins, N];
        elseif f(maxmins(1)) ~= f(1)
            maxmins = [1, maxmins];
        elseif  f(maxmins(end)) ~= f(end)
            maxmins = [maxmins, N];
        end
    end
end

end
