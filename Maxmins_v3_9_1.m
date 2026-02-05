function varargout = Maxmins_v3_9_1(f,tol,Derivative)
% v3_9 Based on version 3_8
% Minor revision: - now we use the option "Derivative" to decide which
%                   derivative of the signal to use for the extrema
%                   identification. If set to 0 the method will use the
%                   signal itself.
%                 - We fix the periodical boundaries. f is assumed
%                   periodical
%
% v3_8 Based on version 3_7
% Minor revision: fixed a little bug when df==0
%
% v3_7 Based on version 3_6
% Minor revision: we extend f as df in order to make possible to compute
%                 the relative derivative like df(i)/abs(f(i)) for every i=1:N
%
% v3_6 Based on version 3_5
% Minor revision: we consider relative differences to avoid having problems with small or big values in a signal
%
% v3_5 Based on version 3_4
% Minor revision: we consider only periodical extension
%
% v3_4 Based on version 3
% Minor revisions: 1) added for constant extention the checking for Mins and
%                     Maxs emptiness
%                  2) completed the code for the periodical case
%
% v3 is Based on Version 2.
% Modified the way zero-derivative regions are handled.
%
% Identify the maxima and minima of a signal f

N = length(f);

if nargin<3
    Derivative=0;
end

if size(f,1)>size(f,2)
    f=f.';
end

Maxs = zeros(1,N);
Mins = zeros(1,N);

df=diff(f);
for i=1:Derivative
    df = diff(df);
end

h = 1;

while h<N && abs(df(h)/f(h)) <= tol
    h=h+1;
end
if h==N
    if nargout<=1
        varargout{1}=[];
    elseif nargout==2
        varargout{1}=[];
        varargout{2}=[];
    end
    return
end

cmaxs=0;
cmins=0;

c = 0;

N_old=N;
f=[f f(2:h)];
df=diff([f f(2:h+1)]);

for i=1:Derivative
    df = diff(df);
end

N=N+h-1;
last_df=[];
for i=h:N-1-Derivative
    if   df(i)*df(i+1)/abs(f(i))^2 <= tol && df(i)*df(i+1)/abs(f(i))^2 >= -tol
        if df(i)/abs(f(i)) < -tol
            last_df=-1;
            posc = i;
        elseif df(i)/abs(f(i)) > tol
            last_df=+1;
            posc = i;
        else
            last_df=0;
            posc = i;
        end
        c = c + 1;
        if df(i+1)/abs(f(i)) < -tol
            if last_df==+1 || last_df==0
                cmaxs=cmaxs+1;                
                Maxs(cmaxs)=mod(posc+floor((c-1)/2)+1,N_old);
            end
            c=0;
        end
        if df(i+1)/abs(f(i)) > tol
            if last_df==-1 || last_df==0
                cmins=cmins+1;
                Mins(cmins)=mod(posc+floor((c-1)/2)+1,N_old);
            end
            c=0;
        end
        
    end
    if   df(i)*df(i+1)/abs(f(i))^2 < -tol
        if df(i)/abs(f(i)) < -tol && df(i+1)/abs(f(i)) > tol
            cmins=cmins+1;
            Mins(cmins)=mod(i+1,N_old);
            if Mins(cmins)==0
                Mins(cmins)=1;
            end
            last_df=-1;
        elseif df(i)/abs(f(i)) > tol && df(i+1)/abs(f(i)) < -tol
            cmaxs=cmaxs+1;
            Maxs(cmaxs)=mod(i+1,N_old);
            if Maxs(cmaxs)==0
                Maxs(cmaxs)=1;
            end
            last_df=+1;
        end
    end
end
if c > 0
    %         % we deal with the boundary
    %         df_0=f(N)-f(1);
    %         if df_0==0
    %             if Initial_df < 0
    %                 if last_df==+1
    %                     cmaxs=cmaxs+1;
    %                     Maxs(cmaxs)=mod(posc+floor((c+cIn-1)/2)+1,N);
    %                 end
    %             elseif Initial_df > 0
    %                 if last_df==-1
    %                     cmins=cmins+1;
    %                     Mins(cmins)=mod(posc+floor((c+cIn-1)/2)+1,N);
    %                 end
    %             end
    %         else
    %             disp('Code missing!')
    %         end
    if cmins>0 && Mins(cmins)==0
        Mins(cmins)=N
    end
    if cmaxs>0 && Maxs(cmaxs)==0
        Maxs(cmaxs)=N
    end
end

Maxs=Maxs(1:cmaxs);
Mins=Mins(1:cmins);
%Maxs=Maxs(Maxs<=N_old);
%Mins=Mins(Mins<=N_old);
maxmins=sort([Maxs Mins]);
%     disp('Code to be completed')
%     if isempty(maxmins)
%         maxmins = 1;
%     else
%         if maxmins(1)~=1 && maxmins(end)~=N
%             if (f(maxmins(end)) > f(maxmins(end)+1) && f(maxmins(1)) > f(maxmins(1)-1)) || (f(maxmins(end)) < f(maxmins(end)+1) && f(maxmins(1)) < f(maxmins(1)-1))
%                 maxmins=[1 maxmins];
%             end
%         end
%     end

if sum(Maxs==0)>0
    Maxs(Maxs==0)=1;
end
if sum(Mins==0)>0
    Mins(Mins==0)=1;
end
if sum(maxmins==0)>0
    maxmins(maxmins==0)=1;
end

if nargout<=1
    varargout{1}=maxmins;
elseif nargout==2
    varargout{1}=Maxs;
    varargout{2}=Mins;
end

end