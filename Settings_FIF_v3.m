function options = Settings_FIF_v3(varargin)
% 
% Settings_FIF_v3 Constructs option structure for the algorithms
%
% EXAMPLES
%
% OPTIONS = Settings_FIF_v3 with no input arguments returns
%   setting structure with default values
%
% OPTIONS = Settings_FIF_v3('NAME1',VALUE1,'NAME2',VALUE2,...) creates a
%   solution options structure OPTIONS in which the named properties have
%   the specified values.  Any unspecified properties have default values.
%   It is sufficient to type only the leading characters that uniquely
%   identify the property.  Case is ignored for property names.
%
% OPTIONS = Settings_FIF_v3(OLDOPTIONS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTIONS.
%
%
% Settings_FIF_v3 PROPERTIES :
%
% GENERAL
% 
%   saveEnd          - (0) If value is 1 save workspace at the end
%
%   saveInter        - (0) If value is 1 save intermediate results every
%                      time a new IMF is generated.
%
%   verbose          - (1) 0 the method is silent, 1 normal, >1 loud, 
%
%       
%   plots            - (0) the algorithm does not produce plots,
%                      1 it produces them 
%   
%   saveplots        - (0) Set to 1 to save automatically each plot
%
%
% SPECIFIC PARAMETERS
%
%
%  delta       (0.001) Stopping criterion
%  ExtPoints   (3)     Number of extrema allowed in the remainder
%  NIMFs       (200)   Number of IMFs we want to produce, not counting
%                           the remainder
%  Xi          (1.6)   Parameter we use to tune the mask length
%  alpha       ('ave') Parameter used for the mask length computation.
%                           Allowed values 'ave', [0,100] or 'Almost_min'.
%                           If set to 'ave' the mask length equals 
%                           round(2*Xi*(length of the signal)/(number of extrema)).
%                           If set to 0 the mask length is proportional to the
%                           0 percentile of the distances between two subsequent extrema.
%                           If set to 100 then it is proportional to the 
%                           100 percentile of the distances between two subsequent extrema. 
%                           Finally if set to 'Almost_min' it is set to 
%                           a value close to the minimum (30-th percentile).
%                           
%  MaxInner    (200)   Maximum number of inner steps allowed.
%  MonotoneMaskLength (true) Boolean: if true when the algorithm compute a new mask length that is smaller or equal 
%                            to the previous one then automatically it increases to 1.1 of the previous mask length. 
%                            If false it allows smaller mask lengths.
%  NumSteps       (1)   Number of inner steps the algorithm is going to jump
%  
%
% ------------------------------------------------------
% EXAMPLE
%          
%   >> options = Settings_FIF_v3('delta',0.08,'NIMFs',5,'plots',1) 
%   >> IMF = FIF_v2_12(x,options)
%              
%  Executes algorithm FIF_v2_12 with delta = 0.08, stops after we have at most 5 IMFs and a trend, it produces plots.                            
% ------------------------------------------------------      
%
% See also FIF_v2_12 and MvFIF_v7
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
%
% (Ripped from sdpsettings.m by Johan Lufberg)


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help Settings_FIF_v3
    return;
end


Names = {
    % General   
    'saveEnd'
    'saveInter'
    'verbose'    
    'plots'
    'saveplots'   
             
    % FIF
    'delta'
    'ExtPoints'
    'NIMFs'
    'Xi'
    'alpha'
    'MaxInner'
    'MonotoneMaskLength'
    'NumSteps'
};

obsoletenames ={ % Use when options have become obsolete
};

[m,n] = size(Names);
names = lower(Names);

if (nargin>0) && isstruct(varargin{1})
    options = varargin{1};
    paramstart = 2;
else
    paramstart = 1;
    
    % General 
    
    options.saveEnd = 0;
    options.saveInter = 0;
    options.verbose = 1;    
    options.plots = 0.0;
    options.saveplots = 0;     
        
    % FIF
    options.delta = 0.001;
    options.ExtPoints=3;
    options.NIMFs=200;
    options.Xi=1.6;
    options.alpha='ave';
    options.MaxInner=200;
    options.MonotoneMaskLength=true;
    options.NumSteps=1;
end

i = paramstart;
% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;       % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string property name.', i));
        end
        
        lowArg = lower(arg);
        
        j_old = strmatch(lowArg,obsoletenames);
        if ~isempty(j_old)
            % For compability... No need yet
        end
        
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized property name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous property name ''%s'' ', arg);
                msg = [msg '(' deblank(Names{j(1)})];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names{k})];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;    % we expect a value next
    else
        eval(['options.' Names{j} '= arg;']);
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for property ''%s''.', arg));
end

end

