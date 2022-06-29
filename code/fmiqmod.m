function y = fmiqmod(x,Fc,Fs,freqdev,varargin)
%FMMOD Frequency modulation.
%   Y = FMMOD(X,Fc,Fs,FREQDEV) uses the message signal X to modulate the
%   carrier frequency Fc (Hz) and sample frequency Fs (Hz),  where Fs >
%   2*Fc. FREQDEV (Hz) is the frequency deviation of the modulated signal.
%
%   Y = FMMOD(X,Fc,Fs,FREQDEV,INI_PHASE) specifies the initial phase of
%   the modulation.
%
%   See also FMDEMOD, PMMOD, PMDEMOD.

%    Copyright 1996-2018 The MathWorks, Inc.

narginchk(1,5);

if(~isreal(x))
    error(message('comm:fmmmod:NonRealX'));
end

if(~isreal(Fs) || ~isscalar(Fs) || Fs<=0 )
    error(message('comm:fmmod:InvalidFs'));
end

if(~isreal(Fc) || ~isscalar(Fc) || Fc<0 )
    error(message('comm:fmmod:InvalidFc'));
end

if(~isreal(freqdev) || ~isscalar(freqdev) || freqdev<=0)
    error(message('comm:fmmod:InvalidFreqdev'));
end

% check that Fs must be greater than 2*Fc
if(Fs<2*Fc)
    error(message('comm:fmmod:FsLessThan2Fc'));
end

ini_phase = 0;
if(nargin == 5)
    ini_phase = varargin{1};
    if(isempty(ini_phase) || (isStringScalar(ini_phase) && strlength(ini_phase)==0))
        ini_phase = 0;
    end
    if(~isreal(ini_phase) || ~isscalar(ini_phase) )
        error(message('comm:fmmod:InvalidIni_Phase'));
    end
end
    

% --- Assure that X, if one dimensional, has the correct orientation --- %
len = size(x,1);
if(len == 1)
    x = x(:);
end
   
t = (0:1/Fs:((size(x,1)-1)/Fs))';
t = t(:,ones(1,size(x,2)));

int_x = cumsum(x)/Fs;
a = cos(2*pi*freqdev*int_x + ini_phase);   
b = sin(2*pi*freqdev*int_x + ini_phase);   
y = a+1i*b;
y = y.*exp(1i*2*pi*t*Fc);
% --- restore the output signal to the original orientation --- %
if(len == 1)
    y = y';
end
    
% EOF
