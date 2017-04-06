function logarea_5_5 = highFreqPower(y)
%-------------------------------------------------------------------------------
% Computes the power in upper fifth of sampled frequencies of y
%-------------------------------------------------------------------------------
% HCTSA FEATURE: SP_Summaries_welch_rect.logarea_5_5
%-------------------------------------------------------------------------------
% Can reproduce as:
% out = SP_Summaries(y,'welch','rect',[],0); -> upperFifth = out.logarea_5_5;
%-------------------------------------------------------------------------------
% Reduced (faster) version below
%-------------------------------------------------------------------------------

psdMeth = 'welch';
nf = []; % num frequency components

% ------------------------------------------------------------------------------
% Check inputs, set defaults:
% ------------------------------------------------------------------------------
if size(y,2) > size(y,1);
    y = y'; % Time series must be a column vector
end
Ny = length(y); % time-series length
window = rectwin(Ny);

% ------------------------------------------------------------------------------
% Compute the Fourier Transform
% ------------------------------------------------------------------------------
% Welch power spectral density estimate:
Fs = 1; % sampling frequency
N = 2^nextpow2(Ny);
[S, f] = pwelch(y,window,[],N,Fs);
w = 2*pi*f'; % angular frequency
S = S/(2*pi); % adjust so that area remains normalized in angular frequency space

if ~any(isfinite(S)) % no finite values in the power spectrum
    % This time series must be really weird -- return NaN (unsuitable operation)...
    warning('NaN in power spectrum? A weird time series.');
    out = NaN; return
end

% Ensure both w and S are row vectors:
if size(S,1) > size(S,2)
    S = S';
end
if size(w,1) > size(w,2)
    w = w';
end

N = length(S); % = length(w)
dw = w(2) - w(1); % spacing increment in w

%-------------------------------------------------------------------------------
% Area (in top fifth of frequency bands)
%-------------------------------------------------------------------------------
split = buffer(S,floor(N/5));
if size(split,2) > 5,
    split = split(:,1:5);
end
logarea_5_5 = sum(log(split(:,5)))*dw;

end
