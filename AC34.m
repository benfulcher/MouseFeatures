function ac34 = AC34(y)
%-------------------------------------------------------------------------------
% Computes the autocorrelation of time-series y, at a time-lag of 34 samples
% (note: therefore depends highly on sampling rate of data)
%-------------------------------------------------------------------------------
% HCTSA FEATURE: "AC_34"
% Can reproduce as:
% out = CO_AutoCorr(x,34,'Fourier');
% Reduced code below
%-------------------------------------------------------------------------------

nFFT = 2^(nextpow2(length(y))+1);
F = fft(y-mean(y),nFFT);
F = F.*conj(F);
acf = ifft(F);
acf = acf./acf(1); % Normalize
acf = real(acf);

ac34 = acf(34+1);

end
