function ac34 = AC34(y)
%-------------------------------------------------------------------------------
% Computes the autocorrelation of time-series y, at a time-lag of 34 samples
% (note: therefore depends highly on sampling rate of data)
%-------------------------------------------------------------------------------
% Can reproduce as:
% out = CO_AutoCorr(x,34,'Fourier');
%-------------------------------------------------------------------------------

nFFT = 2^(nextpow2(length(y))+1);
F = fft(y-mean(y),nFFT);
F = F.*conj(F);
acf = ifft(F);
acf = acf./acf(1); % Normalize
acf = real(acf);

acf = acf(1:length(y));

if isempty(tau) % return the full function
    out = acf;
else % return a specific set of values
    out = zeros(length(tau),1);
    for i = 1:length(tau)
        if (tau(i) > length(acf)-1) || (tau(i) < 0)
            out(i) = NaN;
        else
            out(i) = acf(tau(i)+1);
        end
    end
end

end
