function psd = psd_estimate(x, avg, fs)
N = length(x);
Nfft = N/avg;
psd = zeros(Nfft,1);
for i = 1:avg
    psd = abs(fftshift(fft(x((i-1)*Nfft+1:i*Nfft), Nfft))).^2 + psd;
end
psd = psd/N/fs;