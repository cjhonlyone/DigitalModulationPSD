function plotPSD()

fs = 4e6;
fc = 0.1e6;
T = 480e-3;
t = (0:1/fs:T-1/fs)';
bps = 10e3;
avg_num = 10;
yc = exp(1i*2*pi*t*fc);

dF = 10e3;
N = length(t);
idx = [0:N-1];

Nfft = N/avg_num;
fftaxis = ([0:Nfft-1]-Nfft/2)/Nfft*fs/1e6;
fftwilderspec = idx(floor((fc-5*dF)/fs*Nfft+Nfft/2):floor((fc+5*dF)/fs*Nfft+Nfft/2)-1);
ffthalfwilderspec = idx(floor((fc)/fs*Nfft+Nfft/2):floor((fc+2*dF)/fs*Nfft+Nfft/2)-1);

figure
for nASK = [2,4,8]
semilogy(fftaxis(ffthalfwilderspec), nASK_psd(nASK, bps, fc, fftaxis(ffthalfwilderspec)*1e6) ,'-' ,...
    'LineWidth',1)
hold on 
end
ylim([1e-8,1e-3])

figure
for nPSK = [4,8,16]
semilogy(fftaxis(ffthalfwilderspec), nPSK_psd(nPSK, bps, fc, fftaxis(ffthalfwilderspec)*1e6) ,'-' ,...
    'LineWidth',1)
hold on 
end
ylim([1e-8,1e-3])

figure
for nQAM = [4,16,64]
semilogy(fftaxis(ffthalfwilderspec), nQAM_psd(nQAM, bps, fc, fftaxis(ffthalfwilderspec)*1e6) ,'-' ,...
    'LineWidth',1)
hold on 
end
ylim([1e-8,1e-3])