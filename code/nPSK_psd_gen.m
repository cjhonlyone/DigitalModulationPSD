% function nPSK_psd_gen()
nPSK = 8;

fs = 4e6;
fc = 0.1e6;
T = 480e-3;
t = (0:1/fs:T-1/fs)';
bps = 10e3;
avg_num = 50;
yc = exp(1i*2*pi*t*fc);

dF = 10e3;
N = length(t);
idx = [0:N-1];

Nfft = N/avg_num;
fftaxis = ([0:Nfft-1]-Nfft/2)/Nfft*fs/1e6;
fftwilderspec = idx(floor((fc-5*dF)/fs*Nfft+Nfft/2):floor((fc+5*dF)/fs*Nfft+Nfft/2)-1);

pskset = exp(1i*([1:nPSK]-1)./nPSK.*2*pi);
digitalwaveformidx = randi([1 nPSK],T*bps/log2(nPSK),1);
digitalwaveformamp = pskset(digitalwaveformidx);
digitalwaveformamp = repmat(digitalwaveformamp', 1, fs/(bps/log2(nPSK)));
digitalwaveformamp = reshape(digitalwaveformamp',[],1);


ybaseband = digitalwaveformamp;
ybaseband = awgn(ybaseband,20,'measured');
% plot(real(ybaseband),imag(ybaseband),'.')
scatterplot(ybaseband);
%% PSK
yPSK = yc.*ybaseband;

% plot(t(bpsspec).*1e6, real(yASK(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% ylim([-2,2])
% xlabel('Time (us)') 
% title(sprintf('ASK Waveform (fc = %dKHz, %dKbps)', round(fc/1e3), round(bps/1e3)))
% saveas(h,'ASK_Waveform','jpg');
% clf(h)




figure
yPSK_fft = psd_estimate(yPSK, avg_num, fs);
semilogy (fftaxis(fftwilderspec),yPSK_fft(fftwilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
hold on
semilogy(fftaxis(fftwilderspec), nPSK_psd(nPSK, bps, fc, fftaxis(fftwilderspec)*1e6) ,'-' ,...
    'LineWidth',1, 'color','red')

ylim([1e-8,1e-3])
xlabel('Frequency (MHz)') 
title(sprintf('PSK Spectrum (fc = %dKHz, %dKbps)', round(fc/1e3), round(bps/1e3)))



% saveas(h,'ASK_Spectrum','jpg');
% clf(h)


