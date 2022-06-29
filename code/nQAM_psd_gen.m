% function nQAM_psd_gen()
nQAM = 256;

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

qamset = [(-sqrt(nQAM)/2+1):(sqrt(nQAM)/2)]*2-1;
qamset = repmat(qamset, sqrt(nQAM),1) + 1i*repmat(qamset', 1,sqrt(nQAM));
qamset = reshape(qamset, 1, [])./sqrt(2*(sqrt(nQAM)-1)^2);
digitalwaveformidx = randi([1 nQAM],T*bps/log2(nQAM),1);
digitalwaveformamp = qamset(digitalwaveformidx);
digitalwaveformamp = repmat(digitalwaveformamp', 1, fs/(bps/log2(nQAM)));
digitalwaveformamp = reshape(digitalwaveformamp',[],1);

ybaseband = digitalwaveformamp;
ybaseband = awgn(ybaseband,20,'measured');
% plot(real(ybaseband),imag(ybaseband),'.')
% scatterplot(ybaseband);


%% QAM
yQAM = yc.*ybaseband;

% plot(t(bpsspec).*1e6, real(yASK(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% ylim([-2,2])
% xlabel('Time (us)') 
% title(sprintf('ASK Waveform (fc = %dKHz, %dKbps)', round(fc/1e3), round(bps/1e3)))
% saveas(h,'ASK_Waveform','jpg');
% clf(h)

figure
yQAM_fft = psd_estimate(yQAM, avg_num, fs);
semilogy (fftaxis(fftwilderspec),yQAM_fft(fftwilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
hold on
semilogy(fftaxis(fftwilderspec), nQAM_psd(nQAM, bps, fc, fftaxis(fftwilderspec)*1e6) ,'-' ,...
    'LineWidth',1, 'color','red')

ylim([1e-8,1e-3])
xlabel('Frequency (MHz)') 
title(sprintf('QAM Spectrum (fc = %dKHz, %dKbps)', round(fc/1e3), round(bps/1e3)))



% saveas(h,'ASK_Spectrum','jpg');
% clf(h)


