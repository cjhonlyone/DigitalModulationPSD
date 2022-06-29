clear;clc;
% function [yBPSK_fft,yQPSK_fft,yOPSK_fft] = GenAMDSB_SSB

fs = 4e6;
fc = 0.1e6;

T = 4800e-3;
t = (0:1/fs:T-1/fs)';
dF = 10e3;
% bps = 1e3;
bps = 10e3;

avg_num = 50;

N = length(t);
freqaxis = ([0:N-1]-N/2)/N*fs/1e6;
Nfft = N/avg_num;
fftaxis = ([0:Nfft-1]-Nfft/2)/Nfft*fs/1e6;
idx = [0:N-1];

yc = exp(1i*2*pi*t*fc);

fc1 = fc-2.5*dF;
fc2 = fc+2.5*dF;
yc1 = exp(1i*2*pi*t*fc1);
yc2 = exp(1i*2*pi*t*fc2);

centerspec = idx(floor((fc-0.5*dF)/fs*N+N/2):floor((fc+0.5*dF)/fs*N+N/2)-1);
wilderspec = idx(floor((fc-5*dF)/fs*N+N/2):floor((fc+5*dF)/fs*N+N/2)-1);
halfwilderspec = idx(floor((fc)/fs*N+N/2):floor((fc+5*dF)/fs*N+N/2)-1);
zerospec = idx(floor((0-10*dF)/fs*N+N/2):floor((0+10*dF)/fs*N+N/2)-1);

fftwilderspec = idx(floor((fc-5*dF)/fs*Nfft+Nfft/2):floor((fc+5*dF)/fs*Nfft+Nfft/2)-1);
ffthalfwilderspec = idx(floor((fc)/fs*Nfft+Nfft/2):floor((fc+5*dF)/fs*Nfft+Nfft/2)-1);
Nbit = 10;
bpsspec = 1:floor(Nbit*fs/bps);
waveform = phased.LinearFMWaveform('PulseWidth',T,...
    'SweepBandwidth',dF,'PRF',1/T, 'SampleRate', fs);
analogwaveform = real(waveform());
% analogwaveform = cos(2*pi*dF*t);
% analogwaveform = zeros(N,1);
% for i = 1:100
% analogwaveform = cos(2*pi*i/100*dF*t)+analogwaveform;
% end

% [y,Fs] = audioread('funky-stereo.wav');

Nbit = 4;
bpsspec = 1:floor(Nbit*fs/bps);

% h=figure(1);
% clf(h)
% set(h,'visible','off');
% am_spectrum(analogwaveform, t, dF, fc, fs, bpsspec)


% 
% %% FM
% mf = 1;
% yfm = fmiqmod(analogwaveform,fc,fs,mf*dF);
% yfm_1 = yfm;
% 
% plot(t(bpsspec).*1e6, real(yfm_1(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% hold on
% plot(t(bpsspec).*1e6, 1.5*analogwaveform(bpsspec), 'LineWidth',1, 'color','black')
% ylim([-2,2])
% xlabel('Time (us)') 
% title(sprintf('FM Waveform (mf = %d)', mf))
% saveas(h,'FM-1_Waveform','jpg');
% clf(h)
% 
% yfm_fft = abs(fftshift(fft(yfm)));
% yfm_fft = yfm_fft./max(yfm_fft(centerspec))*abs(besselj(0,mf));
% stem (freqaxis(wilderspec),yfm_fft(wilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% xlabel('Frequency (MHz)') 
% title(sprintf('FM Spectrum (mf = %d)', mf))
% saveas(h,'FM-1_Spectrum','jpg');
% clf(h)
% 
% mf = 5;
% yfm = fmiqmod(analogwaveform,fc,fs,mf*dF);
% yfm_5 = yfm;
% 
% holdposition = h.OuterPosition;
% h.OuterPosition(3) = 2*h.OuterPosition(3);
% plot(t(bpsspec).*1e6, real(yfm_5(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% hold on
% plot(t(bpsspec).*1e6, 1.5*analogwaveform(bpsspec), 'LineWidth',1, 'color','black')
% ylim([-2,2])
% xlabel('Time (us)') 
% title(sprintf('FM Waveform (mf = %d)', mf))
% saveas(h,'FM-5_Waveform_wider','jpg');
% clf(h)
% h.OuterPosition = holdposition;
% 
% plot(t(bpsspec).*1e6, real(yfm_5(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% hold on
% plot(t(bpsspec).*1e6, 1.5*analogwaveform(bpsspec), 'LineWidth',1, 'color','black')
% ylim([-2,2])
% xlabel('Time (us)') 
% title(sprintf('FM Waveform (mf = %d)', mf))
% saveas(h,'FM-5_Waveform','jpg');
% clf(h)
% 
% yfm_fft = abs(fftshift(fft(yfm)));
% yfm_fft = yfm_fft./max(yfm_fft(centerspec))*abs(besselj(0,mf));
% stem (freqaxis(wilderspec),yfm_fft(wilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% xlabel('Frequency (MHz)') 
% title(sprintf('FM Spectrum (mf = %d)', mf))
% saveas(h,'FM-5_Spectrum','jpg');
% clf(h)
% %%
% yrand = double(randn(T*bps,1) > 0.5);
% ys = repmat(yrand, 1, fs/bps);
% digitalwaveform = reshape(ys',[],1);

% % snr = -20;
% 
% Nbit = 5;
% bpsspec = 1:floor(Nbit*fs/bps);
% 
% 
% %% ASK
% yASK = yc.*(0+1i*digitalwaveform);
% 
% plot(t(bpsspec).*1e6, real(yASK(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% ylim([-2,2])
% xlabel('Time (us)') 
% title(sprintf('ASK Waveform (fc = %dKHz, %dKbps)', round(fc/1e3), round(bps/1e3)))
% saveas(h,'ASK_Waveform','jpg');
% clf(h)
% % semilogy (freqaxis,abs(fftshift(fft(yASK))))
% 
% yASK_fft = psd_estimate(yASK, avg_num, fs);
% semilogy (freqaxis(wilderspec),yASK_fft(wilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% ylim([1e2,1e5])
% xlabel('Frequency (MHz)') 
% title(sprintf('ASK Spectrum (fc = %dKHz, %dKbps)', round(fc/1e3), round(bps/1e3)))
% saveas(h,'ASK_Spectrum','jpg');
% clf(h)
% 
%% FSK
% yFSK = yc1.*(0+1i*digitalwaveform) + yc2.*(0+1i*(~digitalwaveform));
% % semilogy (freqaxis,abs(fftshift(fft(yFSK))))
% % yFSK_fft = sum(abs(fftshift(fft(reshape(yFSK, [], avg_num), length(yFSK)))),2);
% % semilogy (freqaxis(wilderspec),yFSK_fft(wilderspec), 'LineWidth',1)
% plot(t(bpsspec).*1e6, real(yFSK(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% ylim([-2,2])
% xlabel('Time (us)') 
% title(sprintf('FSK Waveform (f1 = %dKHz, f2 = %dKHz, %dKbps)', round(fc1/1e3), round(fc2/1e3), round(bps/1e3)))
% saveas(h,'FSK_Waveform','jpg');
% clf(h)
% 
% save_time_fig(t(bpsspec), yFSK(bpsspec), fc, 'FSK')
% 
% yFSK_fft = psd_estimate(yFSK, avg_num, fs);
% semilogy (freqaxis(wilderspec),yFSK_fft(wilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
% ylim([1e2,1e5])
% xlabel('Frequency (MHz)') 
% title(sprintf('FSK Spectrum (f1 = %dKHz, f2 = %dKHz, %dKbps)', round(fc1/1e3), round(fc2/1e3), round(bps/1e3)))
% saveas(h,'FSK_Spectrum','jpg');
% clf(h)



%% BPSK
% nPSK = 2;
% digitalwaveformangle = randi([0,nPSK-1],T*bps/log2(nPSK),1)./nPSK.*2*pi;
% digitalwaveformangle(1:2) = [0:(nPSK-1)]./nPSK.*2*pi;
% digitalwaveformangle = repmat(digitalwaveformangle, 1, fs/(bps/log2(nPSK)));
% digitalwaveformangle = reshape(digitalwaveformangle',[],1);
% % digitalwaveform = awgn(digitalwaveform,snr);
% yBPSK = yc.*exp(1i*digitalwaveformangle);
% save_time_fig(t(bpsspec), yBPSK(bpsspec), fc, 'BPSK')
% yBPSK_fft = psd_estimate(yBPSK, avg_num, fs);
% save_spectrum_fig(fftaxis(fftwilderspec), yBPSK_fft(fftwilderspec), nPSK_psd(nPSK, bps, fc, fftaxis(fftwilderspec)*1e6), fc, 'BPSK')

%% QPSK
% nPSK = 4;
% digitalwaveformangle = randi([0,nPSK-1],T*bps/log2(nPSK),1)./nPSK.*2*pi;
% digitalwaveformangle(1:4) = [0:3]./nPSK.*2*pi;
% digitalwaveformangle = repmat(digitalwaveformangle, 1, fs/(bps/log2(nPSK)));
% digitalwaveformangle = reshape(digitalwaveformangle',[],1);
% % digitalwaveform = awgn(digitalwaveform,snr);
% yQPSK = yc.*exp(1i*digitalwaveformangle);
% bpsspec = 1:floor(Nbit*fs/bps*log2(nPSK));
% save_time_fig(t(bpsspec), yQPSK(bpsspec), fc, 'QPSK')
% yQPSK_fft = psd_estimate(yQPSK, avg_num, fs);
% save_spectrum_fig(fftaxis(fftwilderspec), yQPSK_fft(fftwilderspec), nPSK_psd(nPSK, bps, fc, fftaxis(fftwilderspec)*1e6), fc, 'QPSK')
%% 8PSK
% nPSK = 8;
% digitalwaveformangle = randi([0,nPSK-1],T*bps/log2(nPSK),1)./nPSK.*2*pi;
% digitalwaveformangle(1:8) = [0:7]./nPSK.*2*pi;
% digitalwaveformangle = repmat(digitalwaveformangle, 1, fs/(bps/log2(nPSK)));
% digitalwaveformangle = reshape(digitalwaveformangle',[],1);
% % digitalwaveform = awgn(digitalwaveform,snr);
% yOPSK = yc.*exp(1i*digitalwaveformangle);
% bpsspec = 1:floor(Nbit*fs/bps*log2(nPSK));
% save_time_fig(t(bpsspec), yOPSK(bpsspec), fc, '8PSK')
% yOPSK_fft = psd_estimate(yOPSK, avg_num, fs);
% save_spectrum_fig(fftaxis(fftwilderspec), yOPSK_fft(fftwilderspec), nPSK_psd(nPSK, bps, fc, fftaxis(fftwilderspec)*1e6), fc, '8PSK')


%%

% set(h,'visible','on');
% nPSK = 2;
% semilogy (fftaxis(ffthalfwilderspec),yBPSK_fft(ffthalfwilderspec),':' ,...
%     'LineWidth',2, 'color','blue')
% hold on 
% semilogy (fftaxis(ffthalfwilderspec),nPSK_psd(nPSK, bps, fc, fftaxis(ffthalfwilderspec)*1e6) ,'-', ...
%     'LineWidth',2, 'color','blue')
% hold on 
% 
% nPSK = 4;
% semilogy (fftaxis(ffthalfwilderspec),yQPSK_fft(ffthalfwilderspec),':' ,...
%     'LineWidth',2, 'color','black')
% hold on 
% semilogy (fftaxis(ffthalfwilderspec),nPSK_psd(nPSK, bps, fc, fftaxis(ffthalfwilderspec)*1e6) ,'-', ...
%     'LineWidth',2, 'color','black')
% hold on 
% 
% 
% nPSK = 8;
% semilogy (fftaxis(ffthalfwilderspec),yOPSK_fft(ffthalfwilderspec),':' ,...
%     'LineWidth',2, 'color','red')
% hold on 
% semilogy (fftaxis(ffthalfwilderspec), nPSK_psd(nPSK, bps, fc, fftaxis(ffthalfwilderspec)*1e6) ,'-' ,...
%     'LineWidth',2, 'color','red')
% hold on 
% 
% 
% ylim([1e-8,1e-3])
% xlim([0.1,0.15])
% 
% legend('BPSK Sim','BPSK Theory','QPSK Sim','QPSK Theory','8PSK Sim','8PSK Theory')
% 
% xlabel('Frequency (MHz)') 
% title(sprintf('MPSK Spectrum (fc = %dKHz, %dKbps)', round(fc/1e3), round(bps/1e3)))
% saveas(h,'MPSK_Spectrum','jpg');
% clf(h)
% hold off
%% QAM
% y = modulate(i,70,fs,'qam',q);
nQAM = 16;
qamset = [(-sqrt(nQAM)/2+1):(sqrt(nQAM)/2)]*2-1;
digitalwaveformidxI = randi([1 sqrt(nQAM)],T*bps/log2(nQAM),1);
digitalwaveformidxQ = randi([1 sqrt(nQAM)],T*bps/log2(nQAM),1);
digitalwaveformamp = (qamset(digitalwaveformidxI) + 1j*qamset(digitalwaveformidxQ))./sqrt(2*(sqrt(nQAM)-1)^2);
digitalwaveformamp = repmat(digitalwaveformamp', 1, fs/(bps/log2(nQAM)));
digitalwaveformamp = reshape(digitalwaveformamp',[],1);
% plot(lags,c)

% digitalwaveform = awgn(digitalwaveform,snr);
yQAM = yc.*digitalwaveformamp;
yQAM_fft = psd_estimate(yQAM, avg_num, fs);
semilogy(fftaxis(fftwilderspec), yQAM_fft(fftwilderspec))
hold on
semilogy(fftaxis(fftwilderspec), nQAM_psd(nQAM, bps, fc, fftaxis(fftwilderspec)*1e6) ,'-' ,...
    'LineWidth',1, 'color','red')
% semilogy(fftaxis, nQAM_psd(nQAM, bps, fc, fftaxis*1e6) ,'-' ,...
%     'LineWidth',1, 'color','red')
ylim([1e-8,1e-3])


% plot(real(yQAM),imag(yQAM),'.')
% GenIQData(yam_am, 'yam_am.txt', fc,fs);
% GenIQData(yam_dsb, 'yam_dsb.txt', fc,fs);
% GenIQData(yam_lsb, 'yam_lsb.txt', fc,fs);
% GenIQData(yam_usb, 'yam_usb.txt', fc,fs);
% GenIQData(yfm_1, 'yfm_1.txt', fc,fs);
% GenIQData(yfm_5, 'yfm_5.txt', fc,fs);
% GenIQData(yASK, 'yASK.txt', fc,fs);
% GenIQData(yFSK, 'yFSK.txt', fc,fs);
% GenIQData(yBPSK, 'yBPSK.txt', fc,fs);
% GenIQData(yQPSK, 'yQPSK.txt', fc,fs);
% GenIQData(yOPSK, 'yOPSK.txt', fc,fs);
