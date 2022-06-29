function am_spectrum(analogwaveform, t, dF, fc, fs, bpsspec)



%% AM
ma=2;
yam_am = amiqmod(analogwaveform,fc,fs,0,ma);
yam_am_fft = psd_estimate(yam_am, avg_num, fs);

yam_dsb = amiqmod(analogwaveform,fc,fs);
yam_dsb_fft = abs(fftshift(fft(yam_dsb)));
yam_lsb = ssbiqmod(analogwaveform,fc,fs);
yam_lsb_fft = abs(fftshift(fft(yam_lsb)));
yam_usb = ssbiqmod(analogwaveform,fc,fs,0,'upper');
yam_usb_fft = abs(fftshift(fft(yam_usb)));

h=figure(1);
clf(h)
set(h,'visible','off');
%%
holdposition = h.OuterPosition;
h.OuterPosition(3) = 2*h.OuterPosition(3);
plot(t(bpsspec).*1e6, real(yam_am(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
hold on
plot(t(bpsspec).*1e6, 1.5*analogwaveform(bpsspec), 'LineWidth',1, 'color','black')
ylim([-2,2])
xlabel('Time (us)') 
title(sprintf('AM Waveform (fc = %dKHz, %dKbps)', round(fc/1e3), round(dF/1e3)))
saveas(h,'AM_Waveform_wider','jpg');
clf(h)
h.OuterPosition = holdposition;

plot(t(bpsspec).*1e6, real(yam_am(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
hold on
plot(t(bpsspec).*1e6, 1.5*analogwaveform(bpsspec), 'LineWidth',1, 'color','black')
ylim([-2,2])
xlabel('Time (us)') 
title(sprintf('AM Waveform (fc = %dKHz, %dKbps)', round(fc/1e3), round(dF/1e3)))
saveas(h,'AM_Waveform','jpg');
clf(h)

stem (freqaxis(wilderspec),yam_am_fft(wilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
xlabel('Frequency (MHz)') 
title(sprintf('AM Spectrum (fc = %dKHz, %dKbps)', round(fc/1e3), round(dF/1e3)))
saveas(h,'AM_Spectrum','jpg');
clf(h)
%%
plot(t(bpsspec).*1e6, real(yam_dsb(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
hold on
plot(t(bpsspec).*1e6, 1.5*analogwaveform(bpsspec), 'LineWidth',1, 'color','black')
ylim([-2,2])
xlabel('Time (us)') 
title(sprintf('AM-DSB Waveform (fc = %dKHz, %dKbps)', round(fc/1e3), round(dF/1e3)))
saveas(h,'AM-DSB_Waveform','jpg');
clf(h)
stem (freqaxis(wilderspec),yam_dsb_fft(wilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
xlabel('Frequency (MHz)') 
title(sprintf('AM-DSB Spectrum (fc = %dKHz, %dKbps)', round(fc/1e3), round(dF/1e3)))
saveas(h,'AM-DSB_Spectrum','jpg');
clf(h)
%%
plot(t(bpsspec).*1e6, real(yam_lsb(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
ylim([-2,2])
xlabel('Time (us)') 
title(sprintf('AM-LSB Waveform (fc = %dKHz, %dKbps)', round(fc/1e3), round(dF/1e3)))
saveas(h,'AM-LSB_Waveform','jpg');
clf(h)
stem (freqaxis(wilderspec),yam_lsb_fft(wilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
xlabel('Frequency (MHz)') 
title(sprintf('AM-LSB Spectrum (fc = %dKHz, %dKbps)', round(fc/1e3), round(dF/1e3)))
saveas(h,'AM-LSB_Spectrum','jpg');
clf(h)
%%
plot(t(bpsspec).*1e6, real(yam_usb(bpsspec)), 'LineWidth',1, 'color',[239/255 143/255 38/255])
ylim([-2,2])
xlabel('Time (us)') 
title(sprintf('AM-USB Waveform (fc = %dKHz, %dKbps)', round(fc/1e3), round(dF/1e3)))
saveas(h,'AM-USB_Waveform','jpg');
clf(h)
stem (freqaxis(wilderspec),yam_usb_fft(wilderspec), 'LineWidth',1, 'color',[239/255 143/255 38/255])
xlabel('Frequency (MHz)') 
title(sprintf('AM-USB Spectrum (fc = %dKHz, %dKbps)', round(fc/1e3), round(dF/1e3)))
saveas(h,'AM-USB_Spectrum','jpg');
clf(h)