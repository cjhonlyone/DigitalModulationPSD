function save_time_fig(time_spec, time, fc, titlein)

h=figure(1);
clf(h)
set(h,'visible','off');

subplot(211)
plot (time_spec.*1e6, real(time), 'LineWidth',1, 'color',[239/255 143/255 38/255])
ylim([-2,2])
xlabel('Time (us)') 
title(sprintf('%s Waveform I (fc = %dKHz)', titlein, round(fc/1e3)))
subplot(212)
plot (time_spec.*1e6, imag(time), 'LineWidth',1, 'color',[239/255 143/255 38/255])
ylim([-2,2])
xlabel('Time (us)') 
title(sprintf('%s Waveform Q (fc = %dKHz)', titlein, round(fc/1e3)))

saveas(h,sprintf('%s_Waveform', titlein),'jpg');
clf(h)