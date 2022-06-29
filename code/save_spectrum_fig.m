function save_spectrum_fig(psd_spec, psd, psd_theory, fc, titlein)

h=figure(1);
clf(h)
set(h,'visible','off');
semilogy (psd_spec,psd, 'LineWidth',1, 'color',[239/255 143/255 38/255])
hold on 
semilogy (psd_spec,psd_theory, 'LineWidth',1, 'color','black')
xlabel('Frequency (MHz)') 
legend(sprintf('%s Sim', titlein),sprintf('%s Theory', titlein))

ylim([1e-8,1e-3])
xlim([0.05,0.15])

title(sprintf('%s Spectrum (fc = %dKHz)', titlein, round(fc/1e3)))
saveas(h,sprintf('%s_Spectrum', titlein),'jpg');
clf(h)