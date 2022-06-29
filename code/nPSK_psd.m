function psd = nPSK_psd(nPSK, bps, fc, fftaxis)
pskset = exp(1i*([1:nPSK]-1)./nPSK.*2*pi);
Tb = 1/bps;
T = log2(nPSK)*Tb;
EavgI = sum(real(pskset).^2)/nPSK;
AavgI = sum(real(pskset))/nPSK;
EavgQ = sum(imag(pskset).^2)/nPSK;
AavgQ = sum(imag(pskset))/nPSK;
fnorm = fftaxis-fc;
psd = (EavgI-AavgI^2)*T*(sinc(fnorm*T).^2)+(EavgQ-AavgQ^2)*T*(sinc(fnorm*T).^2);