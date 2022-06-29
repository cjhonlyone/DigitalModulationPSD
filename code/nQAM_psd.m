function psd = nQAM_psd(nQAM, bps, fc, fftaxis)
qamset = [(-sqrt(nQAM)/2+1):(sqrt(nQAM)/2)]*2-1;
qamset = repmat(qamset, sqrt(nQAM),1) + 1i*repmat(qamset', 1,sqrt(nQAM));
qamset = reshape(qamset, 1, [])./sqrt(2*(sqrt(nQAM)-1)^2);
Tb = 1/bps;
T = log2(nQAM)*Tb;
Enorm = abs(sum(sum(abs(qamset + 1j*qamset').^2)))./nQAM./max(max(abs(qamset + 1j*qamset').^2));
fnorm = fftaxis-fc;
% psd = Enorm*T*(sinc(fnorm*T).^2);


EavgI = sum(real(qamset).^2)/nQAM;
AavgI = sum(real(qamset))/nQAM;
EavgQ = sum(imag(qamset).^2)/nQAM;
AavgQ = sum(imag(qamset))/nQAM;
fnorm = fftaxis-fc;
% psd = EavgI*T*(sinc(fnorm*T).^2)+EavgQ*T*(sinc(fnorm*T).^2);

psd = (EavgI-AavgI^2)*T*(sinc(fnorm*T).^2)+(EavgQ-AavgQ^2)*T*(sinc(fnorm*T).^2);