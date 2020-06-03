clear all; close all; clc;

nFFT = 32;
cpLen = 20;
numRB = 50;
bitsPerSubCarrier = 2;
snrdB = 18;
toneOffset = 7;        % Tone offset or excess bandwidth (in subcarriers)
L = 513;                 % Filter length (=filterOrder+1), odd

numDataCarriers = numRB*nFFT;    % number of data subcarriers in subband
halfFilt = floor(L/2);
n = -halfFilt:halfFilt;

% Sinc function prototype filter
pb = sinc((numDataCarriers+2*toneOffset).*n./nFFT);

% Sinc truncation window
w = (0.5*(1+cos(2*pi.*n/(L-1)))).^0.6;

% Normalized lowpass filter coefficients
fnum = (pb.*w)/sum(pb.*w);

bitsIn = randi([0 1],numRB*nFFT*bitsPerSubCarrier,1);
symbolsIn = my_QPSK(bitsIn);

% apply IFFT
ifftOut = [];
for i=0:numRB-1
    blockIFFT = ifft(symbolsIn((i*nFFT+1):(i+1)*nFFT));
    CP = blockIFFT(1:cpLen);
    ifftOut = [ifftOut blockIFFT CP]; %CP added
end
ifftOut = filter2(fnum,ifftOut); % Baseband filtering applied

% Plot power spectral density (PSD)
[psd,f] = periodogram(ifftOut);
plot(f,10*log10(psd),'k')
title('Power Spectrum of Filtered OFDM transmitted Signal')
xlabel('Normalized Frequency'),ylabel('PSD (dBW/Hz)')
grid on

% signal transmitted through channel

channel = randn(1,10);
receivedSig = conv(channel,ifftOut);
receivedSig = awgn(receivedSig, snrdB,'measured');
receivedSig = receivedSig(length(channel):end);

receivedSig = filter2(fnum,receivedSig);

%removing CP and performing FFT
rxSig = [];
for i=0:numRB-1
    temp = receivedSig((i*(nFFT+cpLen)+1):((i+1)*(nFFT+cpLen)));
    cpRemovedSig = temp(1:(end-cpLen));
    rxSig = [rxSig fft(cpRemovedSig)];
end

%Symbol detection
G = fft([channel(2:length(channel)) zeros(1,nFFT-length(channel)) channel(1)]);
detectedSig = [];
for i=0:numRB-1
    temp = rxSig((i*nFFT+1):((i+1)*nFFT));
    detectedSig = [detectedSig temp./G];
end
scatterplot(detectedSig,4)
grid on
scatterplot(symbolsIn,4)
grid on


