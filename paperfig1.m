% paperfig1.m
% 
% Matlab to generate the panes for figure 1 in the Interspeech14
% paper.
%
% 2014-03-15 Dan Ellis dpwe@ee.columbia.edu

% Load audio
sr = 8000;
fn = 'irdial/tcp_d1_02_counting_cia_irdial.mp3';
[d,sr] = audioread(fn,sr);

nr = 3;
nc = 4;

% Specgram
subplot(nr, nc, 1)
specgram(d,512,sr)
caxis([-40 40])
% Center time for analysis
tc = 56.1;
% Width of time context to show
tw = 3.0;
% frequency range to show
fmin = 0; fmax = 3000;
axis([tc-tw/2 tc+tw/2 fmin fmax])
[p,n,e] = fileparts(fn);
title(n, 'interpreter', 'none');

% Add line showing focus time
hold on; plot([tc;tc], [fmin fmax], '-r'); hold off

% Calculate lpc for that frame
W = 4096; H = 128; P = 20;
Xf = frame(X,W,H);
Af = lpc(Xf, P);

% Center times of each frame
TT = (W/2+H*[0:size(Af,1)-1])/sr;
% Find one closest to tc
[vv,frm] = min(abs(TT-tc));
% The actual LPC fit to this frame
Aff = Af(frm, :);

% Second pane is spec of segment and LPC fit at tc
subplot(nr, nc, 2)
fftlen = size(Xf,1);
XX = fft(hanning(1,fftlen).*Xf(:,frm));
XX = XX(1:(fftlen/2+1));
ff = [0:(length(XX)-1)]*sr/fftlen;
% LPC spectrum
[HH,WW] = freqz(1,Aff);
% plot
plot(ff, 20*log10(abs(XX)), WW*sr/(2*pi), 20+20*log10(abs(HH)),'-r');
f2min = 500; f2max = 1500;
axis([f2min f2max -20 60])
grid
title('Spectrum and LPC fit');

% Third pane is poles
subplot(nr, nc, 3)
zplane(1, Aff);
zpax = [-1.1 1.1 -1.1 1.1];
axis(zpax);
zpzoom = [0.68 0.72 0.68 .72];
hold on; plot(zpzoom([1 2 2 1 1]), zpzoom([3 3 4 4 3]), '-r'); hold off
title('LPC poles')

% fourth pane zooms in poles
subplot(nr, nc, 4)
zplane(1, Aff);
hold on; plot(zpzoom([1 2 2 1 1]), zpzoom([3 3 4 4 3]), '-r'); hold off
axis(zpzoom);
grid
title('LPC poles detail')

% Now process to build notch filter
Rf = roots(Aff);
% Magnitudes of each pole
Mf = abs(Rf);
% Corresponding pure-phase parts
Cf = Rf./Mf;

% Use a sigmoid to make poles close to 1 even closer
poleradthresh = 0.97; % 0.98;
poleradtrans =  0.002;% 0.001;
% (supply local fn def for sigmoid)
sigmoid = @(X) 1./(1+exp(-X));
Mfm = sigmoid( (Mf - poleradthresh)/poleradtrans );
% New polynomial with exaggerated pole radii
Bfe = poly( Mfm' .* Cf' );
% Compensatory poles - just inside zeros
polerad = 0.8;  % 0.98;
Afe = poly( polerad*Mfm' .* Cf');

% Plot the modified pole-zero diagram
subplot(nr, nc, nc+3)
zplane(Bfe, Afe)
axis(zpax);
hold on; plot(zpzoom([1 2 2 1 1]), zpzoom([3 3 4 4 3]), '-r'); hold off
title('Transformed filter poles/zeros')

% fourth pane zooms in poles
subplot(nr, nc, nc+4)
zplane(Bfe, Afe);
hold on; plot(zpzoom([1 2 2 1 1]), zpzoom([3 3 4 4 3]), '-r'); hold off
axis(zpzoom);
grid
title('Transformed filter detail')

% Plot the resulting frequency response
XXf = fft(hanning(1,fftlen).*filter(Bfe, Afe, Xf(:,frm)));
XXf = XXf(1:(fftlen/2+1));
% LPC spectrum
[HHf,WWf] = freqz(Bfe,Afe);
% plot
subplot(nr, nc, nc+2)
plot(ff, 20*log10(abs(XXf)), WWf*sr/(2*pi), 40+20*log10(abs(HHf)),'-r');
axis([f2min f2max -20 60])
grid
title('Filter response and filtered spectrum')

% Spectrogram of result of running the process
[y,E,F,R,T] = dedtmf(d, P, W, H);
subplot(nr, nc, nc+1)
specgram(y,512,sr)
caxis([-40 40])
axis([tc-tw/2 tc+tw/2 fmin fmax])
title('Filtered signal')

colormap(1-gray);

% Add line showing focus time
hold on; plot([tc;tc], [fmin fmax], '-r'); hold off





