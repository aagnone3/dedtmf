%% DEDTMF - Tool to remove stationary DTMF-like tones from audio
%
% |dedtmf| is a Matlab script that attempts to suppress stationary tones 
% an audio file, while leaving more dynamic components unchanged.  This 
% could be useful, for instance, to suppress DTMF tones mixed in
% with speech.

%% Principle
%
% Stationary tones (i.e., locally periodic signals) will lead to
% LPC fits with very sharp poles (i.e., z-domain magnitudes close
% to the unit circle).  The routine fits LPC models to fairly long
% (e.g., 250ms) windows of the signal.  A signal whose periodicity
% is stable over that window is detected by its large-radius poles;
% these poles are then converted to zeros on the unit circle along
% with compensatory poles just inside them, leading to a set of
% narrow notches, which are used to filter the tone out of that
% segment.  A small segment from the middle of the block is used in
% overlap-add to reconstruct the original signal with the steady
% tones attenuated by at least 40 dB (by eye).

%% Example usage
%
% The example waveform includes dial tone, the DTMF tones during
% dialing, then some speech mixed with a repeating tone.

[d,sr] = wavread('phonexample.wav');
soundsc(d,sr);

% Now we suppress it using a 40-pole fit, over 4096 point windows
% every 256 samples
H = 256;
W = 4096;
[y,E,F,R,T] = dedtmf(d, 40, W, H);

% Listening to result, tones are largely suppressed
soundsc(y,sr);

% Compare spectrograms
tlim = 3.0;
axs(1) = subplot(311);
specgram(d,256,sr);
caxis([-50 30]);
axis([0 tlim 0 4000]);
title('Original');
axs(2) = subplot(312);
specgram(y,256,sr);
caxis([-50 30]);
axis([0 tlim 0 4000]);
title('Tones suppressed')

% Plot energy ratio
axs(3) = subplot(313);
plot(([0:length(E)-1]*H+W/2)/sr, E);
%axis([0 length(d)/sr 0 1.1]);
axis([0 tlim 0 1.1]);
title('Frame energy ratio - filtered / original')

linkaxes(axs, 'x');

% Show where the main detected tones were
mainpoles = find(R>.99);
subplot(311)
hold on;
plot(T(mainpoles)/sr, F(mainpoles)/pi*sr/2, '.w');
hold off;

% Notice how the DTMF tones are completely eliminated, leading to
% energy ratios close to zero.  The speech, however, is largely
% untouched.  Some of the other tones have complex AM that makes
% some of their components not be detected as steady tones, and
% hence not removed.

%% Limitations
%
% The tones are removed by stationary notches, albeit applied on
% short blocks (256 samples = 16 ms in this eample).  However, if
% the tones are not entirely stable, they will not be completely
% eliminated.  You can hear that with the ascending tones in this
% example, which have a bit of tape flutter and thus are not
% exactly fixed frequency.

%% Installation
% 
% The Matlab code is available at
% 
%   <http://labrosa.ee.columbia.edu/projects/dedtmf/>
%
% All sources are in the package <dedtmf.zip>.
%
% Feel free to contact me with any problems.
%

%% Python Port
%
% I have ported this code to Python.  See
% <https://github.com/dpwe/dedtmf dedtmf.py>.

%% Changelog

% 2013-11-26 v0.1 Initial release

%% Acknowledgment
%
% This work was supported by DARPA under the RATS program via a 
% subcontract from the SRI-led team SCENIC.  My work was on behalf 
% of ICSI.
%
% $Header: /Users/dpwe/docs/grants/2010-01-DARPA-RATS/code/dedtmf/RCS/demo_dedtmf.m,v 1.2 2011/05/03 20:32:37 dpwe Exp dpwe $
