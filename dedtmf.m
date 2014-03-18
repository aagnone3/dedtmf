function [Y,ER,Fs,Rs,Ts,Bfe,Afe] = dedtmf(X,P,W,H,params)
% [Y,ER,Fs,Rs,Ts,Bf,Af] = dedtmf(X,P,W,H,params)
%   Remove DTMF (and other stationary tones) from signal
%   Inputs:
%    X  - input waveform
%    P  - order of LPC models to fit 
%    W  - length of blocks for LPC fitting
%    H  - hop between successive blocks
%    params - additional parameters:
%     params.poleradthreash (0.98) - radius threshold for mapping poles to zeros
%     params.poleradtrans (0.002) - radius transition width for mapping
%     params.polerad (0.98) - radius for added compensatory poles
%     params.noiseexcite (0) - use white noise excitation (to debug)
%   Returns:
%    Y  - audio filtered to remove steady tones
%    ER - block-wise energy ratio of original to filtered (close to
%         zero where significant tones have been removed.
%    Fs, Rs, Ts - frequencies (0..pi), radii (0..1) and times (in
%         samples) of the main components found.
%    Bf, Af - per-frame coefficients of the tone-removal filters used.
%
% 2013-11-22 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; P = 16; end
if nargin < 3; W = 1024; end
if nargin < 4; H = W/2; end
if nargin < 5; params.dummy = []; end  % to make it a struct

VERSION=0.1;
DATE=20131126;

% default parameters
if isfield(params, 'poleradthresh') == 0; params.poleradthresh = 0.98; end
if isfield(params, 'poleradtrans') == 0; params.poleradtrans = 0.002; end
if isfield(params, 'polerad') == 0; params.polerad = 0.98; end
if isfield(params, 'noiseexcite') == 0; params.noiseexcite = 0; end

% Break x into w-length blocks stepped by h
Xf = frame(X,W,H);

% Fit LPC models to every frame
Af = lpc(Xf, P);

% Find all the poles of the per-frame LPC filters
Rf = rootsbyrow(Af);
% Magnitudes of each pole
Mf = abs(Rf);
% Corresponding pure-phase parts
Cf = Rf./Mf;

% Build times, magnitudes, and angles of all positive freq roots 
% to return in diagnostics
Ts = repmat(W/2+H*[0:size(Mf,1)-1]',1,size(Mf,2));
Ts = Ts(:);
Rs = Mf(:);
Fs = angle(Cf(:));
% Keep just +ve freq complex poles
keep = find(Fs>0 & Fs <pi);
Ts = Ts(keep);
Rs = Rs(keep);
Fs = Fs(keep);

% Modify magnitude; ensure smaller than 1
% K = 0.25;
%Mfm = (1-abs(1-Mf)).^K;
% Use a sigmoid to make poles close to 1 even closer
%Mfm = (Mf > poleradthresh);
Mfm = sigmoid( (Mf - params.poleradthresh)/params.poleradtrans );
% New polynomial with exaggerated pole radii
Bfe = polybyrow( Mfm .* Cf );
% Compensatory poles - just inside zeros
Afe = polybyrow( params.polerad*Mfm .* Cf);

% Apply inverse filters to each block of X
if params.noiseexcite
  % Use white noise excitation of output to debug filtering
  E = randn(length(X), 1);
  Xf = frame(E,W,H);
end
[nr,nc] = size(Xf);
Xff = zeros(nr, nc);
for i = 1:nc
  Xff(:,i) = filter(Bfe(i,:),Afe(i,:),Xf(:,i)')';
end

% Reconstruct with 50%-overlapped hanning window
L = 2 * H;
pad = W - L;
win = [zeros(1,floor(pad/2)), hanning(L)', zeros(1,ceil(pad/2))];

% Figure per-block energy ratio
ER = sum(bsxfun(@times, Xff, win').^2)./sum(bsxfun(@times, Xf, win').^2);
% Threshold
ER = min(ER, 1.0);
% scale
%Xff = repmat(ER, size(Xff,1), 1).* Xff;

%Y = ola(diag(win)*Xff, H);
Y = ola(bsxfun(@times, Xff, win'), H);
%Y = ola(diag(win)*Xf, H);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = rootsbyrow(A)
% Calculate the roots for the polynomials on each row of A
[nr,nc] = size(A);
R = zeros(nr, nc-1);
for i = 1:nr
  R(i,:) = roots(A(i,:))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = polybyrow(R)
% Calculate the polynomial coefficients for rows of roots in R
[nr,nc] = size(R);
A = zeros(nr, nc+1);
for i = 1:nr
  A(i,:) = poly(R(i,:)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = sigmoid(X)
Y = 1./(1+exp(-X));
