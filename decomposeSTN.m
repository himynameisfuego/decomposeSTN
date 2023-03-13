% -------------------------------------------------------------------------
% [xs,xt,xn] = decomposeSTN(x,Fs,nWin)
% 
% Decomposition of time-domain input into tonal, transient and noise
% components.
%
% Inputs:
% x = input audio signal;
% Fs = sample rate
% nWin = [nWin1 nWin2] window lengths for each decomposition round
%
% Outputs:
% xs = tonal component; xt = transients; xn = noise.
%
% L. Fierro, Acoustics Lab, Aalto University, 2020
% Last modified: 15.04.2021 - L. Fierro 
% -------------------------------------------------------------------------

function [xs,xt,xn] = decomposeSTN(x,Fs,nWin)

if nargin < 3
    nWin = [8192 512];
end

nHop = nWin / 8; NFFT = nWin; 
filter_length_t = 200e-3; % in ms
filter_length_f = 500; % in Hz

% Round 1
win = hann(nWin(1),'periodic');
[X,~] = stft(x,win,nHop(1),NFFT(1));
nMedianH = round(filter_length_t * Fs / nHop(1));
nMedianV = round(filter_length_f * NFFT(1) / Fs);
Rt = transientness(X,nMedianH,nMedianV);

[S,T,N] = decSTN(Rt,0.7,0.8);

xs = istft(S.*X, nHop(1), win, win);
xres = istft((T+N).*X, nHop(1), win, win);

% Round 2
win = hann(nWin(2),'periodic');

[X,~] = stft(xres,win,nHop(2),NFFT(2));
nMedianH = round(filter_length_t * Fs / nHop(2) );
nMedianV = round(filter_length_f * NFFT(2) / Fs);
Rt = transientness(X,nMedianH,nMedianV);

[S,T,N] = decSTN(Rt,0.75,0.85);

xt = istft(T.*X, nHop(2), win, win);
xn = istft((S+N).*X, nHop(2), win, win);

end


function [S,T,N] = decSTN(Rt, G2, G1)

if nargin < 2
    G1 = 0.9;    % Upper threshold    
    G2 = 0.75;    % Lower threshold
end

Rs = 1-Rt;  

S = sin(pi*(Rs-G2)/(2*(G1-G2))).^2;    
S(Rs>=G1) = 1; S(Rs<G2) = 0;

T = sin(pi*(Rt-G2)/(2*(G1-G2))).^2; 
T(Rt>=G1) = 1; T(Rt<G2) = 0;

N = zeros(size(S)); N = 1-S-T;
end

function Y = transientness(X,nMedianH,nMedianV)

X_v_median = medfilt1(abs(X),nMedianV,[],1);
X_h_median = medfilt1(abs(X),nMedianH,[],2);

Y = X_v_median ./ (X_v_median + X_h_median);
Y(isnan(Y)) = 0;

end

function [Y,T] = stft(x,win,nHop,NFFT)

nWin = length(win);
L = length(x);

nFrames = floor((L-nWin)/nHop+1);
nBins = NFFT/2+1;
Y = zeros(nBins,nFrames);
T = zeros(1,nFrames);

pin = 0;
x = [x;zeros(nWin/2,1)];
for n = 1:nFrames
    grain = x(pin+1:pin+nWin).*win;
    f = fft(grain,NFFT);
    Y(:,n) = f(1:nBins);
    T(n) = pin + 1;
    pin = pin + nHop;
end
end

function y = istft(X,nHop,win,win_analysis)

nWin = length(win);
[nBins,nFrames] = size(X);
NFFT = (nBins-1)*2;

% Window function
if nargin < 3
    % Default synthesis window
    win = rectwin(NFFT);
end
nWin = length(win);
if nargin < 4
    % Default analysis window
    win_analysis = rectwin(nWin);
end

% Length of output
L = (nFrames-1)*nHop + nWin;
y = zeros(L,1);

% OLA normalization
norm_coef = ola_norm_coef(win_analysis,win,nHop);

% Compute two-sided spectrogram
XF = zeros(NFFT,nFrames);
XF(1:nBins,:) = X(1:nBins,:);
XF(nBins+1:end,:) = conj(flipud(X(2:end-1,:)));

% Overlap-add synthesis
p = 0;
for n = 1:nFrames
    grain = real(ifft(XF(:,n)));
    grain = grain(1:nWin).*win ./ norm_coef;
    y(p+1:p+nWin) = y(p+1:p+nWin) + grain;
    p = p + nHop;
end
end

function y = ola_norm_coef(win_analysis,win_synthesis,nHop_synthesis)

nWin = length(win_analysis);
nHop = nHop_synthesis;

win = win_analysis .* win_synthesis;
idx = nWin / 2 + 1;
y = win(idx);

m = 1;
i = idx - m * nHop;
j = idx + m * nHop;
while i > 0 &&  j <= nWin
    y = y + win(i) + win(j);
    m = m + 1;
    i = idx - m * nHop;
    j = idx + m * nHop;
end
end