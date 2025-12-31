%% 函数: GCC-PHAT TDOA
function [tauHat, Rtau, lags] = gcc_phat_tdoa(Z1, Z2, fs)
%STFT-GCC_PHAT_TDOA 使用 STFT 复矩阵计算时延差 (已修正倍数和符号方向)
% Z1, Z2 : 频率 × 时间 复矩阵
% fs : 采样率 (Hz)
% 返回:
% tauHat : 估计时延 (s)，取互相关最大峰值处
% Rtau : 双边互相关序列
% lags : 与 Rtau 对应的滞后轴 (s)
% 1) 交叉功率谱并 PHAT 归一化
S12 = Z1 .* conj(Z2);
S12 = S12 ./ (abs(S12) + eps);
% 2) 沿时间帧平均
G = mean(S12, 2); % [freqBins × 1]
% 3) 逆 FFT → 互相关 (注意：直接使用 STFT 输出长度，不扩展)
Nfft = size(Z1,1);
R = ifft(G, Nfft);

Rtau = ifftshift(R);

half = floor(Nfft/2);
% lags = 0:length(R)-1;


if mod(Nfft,2) ==0
   lags = (-half:half-1).' / fs;
else
   lags = (0:Nfft-1)/fs - (Nfft+1)/fs/2;
end
   
   % lags = (1:length(R))/fs - (length(R))/fs/2;


% 4) 生成双边互相关序列
% if mod(length(R),2) == 0
% 
%     half = floor(Nfft/2);
%     Rtau = [R(half+1:end); R(1:half)];
%     lags = (-half:half-1).' / fs;
% else
%     half = floor(Nfft/2);
%     Rtau1 = [R(half+1:end); R(1:half)];
%     lags = (-half:half).' / fs;
% 
%     Rtau = ifftshift(R);
%     aa = 1;
% end
% 5) 找最大峰值
[~, idx] = max(abs(Rtau));
tauHat = lags(idx);
end