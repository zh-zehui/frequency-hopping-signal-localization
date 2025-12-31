clear; close all; clc;

%% 1. 参数设置模块 - 基带等效模型
fprintf('设置系统参数（基带等效模型）...\n');

% 跳频信号参数
T_hop = 10e-3;          % 每跳持续时间 10ms
SymbolRate = 50e3;      % 符号速率 50k Baud
Modulation = 'QPSK';
NumHops = 4;            % 用于定位的跳数 (40ms / 10ms)
NumChannels = 20;       % 频点数量

% 跳频频点定义（基带等效频率，相对于中心频率的偏移）
IF_CenterFreqs_MHz = [48.16, 48.48, 48.80, 49.12, 49.44, 49.76, 50.08, 50.40, ...
                      50.72, 51.04, 51.36, 51.68, 52.00, 52.32, 52.64, 52.96, ...
                      53.28, 53.60, 53.92, 54.24];
                  
% 计算基带等效频率（相对于中心频率 51.20 MHz 的偏移）
IF_Center = 51.20e6;    % 中频中心频率
Baseband_Freqs = (IF_CenterFreqs_MHz * 1e6) - IF_Center; % 基带频率偏移

fprintf('基带等效频率范围: %.2f MHz 到 %.2f MHz\n', ...
    min(Baseband_Freqs)/1e6, max(Baseband_Freqs)/1e6);

% 接收机参数
Fs = 6.4e6;             % 基带采样率 6.4MHz（满足信号带宽需求）
NumStations = 5;        % 观测站数量

% 定位场景参数
% === 修正：非对称接收站布局 ===
TargetPos = [100, -100];
StationPos = [
    15e3, 2e3;      % 东北方向，距离约15.1km
    -12e3, 8e3;     % 西北方向，距离约14.4km  
    5e3, -18e3;     % 东南方向，距离约18.7km
    -8e3, -6e3      % 西南方向，距离约10.0km
    10e3, -15e3
];

% 计算理论时延（以第1个站为参考站）
TrueTDOA = zeros(NumStations, 1);
c = 3e8;
distances = vecnorm(StationPos - TargetPos, 2, 2);
fprintf('各站到目标的距离:\n');
for l = 1:NumStations
    TrueTDOA(l) = (distances(l) - distances(1)) / c;
    fprintf('  站%d: %.1f km, 相对时延: %.1f μs\n', ...
        l, distances(l)/1000, TrueTDOA(l)*1e6);
end


% 算法参数
SNR_dB = -15;
% 设置固定的噪声功率
FixedNoisePower = 1;  % 固定噪声功率为1（可根据需要调整）
fprintf('固定噪声功率: %.2f\n', FixedNoisePower);

fprintf('参数设置完成.\n');

%% 2. 信号生成模块 - 基带等效信号
fprintf('生成基带等效跳频QPSK信号...\n');

% 每跳符号数
SamplesPerSymbol = floor(Fs / SymbolRate);
SymbolsPerHop = floor(T_hop * SymbolRate);

% 总采样点数
TotalSamples = NumHops * T_hop * Fs;

% 生成跳频图案
HopPattern = mod(0:NumHops-1, NumChannels) + 1;

% 生成基带等效信号
TimeVector = (0:TotalSamples-1) / Fs;
TxSignal_Baseband = zeros(1, TotalSamples);
SamplesPerHop = round(T_hop * Fs);
for h = 1:NumHops
    % 生成随机QPSK符号
    DataSymbols = (randi(4, 1, SymbolsPerHop) - 1);
    QPSKSymbols = exp(1j * (pi/4 + DataSymbols * pi/2));
    
    % 上采样成型（矩形脉冲）
    SymbolWaveform = reshape(repmat(QPSKSymbols, SamplesPerSymbol, 1), 1, []);
    
    

    % 确定当前跳的时间范围
    HopStartSample = (h-1) * SamplesPerHop + 1;
    HopEndSample = h * SamplesPerHop;
    HopSamples = HopStartSample:HopEndSample;
    
    % 调整长度匹配
    CurrentLength = length(HopSamples);
    if length(SymbolWaveform) > CurrentLength
        SymbolWaveform = SymbolWaveform(1:CurrentLength);
    elseif length(SymbolWaveform) < CurrentLength
        SymbolWaveform = [SymbolWaveform, zeros(1, CurrentLength - length(SymbolWaveform))];
    end
    
    % 调制到当前跳频点（基带等效频率）
    CurrentFreq = Baseband_Freqs(HopPattern(h));  % 使用基带频率偏移
    Carrier = exp(1j * 2 * pi * CurrentFreq * TimeVector(HopSamples));
    TxSignal_Baseband(HopSamples) = SymbolWaveform .* Carrier;
end

TxSignal = TxSignal_Baseband;
fprintf('基带等效信号生成完成. 跳频图案: %s\n', mat2str(HopPattern));

%% 3. 信道模型模块 - 基带等效
fprintf('模拟信号传播（基带等效）...\n');

RxSignals = zeros(NumStations, TotalSamples);

for l = 1:NumStations
    % 计算时延对应的采样点数
    DelaySamples = round(TrueTDOA(l) * Fs);
    
    % 应用时延（基带等效模型中的时延表现为相位旋转）
    % 对于窄带信号，时延 ≈ 相位旋转
    DelayedSignal = circshift(TxSignal, DelaySamples);
    

    % === 修正：保持噪声功率不变，调节信号功率 ===
    % 计算当前信号功率
    CurrentSignalPower = mean(abs(DelayedSignal).^2);
    % 根据SNR计算需要的信号功率
    RequiredSignalPower = FixedNoisePower * (10^(SNR_dB/10));
    % 调整信号幅度以达到目标SNR
    ScaleFactor = sqrt(RequiredSignalPower / CurrentSignalPower);
    ScaledSignal = DelayedSignal * ScaleFactor;

    % 添加固定功率的复高斯白噪声
    Noise = sqrt(FixedNoisePower/2) * (randn(1, TotalSamples) + 1j * randn(1, TotalSamples));
    
    RxSignals(l, :) = ScaledSignal + Noise;

%     % 添加复高斯白噪声
%     SignalPower = mean(abs(DelayedSignal).^2);
%     NoisePower = SignalPower / (10^(SNR_dB/10));
%     Noise = sqrt(NoisePower/2) * (randn(1, TotalSamples) + 1j * randn(1, TotalSamples));
%     
%     RxSignals(l, :) = DelayedSignal + Noise;
    % 验证实际的SNR
    ActualSignalPower = mean(abs(ScaledSignal).^2);
    ActualNoisePower = mean(abs(Noise).^2);
    ActualSNR = 10 * log10(ActualSignalPower / ActualNoisePower);
    fprintf('  站%d: 目标SNR=%.1fdB, 实际SNR=%.1fdB\n', l, SNR_dB, ActualSNR);

end

fprintf('信道模拟完成.\n');

%% 4. 算法一：原始直接定位算法（基带等效版本）
fprintf('执行原始直接定位算法...\n');

% 搜索网格参数
SearchRange = 1000;
GridSize = 21;      % 减小网格大小以提高速度
x_grid = linspace(-SearchRange, SearchRange, GridSize);
y_grid = linspace(-SearchRange, SearchRange, GridSize);


% === 关键修改：只使用当前跳的4个活跃频点 ===
ActiveHopPattern = HopPattern;  % [1, 2, 3, 4] 对应实际的跳频图案
ActiveFreqs = Baseband_Freqs(ActiveHopPattern);  % 只取这4个频点

fprintf('使用的活跃频点: %s\n', mat2str(ActiveHopPattern));
fprintf('对应的基带频率: ');
for i = 1:length(ActiveFreqs)
    fprintf('%.3fMHz ', ActiveFreqs(i)/1e6);
end
fprintf('\n');

% 更稳健的方法：使用MATLAB的频率轴查找
freq_axis = (0:TotalSamples-1) * Fs / TotalSamples;
% 将频率轴调整为从 -Fs/2 到 +Fs/2
freq_axis_shifted = freq_axis;
freq_axis_shifted(freq_axis > Fs/2) = freq_axis(freq_axis > Fs/2) - Fs;

% % === 修正：正确的频率索引计算（处理负频率）===
% k_m  = zeros(1, length(Baseband_Freqs));
% 
% % 找到每个基带频率对应的最近DFT点
% for m = 1:length(Baseband_Freqs)
%     [~, k_m(m)] = min(abs(freq_axis_shifted - Baseband_Freqs(m)));
% end

k_m = zeros(1, length(ActiveFreqs));
for m = 1:length(ActiveFreqs)
    [~, k_m(m)] = min(abs(freq_axis_shifted - ActiveFreqs(m)));
end

% 验证频率对应关系
fprintf('频率对应验证:\n');
for m = 1:length(ActiveFreqs )
    actual_freq = freq_axis_shifted(k_m(m));
    fprintf('  信道%d: 理论=%.3fMHz, 实际=%.3fMHz, 索引=%d\n', ...
            m, ActiveFreqs (m)/1e6, actual_freq/1e6, k_m(m));
end


% % 验证频率对应关系
% fprintf('频率对应验证:\n');
% for m = 1:length(Baseband_Freqs)
%     actual_freq = freq_axis_shifted(k_m(m));
%     fprintf('  信道%d: 理论=%.3fMHz, 实际=%.3fMHz, 索引=%d\n', ...
%             m, Baseband_Freqs(m)/1e6, actual_freq/1e6, k_m(m));
% end

% 确保索引在有效范围内
k_m = max(1, min(TotalSamples, k_m));

% 初始化结果矩阵
ML_Scores = zeros(GridSize, GridSize);

for i = 1:GridSize
    for j = 1:GridSize
        TestPos = [x_grid(i), y_grid(j)];
        
        % 计算到各站的时延
        TestTDOA = zeros(NumStations, 1);
        for l = 1:NumStations
            TestTDOA(l) = (norm(StationPos(l,:) - TestPos) - norm(StationPos(1,:) - TestPos)) / c;
        end
        
        % 构建矩阵 Y
        Y = zeros(length(k_m), NumStations);
        for l = 1:NumStations
            % 计算DFT
            R_l = fft(RxSignals(l, :));
            
            % 提取跳频频点并补偿时延
            for m = 1:length(k_m)
                % 使用基带频率进行相位补偿
                PhaseCompensation = exp(1j * 2 * pi * Baseband_Freqs(m) * TestTDOA(l));
                Y(m, l) = R_l(k_m(m)) * PhaseCompensation;
            end
        end
        
        % 计算互相关矩阵和最大特征值
        B = Y' * Y;
        EigenValues = eig(B);
        ML_Scores(i, j) = max(real(EigenValues));
    end
end

% 找到估计位置
[~, max_idx] = max(ML_Scores(:));
[i_est, j_est] = ind2sub(size(ML_Scores), max_idx);
ML_EstPos = [x_grid(i_est), y_grid(j_est)];

fprintf('原始算法定位结果: (%.1f, %.1f)m\n', ML_EstPos(1), ML_EstPos(2));

%% 5. 算法二：信道化融合定位算法（基带等效版本）
fprintf('执行信道化融合定位算法...\n');

% 信道化滤波器参数
NumChannels_Filterbank = 20;  % 滤波器组信道数（可以大于跳频信道数）
ChannelBandwidth_Total = Fs;  % 总带宽 = 采样率
ChannelBandwidth_Individual = ChannelBandwidth_Total / NumChannels_Filterbank;

% 设计均匀分布的信道化滤波器组（低通原型滤波器）
NumTaps = 128;
PrototypeFilter = fir1(NumTaps-1, 1/NumChannels_Filterbank);  % 低通原型滤波器

fprintf('设计完成 %d 信道滤波器组，每个带宽 %.1f kHz\n', ...
    NumChannels_Filterbank, ChannelBandwidth_Individual/1e3);


% 对每个观测站进行处理
ActiveChannels = cell(1, NumStations);
ChannelOutputs_Baseband = cell(NumStations, NumChannels_Filterbank);  % 基带数据
ChannelFrequencies = linspace(-Fs/2 + ChannelBandwidth_Individual/2, ...
                              Fs/2 - ChannelBandwidth_Individual/2, ...
                              NumChannels_Filterbank);


% 存储每个站每个信道的检测信息
DetectionThresholds = zeros(NumStations, NumChannels_Filterbank);
NoiseFloors = zeros(NumStations, NumChannels_Filterbank);
SignalPowers = zeros(NumStations, NumChannels_Filterbank);
DetectionStatistics = zeros(NumStations, NumChannels_Filterbank);
StationActiveChannels = false(NumStations, NumChannels_Filterbank);  % 每个信道的检测结果

for l = 1:NumStations
    fprintf('处理站%d的信号...\n', l);
    
     for ch = 1:NumChannels_Filterbank
        % === 正确的信道化处理流程 ===
        
        % 1. 下变频到基带（先下变频）
        Fc = ChannelFrequencies(ch);
        t = (0:length(RxSignals(l,:))-1) / Fs;
        Downconverter = exp(-1j * 2 * pi * Fc * t);
        DownconvertedSignal = RxSignals(l, :) .* Downconverter;
        
        % 2. 低通滤波（使用原型滤波器）
        LowpassFiltered = filter(PrototypeFilter, 1, DownconvertedSignal);
        
        % 3. 降采样
        DecimationFactor = NumChannels_Filterbank;
        DecimatedSignal = LowpassFiltered(1:DecimationFactor:end);
        
        % 存储基带信号
        ChannelOutputs_Baseband{l, ch} = DecimatedSignal;

    end   

   % === 使用频域能量检测器进行自适应门限检测 ===
    fprintf('  站%d使用频域能量检测器:\n', l);
    
    for ch = 1:NumChannels_Filterbank
        % 获取当前信道的基带IQ数据
        baseband_iq = ChannelOutputs_Baseband{l, ch};
        
        % 使用频域能量检测器进行检测
        [detection_result, signal_power, noise_power, detection_statistic,threshold] = ...
            frequency_domain_energy_detector(baseband_iq, Fs/NumChannels_Filterbank, SymbolRate);
        
        DetectionThresholds(l, ch) = threshold;
        NoiseFloors(l, ch) = noise_power;
        SignalPowers(l, ch) = signal_power;
        StationActiveChannels(l,ch) = detection_result;
        DetectionStatistics(l, ch) = detection_statistic;
        
        % 输出检测详情
        if detection_result
            fprintf('    信道%d: 检测到信号! 信号功率=%.6f, 噪声=%.6f, 统计量=%.1fdB\n', ...
                ch, signal_power, noise_power, detection_statistic);
        else
            if mod(ch, 5) == 1  % 减少输出，每5个信道输出一个
                fprintf('    信道%d: 无信号, 信号功率=%.6f, 噪声=%.6f, 统计量=%.1fdB\n', ...
                    ch, signal_power, noise_power, 10*log10(detection_statistic));
            end
        end
    end
    
    ActiveChannels{l} = find(StationActiveChannels(l,:));
    fprintf('  站%d检测到 %d 个活跃信道: %s\n', ...
        l, length(ActiveChannels{l}), mat2str(ActiveChannels{l}));

    % 可选：绘制能量检测结果
    if l == 1
        figure;
        plot(1:NumChannels_Filterbank, DetectionStatistics(l,:), 'b.-', 'LineWidth', 1.5);
        hold on;
        plot(1:NumChannels_Filterbank, DetectionThresholds(l,:), 'r--', 'LineWidth', 2);
        plot(ActiveChannels{l}, DetectionStatistics(l,ActiveChannels{l}), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        xlabel('信道编号');
        ylabel('能量 (dB)');
        title(sprintf('站%d - 信道能量检测', l));
        legend('信道能量', '检测门限', '活跃信道', 'Location', 'best');
        grid on;
    end
end

% 确定所有站共有的活跃信道
CommonActiveChannels = ActiveChannels{1};
for l = 2:NumStations
    CommonActiveChannels = intersect(CommonActiveChannels, ActiveChannels{l});
end

fprintf('所有站共有的活跃信道: %s\n', mat2str(CommonActiveChannels));

if isempty(CommonActiveChannels)
    fprintf('警告：未检测到所有站共有的活跃信道，使用第一个站的检测结果\n');
    CommonActiveChannels = ActiveChannels{1};
end

% 获取活跃信道对应的物理频率
ActiveChannelFreqs = ChannelFrequencies(CommonActiveChannels);


fprintf('活跃信道对应的基带频率: ');
for i = 1:length(ActiveChannelFreqs)
    fprintf('%.3fMHz ', ActiveChannelFreqs(i)/1e6);
end
fprintf('\n');

% 验证与期望频率的匹配
ExpectedFreqs = Baseband_Freqs(HopPattern);
fprintf('期望的跳频频率: ');
for i = 1:length(ExpectedFreqs)
    fprintf('%.3fMHz ', ExpectedFreqs(i)/1e6);
end
fprintf('\n');

% 计算频率匹配度
MatchScores = zeros(length(ExpectedFreqs), 1);
for i = 1:length(ExpectedFreqs)
    [~, MatchScores(i)] = min(abs(ActiveChannelFreqs - ExpectedFreqs(i)));
end
fprintf('频率匹配误差: %.3f MHz\n', mean(abs(ActiveChannelFreqs(MatchScores) - ExpectedFreqs))/1e6);

% 初始化结果矩阵
Fusion_Scores = zeros(GridSize, GridSize);
TotalSamples_Decimated = length(ChannelOutputs_Baseband{1,1});
fprintf('开始网格搜索...\n');
for i = 1:GridSize
    if mod(i, 10) == 0
        fprintf('  进度: %d/%d\n', i, GridSize);
    end
    for j = 1:GridSize
        TestPos = [x_grid(i), y_grid(j)];
        
        % 计算到各站的时延
        TestTDOA = zeros(NumStations, 1);
        for l = 1:NumStations
            TestTDOA(l) = (norm(StationPos(l,:) - TestPos) - norm(StationPos(1,:) - TestPos)) / c;
        end
        
        % 构建矩阵 Y（使用基带信道化输出）
        NumActiveChannels = length(CommonActiveChannels);
        Y_fusion = zeros(NumActiveChannels, NumStations);
        
        for l = 1:NumStations
            for idx = 1:NumActiveChannels
                ch = CommonActiveChannels(idx);
                
                % 获取该信道的基带数据（已经是零频）
                BasebandData = ChannelOutputs_Baseband{l, ch};
                
                % 计算复包络
                Z_l_m = sum(BasebandData);
                
                % 时延补偿（使用该信道的原始物理频率）
                PhaseCompensation = exp(1j * 2 * pi * ActiveChannelFreqs(idx) * TestTDOA(l));
                Y_fusion(idx, l) = Z_l_m * PhaseCompensation;
            end
        end
        
        % 计算互相关矩阵和最大特征值
        B_fusion = Y_fusion' * Y_fusion;
        EigenValues_fusion = eig(B_fusion);
        Fusion_Scores(i, j) = max(real(EigenValues_fusion));
    end
end

% 找到估计位置
[~, max_idx_fusion] = max(Fusion_Scores(:));
[i_est_fusion, j_est_fusion] = ind2sub(size(Fusion_Scores), max_idx_fusion);
Fusion_EstPos = [x_grid(i_est_fusion), y_grid(j_est_fusion)];

fprintf('信道化融合定位结果: (%.1f, %.1f)m\n', Fusion_EstPos(1), Fusion_EstPos(2));
%% 6. 性能评估与可视化
fprintf('性能评估与可视化...\n');

% 计算定位误差
ML_Error = norm(ML_EstPos - TargetPos);
Fusion_Error = norm(Fusion_EstPos - TargetPos);

fprintf('定位误差:\n');
fprintf('  原始算法: %.2f m\n', ML_Error);
fprintf('  融合算法: %.2f m\n', Fusion_Error);

% 绘制结果
figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
contourf(x_grid, y_grid, ML_Scores', 20);
hold on;
plot(TargetPos(1), TargetPos(2), 'rp', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '真实位置');
plot(ML_EstPos(1), ML_EstPos(2), 'wx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '估计位置');
plot(StationPos(:,1), StationPos(:,2), '^', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 8, 'DisplayName', '观测站');
title(sprintf('原始直接定位算法 (误差: %.1fm)', ML_Error));
xlabel('X (m)'); ylabel('Y (m)');
legend; colorbar; axis equal;

subplot(1, 2, 2);
contourf(x_grid, y_grid, Fusion_Scores', 20);
hold on;
plot(TargetPos(1), TargetPos(2), 'rp', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '真实位置');
plot(Fusion_EstPos(1), Fusion_EstPos(2), 'wx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '估计位置');
plot(StationPos(:,1), StationPos(:,2), '^', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 8, 'DisplayName', '观测站');
title(sprintf('信道化融合定位算法 (误差: %.1fm)', Fusion_Error));
xlabel('X (m)'); ylabel('Y (m)');
legend; colorbar; axis equal;

% 绘制信号频谱
figure;
R1_spectrum = fft(RxSignals(1, :));
freq_axis = (0:TotalSamples-1) * Fs / TotalSamples / 1e6;
plot(freq_axis, 20*log10(abs(fftshift(R1_spectrum))));
xlabel('频率 (MHz)'); ylabel('幅度 (dB)');
title('基带等效接收信号频谱 (站1)');
grid on;



fprintf('基带等效模型仿真完成!\n');




%% 频域能量检测器函数
function [detection_result, signal_power, noise_power, detection_statistic,threshold] = frequency_domain_energy_detector(iq_data, fs, signal_bandwidth,threshold)
% QPSK信号频域能量检测器
% 输入参数：
%   iq_data: 输入的IQ数据
%   fs: 采样率
%   false_alarm_prob: 虚警概率（默认0.01）
%   signal_bandwidth: 预期信号带宽（默认fs/4）
% 输出参数：
%   detection_result: 检测结果
%   signal_power: 信号功率估计
%   noise_power: 噪声功率估计
%   threshold: 检测门限
%   freq_axis: 频率轴
%   psd_est: 功率谱密度估计

    if nargin < 3
        signal_bandwidth = fs/4;
    end
    if nargin < 4
        threshold  = 2;  % 频域信噪比
    end   

    iq_data = iq_data(:);


    % 计算功率谱密度（PSD）
    [psd_est, freq_axis] = pwelch(iq_data, hamming(256), 128, 1024, fs, 'centered');
    psd_est = psd_est / max(psd_est); % 归一化

    % figure
    % plot(freq_axis/1e3,psd_est)
    % hold on
    % plot(freq_axis(peak_indices(1))/1e3,psd_est(peak_indices(1)),'ro')

    % 寻找信号可能存在的频带（峰值检测）
    [peak_values, peak_indices] = findpeaks(psd_est, 'SortStr', 'descend', 'NPeaks', 3);

    if isempty(peak_indices)
        detection_result = false;
        signal_power = 0;
        noise_power = mean(psd_est);
        threshold = 0;
        return;
    end

    % 确定主峰对应的中心频率
    main_peak_idx = peak_indices(1);
    center_freq = freq_axis(main_peak_idx);

    % 定义信号带宽范围
    freq_resolution = fs / 1024; % 频率分辨率
    signal_bw_bins = ceil(signal_bandwidth / freq_resolution / 2);

    % 信号频带索引
    signal_band_start = max(1, main_peak_idx - signal_bw_bins);
    signal_band_end = min(length(freq_axis), main_peak_idx + signal_bw_bins);
    signal_band_indices = signal_band_start:signal_band_end;

    % 噪声频带索引（除去信号频带）
    noise_band_indices = setdiff(1:length(freq_axis), signal_band_indices);

    % 估计噪声功率（使用噪声频带的PSD）
    noise_psd_values = psd_est(noise_band_indices);
    noise_power = mean(noise_psd_values);

    % 估计信号功率（信号频带总功率减去噪声功率）
    signal_band_power = mean(psd_est(signal_band_indices));
    signal_power = max(0, signal_band_power - noise_power);

    % 计算检测统计量（信号频带平均功率与噪声功率之比）
    detection_statistic = 10*log10(signal_band_power / (noise_power + eps));

    % 判决,频域信噪比大于-5dB
    detection_result = detection_statistic > threshold;
end