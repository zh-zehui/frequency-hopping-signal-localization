clear; close all; clc;

%% 1. 参数设置模块 - 基带等效模型
fprintf('设置系统参数（基带等效模型）...\n');

% 跳频信号参数
T_hop = 10e-3;          % 每跳持续时间 10ms
SymbolRate = 50e3;      % 符号速率 50k Baud
Modulation = 'QPSK';
NumHops = 4;            % 每个时间片的跳数 (40ms / 10ms)
NumChannels = 20;       % 频点数量
NumTimeSegments = 5;    % 时间片段数量（N=6）

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
Fs = 6.4e6;             % 基带采样率 6.4MHz
NumStations = 5;        % 观测站数量

% 定位场景参数
TargetPos = [100, -100];
StationPos = [
    15e3, 2e3;      % 东北方向，距离约15.1km
    -12e3, 8e3;     % 西北方向，距离约14.4km  
    5e3, -18e3;     % 东南方向，距离约18.7km
    -8e3, -6e3;     % 西南方向，距离约10.0km
    10e3, -15e3     % 东南偏东方向
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
FixedNoisePower = 1;  % 固定噪声功率为1

fprintf('参数设置完成.\n');

%% 2. 多时间片段信号生成
fprintf('生成多时间片段基带等效跳频QPSK信号...\n');

% 每跳符号数
SamplesPerSymbol = floor(Fs / SymbolRate);
SymbolsPerHop = floor(T_hop * SymbolRate);
SamplesPerHop = round(T_hop * Fs);

% 每个时间片段的采样点数
SamplesPerSegment = NumHops * SamplesPerHop;

% 生成多时间片段的跳频图案
HopPatterns = zeros(NumTimeSegments, NumHops);
for seg = 1:NumTimeSegments
    % 循环使用所有频点
    start_channel = mod((seg-1)*NumHops, NumChannels) + 1;
    HopPatterns(seg, :) = mod((start_channel:start_channel+NumHops-1)-1, NumChannels) + 1;
end

fprintf('多时间片段跳频图案:\n');
for seg = 1:NumTimeSegments
    fprintf('  片段%d: %s\n', seg, mat2str(HopPatterns(seg, :)));
end

% 生成所有时间片段的信号
AllTxSignals = cell(1, NumTimeSegments);

for seg = 1:NumTimeSegments
    % 当前时间片段的信号
    TxSignal_Segment = zeros(1, SamplesPerSegment);
    TimeVector = (0:SamplesPerSegment-1) / Fs;
    
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
        CurrentFreq = Baseband_Freqs(HopPatterns(seg, h));
        Carrier = exp(1j * 2 * pi * CurrentFreq * TimeVector(HopSamples));
        TxSignal_Segment(HopSamples) = SymbolWaveform .* Carrier;
    end
    
    AllTxSignals{seg} = TxSignal_Segment;
end

fprintf('多时间片段信号生成完成.\n');

%% 3. 多时间片段信道模拟
fprintf('模拟多时间片段信号传播...\n');

AllRxSignals = cell(NumTimeSegments, NumStations);

for seg = 1:NumTimeSegments
    fprintf('处理时间片段 %d/%d...\n', seg, NumTimeSegments);
    
    for l = 1:NumStations
        % 计算时延对应的采样点数
        DelaySamples = round(TrueTDOA(l) * Fs);
        
        % 应用时延
        DelayedSignal = circshift(AllTxSignals{seg}, DelaySamples);
        
        % 调整信号功率以达到目标SNR
        CurrentSignalPower = mean(abs(DelayedSignal).^2);
        RequiredSignalPower = FixedNoisePower * (10^(SNR_dB/10));
        ScaleFactor = sqrt(RequiredSignalPower / CurrentSignalPower);
        ScaledSignal = DelayedSignal * ScaleFactor;

        % 添加固定功率的复高斯白噪声
        Noise = sqrt(FixedNoisePower/2) * (randn(1, SamplesPerSegment) + 1j * randn(1, SamplesPerSegment));
        
        AllRxSignals{seg, l} = ScaledSignal + Noise;

        % 验证实际的SNR
        ActualSignalPower = mean(abs(ScaledSignal).^2);
        ActualNoisePower = mean(abs(Noise).^2);
        ActualSNR = 10 * log10(ActualSignalPower / ActualNoisePower);
        
        if seg == 1
            fprintf('  站%d: 目标SNR=%.1fdB, 实际SNR=%.1fdB\n', l, SNR_dB, ActualSNR);
        end
    end
end

fprintf('多时间片段信道模拟完成.\n');

%% 4. 多时间片段信道化处理
fprintf('执行多时间片段信道化处理...\n');

% 信道化滤波器参数
NumChannels_Filterbank = 20;  % 滤波器组信道数
ChannelBandwidth_Total = Fs;
ChannelBandwidth_Individual = ChannelBandwidth_Total / NumChannels_Filterbank;

% 设计均匀分布的信道化滤波器组
NumTaps = 128;
PrototypeFilter = fir1(NumTaps-1, 1/NumChannels_Filterbank);

% 信道频率定义
ChannelFrequencies = linspace(-Fs/2 + ChannelBandwidth_Individual/2, ...
                              Fs/2 - ChannelBandwidth_Individual/2, ...
                              NumChannels_Filterbank);

% 存储所有时间片段和所有站的信道化输出
AllChannelOutputs = cell(NumTimeSegments, NumStations, NumChannels_Filterbank);
AllActiveChannels = cell(NumTimeSegments, NumStations);

% 对每个时间片段和每个观测站进行处理
for seg = 1:NumTimeSegments
    fprintf('处理时间片段 %d/%d 的信道化...\n', seg, NumTimeSegments);
    
    for l = 1:NumStations
        % 获取当前时间片段和站的接收信号
        RxSignal = AllRxSignals{seg, l};
        
        % 信道化处理
        for ch = 1:NumChannels_Filterbank
            % 下变频到基带
            Fc = ChannelFrequencies(ch);
            t = (0:length(RxSignal)-1) / Fs;
            Downconverter = exp(-1j * 2 * pi * Fc * t);
            DownconvertedSignal = RxSignal .* Downconverter;
            
            % 低通滤波
            LowpassFiltered = filter(PrototypeFilter, 1, DownconvertedSignal);
            
            % 降采样
            DecimationFactor = NumChannels_Filterbank;
            DecimatedSignal = LowpassFiltered(1:DecimationFactor:end);
            
            % 存储基带信号
            AllChannelOutputs{seg, l, ch} = DecimatedSignal;
        end
        
        % 信道检测
        ActiveChannels_Segment = [];
        for ch = 1:NumChannels_Filterbank
            baseband_iq = AllChannelOutputs{seg, l, ch};
            
            % 使用频域能量检测器
            [detection_result, ~, ~, detection_statistic, ~] = ...
                frequency_domain_energy_detector(baseband_iq, Fs/NumChannels_Filterbank, SymbolRate);
            
            if detection_result
                ActiveChannels_Segment = [ActiveChannels_Segment, ch];
            end
        end
        
        AllActiveChannels{seg, l} = ActiveChannels_Segment;
        
        if seg == 1
            fprintf('  站%d检测到 %d 个活跃信道: %s\n', ...
                l, length(ActiveChannels_Segment), mat2str(ActiveChannels_Segment));
        end
    end
end

% 确定每个时间片段所有站共有的活跃信道
CommonActiveChannelsPerSegment = cell(1, NumTimeSegments);
AllActiveChannelFreqs = cell(1, NumTimeSegments);

for seg = 1:NumTimeSegments
    CommonChannels = AllActiveChannels{seg, 1};
    for l = 2:NumStations
        CommonChannels = intersect(CommonChannels, AllActiveChannels{seg, l});
    end
    
    if isempty(CommonChannels)
        fprintf('警告：时间片段%d未检测到所有站共有的活跃信道，使用第一个站的检测结果\n', seg);
        CommonChannels = AllActiveChannels{seg, 1};
    end
    
    CommonActiveChannelsPerSegment{seg} = CommonChannels;
    AllActiveChannelFreqs{seg} = ChannelFrequencies(CommonChannels);
    
    fprintf('时间片段%d的共有活跃信道: %s\n', seg, mat2str(CommonChannels));
end

fprintf('多时间片段信道化处理完成.\n');

%% 5. 多时间片段融合定位算法
fprintf('执行多时间片段融合定位算法...\n');

% 搜索网格参数
SearchRange = 1000;
GridSize = 21;
x_grid = linspace(-SearchRange, SearchRange, GridSize);
y_grid = linspace(-SearchRange, SearchRange, GridSize);

% 初始化结果矩阵
Fusion_Scores_MultiSegment = zeros(GridSize, GridSize);

% 获取降采样后的数据长度
TotalSamples_Decimated = length(AllChannelOutputs{1,1,1});

fprintf('开始多时间片段融合网格搜索...\n');

for i = 1:GridSize
    if mod(i, 5) == 0
        fprintf('  进度: %d/%d\n', i, GridSize);
    end
    
    for j = 1:GridSize
        TestPos = [x_grid(i), y_grid(j)];
        
        % 计算到各站的时延（所有时间片段相同）
        TestTDOA = zeros(NumStations, 1);
        for l = 1:NumStations
            TestTDOA(l) = (norm(StationPos(l,:) - TestPos) - norm(StationPos(1,:) - TestPos)) / c;
        end
        
        % 初始化融合矩阵
        TotalActiveChannels = 0;
        for seg = 1:NumTimeSegments
            TotalActiveChannels = TotalActiveChannels + length(CommonActiveChannelsPerSegment{seg});
        end
        
        Y_fusion_multi = zeros(TotalActiveChannels, NumStations);
        current_row = 1;
        
        % 融合所有时间片段的数据
        for seg = 1:NumTimeSegments
            NumActiveChannels = length(CommonActiveChannelsPerSegment{seg});
            ActiveFreqs = AllActiveChannelFreqs{seg};
            
            for l = 1:NumStations
                for idx = 1:NumActiveChannels
                    ch = CommonActiveChannelsPerSegment{seg}(idx);
                    
                    % 获取该信道的基带数据
                    BasebandData = AllChannelOutputs{seg, l, ch};
                    
                    % 计算复包络（使用整个时间片段的数据）
                    Z_l_m = sum(BasebandData);
                    
                    % 时延补偿（使用该信道的原始物理频率）
                    PhaseCompensation = exp(1j * 2 * pi * ActiveFreqs(idx) * TestTDOA(l));
                    Y_fusion_multi(current_row + idx - 1, l) = Z_l_m * PhaseCompensation;
                end
            end
            current_row = current_row + NumActiveChannels;
        end
        
        % 计算互相关矩阵和最大特征值
        B_fusion = Y_fusion_multi' * Y_fusion_multi;
        EigenValues_fusion = eig(B_fusion);
        Fusion_Scores_MultiSegment(i, j) = max(real(EigenValues_fusion));
    end
end

% 找到估计位置
[~, max_idx_fusion] = max(Fusion_Scores_MultiSegment(:));
[i_est_fusion, j_est_fusion] = ind2sub(size(Fusion_Scores_MultiSegment), max_idx_fusion);
Fusion_EstPos_Multi = [x_grid(i_est_fusion), y_grid(j_est_fusion)];

fprintf('多时间片段融合定位结果: (%.1f, %.1f)m\n', Fusion_EstPos_Multi(1), Fusion_EstPos_Multi(2));

%% 6. 单时间片段定位对比（使用第一个时间片段）
fprintf('执行单时间片段定位对比...\n');

% 使用第一个时间片段进行定位
SingleSegment_Scores = zeros(GridSize, GridSize);
seg = 1;  % 使用第一个时间片段

for i = 1:GridSize
    for j = 1:GridSize
        TestPos = [x_grid(i), y_grid(j)];
        
        % 计算到各站的时延
        TestTDOA = zeros(NumStations, 1);
        for l = 1:NumStations
            TestTDOA(l) = (norm(StationPos(l,:) - TestPos) - norm(StationPos(1,:) - TestPos)) / c;
        end
        
        % 构建矩阵 Y（使用第一个时间片段）
        NumActiveChannels = length(CommonActiveChannelsPerSegment{seg});
        Y_single = zeros(NumActiveChannels, NumStations);
        ActiveFreqs = AllActiveChannelFreqs{seg};
        
        for l = 1:NumStations
            for idx = 1:NumActiveChannels
                ch = CommonActiveChannelsPerSegment{seg}(idx);
                
                % 获取该信道的基带数据
                BasebandData = AllChannelOutputs{seg, l, ch};
                
                % 计算复包络
                Z_l_m = sum(BasebandData);
                
                % 时延补偿
                PhaseCompensation = exp(1j * 2 * pi * ActiveFreqs(idx) * TestTDOA(l));
                Y_single(idx, l) = Z_l_m * PhaseCompensation;
            end
        end
        
        % 计算互相关矩阵和最大特征值
        B_single = Y_single' * Y_single;
        EigenValues_single = eig(B_single);
        SingleSegment_Scores(i, j) = max(real(EigenValues_single));
    end
end

% 找到估计位置
[~, max_idx_single] = max(SingleSegment_Scores(:));
[i_est_single, j_est_single] = ind2sub(size(SingleSegment_Scores), max_idx_single);
Single_EstPos = [x_grid(i_est_single), y_grid(j_est_single)];

fprintf('单时间片段定位结果: (%.1f, %.1f)m\n', Single_EstPos(1), Single_EstPos(2));

%% 7. 性能评估与可视化
fprintf('性能评估与可视化...\n');

% 计算定位误差
Single_Error = norm(Single_EstPos - TargetPos);
Multi_Error = norm(Fusion_EstPos_Multi - TargetPos);

fprintf('定位误差对比:\n');
fprintf('  单时间片段(%d个频点): %.2f m\n', length(CommonActiveChannelsPerSegment{1}), Single_Error);
fprintf('  多时间片段融合(%d个频点): %.2f m\n', TotalActiveChannels, Multi_Error);
fprintf('  改善比例: %.1f%%\n', (Single_Error - Multi_Error) / Single_Error * 100);

% 绘制结果对比
figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
contourf(x_grid, y_grid, SingleSegment_Scores', 20);
hold on;
plot(TargetPos(1), TargetPos(2), 'rp', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '真实位置');
plot(Single_EstPos(1), Single_EstPos(2), 'wx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '估计位置');
plot(StationPos(:,1), StationPos(:,2), '^', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 8, 'DisplayName', '观测站');
title(sprintf('单时间片段定位 (误差: %.1fm)', Single_Error));
xlabel('X (m)'); ylabel('Y (m)');
legend; colorbar; axis equal;

subplot(1, 2, 2);
contourf(x_grid, y_grid, Fusion_Scores_MultiSegment', 20);
hold on;
plot(TargetPos(1), TargetPos(2), 'rp', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '真实位置');
plot(Fusion_EstPos_Multi(1), Fusion_EstPos_Multi(2), 'wx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '估计位置');
plot(StationPos(:,1), StationPos(:,2), '^', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 8, 'DisplayName', '观测站');
title(sprintf('多时间片段融合定位 (误差: %.1fm)', Multi_Error));
xlabel('X (m)'); ylabel('Y (m)');
legend; colorbar; axis equal;

% 绘制使用的频点分布
figure;
colors = lines(NumTimeSegments);
for seg = 1:NumTimeSegments
    freqs = AllActiveChannelFreqs{seg} / 1e6;
    plot(freqs, seg * ones(size(freqs)), 'o', 'Color', colors(seg,:), ...
         'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', sprintf('片段%d', seg));
    hold on;
end
xlabel('频率 (MHz)');
ylabel('时间片段');
title('多时间片段使用的频点分布');
legend('Location', 'best');
grid on;

fprintf('多时间片段融合定位仿真完成!\n');

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