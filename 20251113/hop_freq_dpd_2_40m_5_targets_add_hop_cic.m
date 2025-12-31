clear; close all; clc;

%% 1. 参数设置模块 - 多目标定位
fprintf('设置多目标定位系统参数...\n');

% 跳频信号参数
T_hop = 10e-3;          % 每跳持续时间 10ms
SymbolRate = 50e3;      % 符号速率 50k Baud
Modulation = 'QPSK';
NumHopsPerTarget = 4;   % 每个目标每时间片的跳数
NumTargets = 2;         % 目标数量
NumTimeSegments = 1;    % 时间片段数量

% 跳频频点定义
IF_CenterFreqs_MHz = [48.16, 48.48, 48.80, 49.12, 49.44, 49.76, 50.08, 50.40, ...
                      50.72, 51.04, 51.36, 51.68, 52.00, 52.32, 52.64, 52.96, ...
                      53.28, 53.60, 53.92, 54.24];

IF_CenterRFFreqs_MHz =  IF_CenterFreqs_MHz + (1180-51.2);
fprintf('%.2f\n', IF_CenterRFFreqs_MHz);  % 输出: 1000000.00
IF_Center = 51.20e6;    % 中频中心频率
Baseband_Freqs = (IF_CenterFreqs_MHz * 1e6) - IF_Center;
NumChannels = 20;
 
% % 方法1:定义多目标跳频图案
% HopPatterns = {
%     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];  % 目标1
%     [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1, 2, 3, 4];  % 目标2
%     [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1, 2, 3, 4, 5, 6, 7, 8];  % 目标3
%     [13, 14, 15, 16, 17, 18, 19, 20, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];  % 目标4
%     [17, 18, 19, 20, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]   % 目标5
% };

% 方法2: 或者使用您原来的设计但支持任意时间片段
HopPatterns_Flexible = cell(NumTargets, 1);
for t = 1:NumTargets
    % 生成足够长的跳频序列
    total_hops_needed = NumTimeSegments * NumHopsPerTarget;
    full_pattern = zeros(1, total_hops_needed);
    
    % 使用循环序列填充
    for i = 1:total_hops_needed
        full_pattern(i) = mod((t-1)*4 + i - 1, NumChannels) + 1;
    end
    
    HopPatterns_Flexible{t} = full_pattern;
end
% 选择使用哪种方法（更灵活）
HopPatterns = HopPatterns_Flexible;

fprintf('各目标跳频图案:\n');
for t = 1:NumTargets
    fprintf('  目标%d: %s\n', t, mat2str(HopPatterns{t}(1:NumHopsPerTarget)));
end



% 接收机参数
Fs = 6.4e6;             % 基带采样率 6.4MHz
NumStations = 5;        % 观测站数量

% 多目标位置定义
TargetPositions = [
    100, -100;      % 目标1
    -200, -200;      % 目标2
    200, 200;       % 目标3
    -100, -200;     % 目标4
    0, 0            % 目标5
];

StationPos = [
    15e3, 2e3;      % 站1
    -12e3, 8e3;     % 站2
    5e3, -18e3;     % 站3
    -8e3, -6e3;     % 站4
    10e3, -15e3     % 站5
];

% 计算各目标的TDOA
c = 3e8;
TrueTDOAs = zeros(NumTargets, NumStations);
for t = 1:NumTargets
    distances = vecnorm(StationPos - TargetPositions(t,:), 2, 2);
    for l = 1:NumStations
        TrueTDOAs(t, l) = (distances(l) - distances(1)) / c;
    end
end

% 算法参数
SNR_dB = -15;
FixedNoisePower = 1;

fprintf('多目标定位参数设置完成.\n');

%% 2. 多目标多时间片段信号生成
fprintf('生成多目标多时间片段信号...\n');

% 计算采样参数
SamplesPerSymbol = floor(Fs / SymbolRate);
SymbolsPerHop = floor(T_hop * SymbolRate);
SamplesPerHop = round(T_hop * Fs);
SamplesPerSegment = NumHopsPerTarget * SamplesPerHop;

% 生成所有目标的信号
AllTargetSignals = cell(NumTargets, NumTimeSegments);

for t = 1:NumTargets
    fprintf('生成目标%d的信号...\n', t);
    
    for seg = 1:NumTimeSegments
        % 当前目标在当前时间片段的跳频图案
        PatternStart = (seg-1) * NumHopsPerTarget + 1;
        PatternEnd = seg * NumHopsPerTarget;
        CurrentPattern = HopPatterns{t}(PatternStart:PatternEnd);
        
        % 生成当前时间片段的信号
        TxSignal_Segment = zeros(1, SamplesPerSegment);
        TimeVector = (0:SamplesPerSegment-1) / Fs;
        
        for h = 1:NumHopsPerTarget
            % 生成随机QPSK符号
            DataSymbols = (randi(4, 1, SymbolsPerHop) - 1);
            QPSKSymbols = exp(1j * (pi/4 + DataSymbols * pi/2));
            
            % 上采样成型
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
            
            % 调制到当前跳频点
            CurrentFreq = Baseband_Freqs(CurrentPattern(h));
            Carrier = exp(1j * 2 * pi * CurrentFreq * TimeVector(HopSamples));
            TxSignal_Segment(HopSamples) = SymbolWaveform .* Carrier;
        end
        
        AllTargetSignals{t, seg} = TxSignal_Segment;

        % 验证生成的信号
        if seg == 1 && t <= 3  % 只验证前3个目标的第一个片段
            signal_power = mean(abs(TxSignal_Segment).^2);
            fprintf('  目标%d片段%d: 使用频点%s, 信号功率=%.6f\n', ...
                t, seg, mat2str(CurrentPattern), signal_power);
        end

    end
end

fprintf('多目标信号生成完成.\n');

%% 3. 多目标信号传播模拟
fprintf('模拟多目标信号传播...\n');

% 接收信号是多个目标信号的叠加
AllRxSignals = cell(NumTimeSegments, NumStations);

for seg = 1:NumTimeSegments
    fprintf('处理时间片段 %d/%d...\n', seg, NumTimeSegments);
    
    for l = 1:NumStations
        CombinedSignal = zeros(1, SamplesPerSegment);
        
        % 叠加所有目标的信号
        for t = 1:NumTargets
            % 获取目标信号
            TargetSignal = AllTargetSignals{t, seg};
            
            % 应用时延
            DelaySamples = round(TrueTDOAs(t, l) * Fs);
            DelayedSignal = circshift(TargetSignal, DelaySamples);
            
            % 调整信号功率
            CurrentSignalPower = mean(abs(DelayedSignal).^2);
            RequiredSignalPower = FixedNoisePower * (10^(SNR_dB/10)) ; % 平均分配功率
            ScaleFactor = sqrt(RequiredSignalPower / CurrentSignalPower);
            ScaledSignal = DelayedSignal * ScaleFactor;
            
            CombinedSignal = CombinedSignal + ScaledSignal;
        end
        
        % 添加噪声
        Noise = sqrt(FixedNoisePower/2) * (randn(1, SamplesPerSegment) + 1j * randn(1, SamplesPerSegment));
        AllRxSignals{seg, l} = CombinedSignal + Noise;
    end
end





% 绘制信号频谱
figure;
R1_spectrum = fft(AllRxSignals{1,1}(1, :));
freq_axis = (0:SamplesPerSegment-1) * Fs / SamplesPerSegment / 1e6;
plot(freq_axis, 20*log10(abs(fftshift(R1_spectrum))));
xlabel('频率 (MHz)'); ylabel('幅度 (dB)');
title('基带等效接收信号频谱 (站1)');
grid on;

fprintf('多目标信号传播模拟完成.\n');

%% 4. 多目标信道化处理与信号分离（使用频域能量检测）
fprintf('执行多目标信道化处理与信号分离（使用频域能量检测）...\n');

% 信道化滤波器参数
NumChannels_Filterbank = 20;
ChannelBandwidth_Total = Fs;
ChannelBandwidth_Individual = ChannelBandwidth_Total / NumChannels_Filterbank;

% 设计滤波器组
NumTaps = 128;
PrototypeFilter = fir1(NumTaps-1, 1/NumChannels_Filterbank);
ChannelFrequencies = linspace(-Fs/2 + ChannelBandwidth_Individual/2, ...
                              Fs/2 - ChannelBandwidth_Individual/2, ...
                              NumChannels_Filterbank);

% 存储信道化输出和检测结果
AllChannelOutputs = cell(NumTimeSegments, NumStations, NumChannels_Filterbank);
ChannelDetectionResults = cell(NumTimeSegments, NumStations);

% 信道化处理和能量检测
for seg = 1:NumTimeSegments
    fprintf('时间片段 %d/%d 信道化处理与能量检测...\n', seg, NumTimeSegments);
    
    for l = 1:NumStations
        RxSignal = AllRxSignals{seg, l};
        ActiveChannels = [];
        
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
            
            AllChannelOutputs{seg, l, ch} = DecimatedSignal;
            
            % 频域能量检测
            [detection_result, signal_power, noise_power, detection_statistic, threshold] = ...
                frequency_domain_energy_detector(DecimatedSignal, Fs/NumChannels_Filterbank, SymbolRate);
            
            if detection_result
                ActiveChannels = [ActiveChannels, ch];
                if l == 1 && seg == 1  % 只输出第一个站第一个片段的检测详情
                    fprintf('    信道%d: 检测到信号! 功率=%.6f, 统计量=%.1fdB\n', ...
                        ch, signal_power, detection_statistic);
                end
            end
        end
        
        ChannelDetectionResults{seg, l} = ActiveChannels;
        fprintf('  站%d在片段%d检测到 %d 个活跃信道: %s\n', ...
            l, seg, length(ActiveChannels), mat2str(ActiveChannels));
    end
end

% 确定每个时间片段所有站共有的活跃信道
CommonActiveChannelsPerSegment = cell(1, NumTimeSegments);
for seg = 1:NumTimeSegments
    CommonChannels = ChannelDetectionResults{seg, 1};
    for l = 2:NumStations
        CommonChannels = intersect(CommonChannels, ChannelDetectionResults{seg, l});
    end
    
    if isempty(CommonChannels)
        fprintf('警告：时间片段%d未检测到所有站共有的活跃信道，使用第一个站的检测结果\n', seg);
        CommonChannels = ChannelDetectionResults{seg, 1};
    end
    
    CommonActiveChannelsPerSegment{seg} = CommonChannels;
    fprintf('时间片段%d的共有活跃信道: %s\n', seg, mat2str(CommonChannels));
end

% 基于检测到的活跃信道进行目标分离
fprintf('基于能量检测结果进行目标分离...\n');



% 构建目标-频点映射关系
TargetChannelMapping = cell(NumTargets, NumTimeSegments);
for t = 1:NumTargets
    for seg = 1:NumTimeSegments
        PatternStart = (seg-1) * NumHopsPerTarget + 1;
        PatternEnd = seg * NumHopsPerTarget;
        TargetChannelMapping{t, seg} = HopPatterns{t}(PatternStart:PatternEnd);
    end
end

fprintf('目标-频点映射关系:\n');
for t = 1:NumTargets
    fprintf('  目标%d:\n', t);
    for seg = 1:NumTimeSegments
        expected_channels = TargetChannelMapping{t, seg};
        fprintf('    片段%d: %s\n', seg, mat2str(TargetChannelMapping{t, seg}));

        % 检查频点是否重复（不应该重复）
        if length(unique(expected_channels)) ~= length(expected_channels)
            fprintf('    警告: 频点重复!\n');
        end

    end
end

% 构建目标-频点关联（基于检测结果和跳频图案的匹配）
TargetDetectedChannels = cell(NumTargets, NumTimeSegments);
TargetConfidence = zeros(NumTargets, NumTimeSegments);

for seg = 1:NumTimeSegments
    DetectedChannels = CommonActiveChannelsPerSegment{seg};
    
    % 为每个目标计算信道匹配度
    for t = 1:NumTargets
        ExpectedChannels = TargetChannelMapping{t, seg};
        
        % 计算匹配的信道
        MatchedChannels = intersect(ExpectedChannels, DetectedChannels);
        TargetDetectedChannels{t, seg} = MatchedChannels;
        
        % 计算匹配置信度
        TargetConfidence(t, seg) = length(MatchedChannels) / length(ExpectedChannels);
        
        fprintf('  目标%d在片段%d: 期望%s, 检测%s, 匹配度%.1f%%\n', ...
            t, seg, mat2str(ExpectedChannels), mat2str(MatchedChannels), ...
            TargetConfidence(t, seg) * 100);
    end
end

% 确定可定位的目标（匹配度高于阈值）
LocatableTargets = [];
ConfidenceThreshold = 0.5;  % 50%匹配度阈值

for t = 1:NumTargets
    avg_confidence = mean(TargetConfidence(t, :));
    if avg_confidence >= ConfidenceThreshold
        LocatableTargets = [LocatableTargets, t];
        fprintf('目标%d平均匹配度%.1f%%，可进行定位\n', t, avg_confidence * 100);
    else
        fprintf('目标%d平均匹配度%.1f%%，匹配度不足\n', t, avg_confidence * 100);
    end
end

if isempty(LocatableTargets)
    fprintf('警告：没有检测到可定位的目标，使用所有检测到的信道进行定位\n');
    LocatableTargets = 1:NumTargets;
end

fprintf('可定位的目标: %s\n', mat2str(LocatableTargets));

%% 5. 多目标定位算法（基于能量检测结果）
fprintf('执行基于能量检测的多目标定位...\n');

SearchRange = 1000;
GridSize = 21;
x_grid = linspace(-SearchRange, SearchRange, GridSize);
y_grid = linspace(-SearchRange, SearchRange, GridSize);

EstimatedPositions = zeros(NumTargets, 2);
TargetScores = cell(1, NumTargets);
TargetLocatable = false(1, NumTargets);

for t = 1:NumTargets
    if ~ismember(t, LocatableTargets)
        fprintf('跳过目标%d（匹配度不足）\n', t);
        TargetLocatable(t) = false;
        continue;
    end
    
    fprintf('定位目标%d（基于能量检测结果）...\n', t);
    TargetLocatable(t) = true;
    
    Target_Scores = zeros(GridSize, GridSize);
    
    for i = 1:GridSize
        if mod(i, 5) == 0
            fprintf('  目标%d进度: %d/%d\n', t, i, GridSize);
        end
        
        for j = 1:GridSize
            TestPos = [x_grid(i), y_grid(j)];
            
            % 计算到各站的时延
            TestTDOA = zeros(NumStations, 1);
            for l = 1:NumStations
                TestTDOA(l) = (norm(StationPos(l,:) - TestPos) - norm(StationPos(1,:) - TestPos)) / c;
            end
            
            % 构建该目标的观测矩阵（使用检测到的活跃信道）
            TotalChannelsForTarget = 0;
            for seg = 1:NumTimeSegments
                TotalChannelsForTarget = TotalChannelsForTarget + length(TargetDetectedChannels{t, seg});
            end
            
            if TotalChannelsForTarget == 0
                Target_Scores(i, j) = 0;
                continue;
            end
            
            Y_target = zeros(TotalChannelsForTarget, NumStations);
            current_row = 1;
            
            for seg = 1:NumTimeSegments
                DetectedChannels = TargetDetectedChannels{t, seg};
                
                for l = 1:NumStations
                    for idx = 1:length(DetectedChannels)
                        ch = DetectedChannels(idx);
                        
                        % 获取该信道的基带数据
                        BasebandData = AllChannelOutputs{seg, l, ch};
                        
                        % 计算复包络
                        Z_l_m = sum(BasebandData);
                        
                        % 时延补偿
                        FreqIndex = ch;  % 信道编号对应频率索引
                        PhaseCompensation = exp(1j * 2 * pi * Baseband_Freqs(FreqIndex) * TestTDOA(l));
                        Y_target(current_row + idx - 1, l) = Z_l_m * PhaseCompensation;
                    end
                end
                current_row = current_row + length(DetectedChannels);
            end
            
            % 计算互相关矩阵和最大特征值
            if rank(Y_target) >= NumStations - 1
                B_target = Y_target' * Y_target;
                EigenValues_target = eig(B_target);
                Target_Scores(i, j) = max(real(EigenValues_target));
            else
                Target_Scores(i, j) = 0;
            end
        end
    end
    
    % 找到估计位置
    [max_score, max_idx] = max(Target_Scores(:));
    [i_est, j_est] = ind2sub(size(Target_Scores), max_idx);
    EstimatedPositions(t, :) = [x_grid(i_est), y_grid(j_est)];
    
    TargetScores{t} = Target_Scores;
    
    fprintf('  目标%d定位完成: (%.1f, %.1f)m, 最大得分: %.2f\n', ...
        t, EstimatedPositions(t, 1), EstimatedPositions(t, 2), max_score);
end


%% 6. 性能评估与可视化
fprintf('多目标定位性能评估...\n');

% 计算各目标的定位误差
PositionErrors = zeros(NumTargets, 1);
for t = 1:NumTargets
    PositionErrors(t) = norm(EstimatedPositions(t, :) - TargetPositions(t, :));
    fprintf('目标%d: 真实位置(%.1f, %.1f), 估计位置(%.1f, %.1f), 误差=%.2fm\n', ...
        t, TargetPositions(t, 1), TargetPositions(t, 2), ...
        EstimatedPositions(t, 1), EstimatedPositions(t, 2), PositionErrors(t));
end

fprintf('平均定位误差: %.2fm\n', mean(PositionErrors));

% 绘制多目标定位结果
figure('Position', [100, 100, 1200, 800]);

% 绘制所有目标的定位结果
colors = lines(NumTargets);
for t = 1:NumTargets
    subplot(2, 3, t);
    contourf(x_grid, y_grid, 10*log10(TargetScores{t}'), 20);
    hold on;
    
    % 真实位置
    plot(TargetPositions(t, 1), TargetPositions(t, 2), 'ro', ...
        'MarkerSize', 10, 'LineWidth', 3, 'DisplayName', '真实位置');
    
    % 估计位置
    plot(EstimatedPositions(t, 1), EstimatedPositions(t, 2), 'wx', ...
        'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', '估计位置');
    
    % 观测站
    plot(StationPos(:,1), StationPos(:,2), '^', 'Color', [0.5, 0.5, 0.5], ...
        'MarkerSize', 8, 'DisplayName', '观测站');
    
    title(sprintf('目标%d定位 (误差: %.1fm)', t, PositionErrors(t)));
    xlabel('X (m)'); ylabel('Y (m)');
    legend; colorbar; axis equal;
    xlim([-1000 1000])
    ylim([-1000 1000])
    grid on;
end

% 绘制总体定位结果
subplot(2, 3, 6);
for t = 1:NumTargets
    % 真实位置
    plot(TargetPositions(t, 1), TargetPositions(t, 2), 'o', 'Color', colors(t,:), ...
        'MarkerSize', 10, 'LineWidth', 3, 'DisplayName', sprintf('目标%d真实', t));
    hold on;
    
    % 估计位置
    plot(EstimatedPositions(t, 1), EstimatedPositions(t, 2), 'x', 'Color', colors(t,:), ...
        'MarkerSize', 12, 'LineWidth', 3, 'DisplayName', sprintf('目标%d估计', t));
    
    % 连接真实和估计位置
    plot([TargetPositions(t, 1), EstimatedPositions(t, 1)], ...
         [TargetPositions(t, 2), EstimatedPositions(t, 2)], '--', ...
         'Color', colors(t,:), 'LineWidth', 1);
end

% 观测站
plot(StationPos(:,1), StationPos(:,2), 'k^', 'MarkerSize', 10, ...
    'LineWidth', 2, 'DisplayName', '观测站');

title('多目标定位总体结果');
xlabel('X (m)'); ylabel('Y (m)');
legend('Location', 'best'); grid on; axis equal;

% 绘制频点使用情况
figure;
for t = 1:NumTargets
    all_channels = [];
    for seg = 1:NumTimeSegments
        all_channels = [all_channels, TargetChannelMapping{t, seg}];
    end
    plot(all_channels, t * ones(size(all_channels)), 's-', 'Color', colors(t,:), ...
        'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', sprintf('目标%d', t));
    hold on;
end
xlabel('频点编号');
ylabel('目标编号');
title('多目标频点使用情况');
legend('Location', 'best');
grid on;

fprintf('多目标定位仿真完成!\n');

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
        threshold  = 0;  % 频域信噪比
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