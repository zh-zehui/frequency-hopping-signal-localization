function multi_target_signal_sorting(RxSignals, StationPos, Fs, T_hop)
    % 多目标信号分选主函数
    % 输入：多站接收信号，观测站位置，采样率，跳频时间
    % 输出：分选后的信号，目标分配结果
    
    [NumStations, TotalSamples] = size(RxSignals);
    NumHops = TotalSamples / (T_hop * Fs);
    
    %% 步骤1：信道化处理 + 自适应能量检测
    fprintf('=== 步骤1: 信道化处理与能量检测 ===\n');
    
    NumChannels_Filterbank = 20;
    ChannelOutputs_Baseband = cell(NumStations, NumChannels_Filterbank);
    DetectionStatistics = zeros(NumStations, NumChannels_Filterbank);
    StationActiveChannels = false(NumStations, NumChannels_Filterbank);
    for l = 1:NumStations
        [ChannelOutputs_Baseband(l, :), DetectionStatistics(l, :),StationActiveChannels(l, :),ChannelFrequencies] = ...
            channelize_signal(RxSignals(l, :), Fs, NumChannels_Filterbank);
    end
    
    %% 步骤2. 多辐射源信号分选
    fprintf('执行多辐射源信号分选...\n');
    
    [source_clusters, source_features] = multi_source_sorting(...
        ChannelOutputs_Baseband, StationActiveChannels, ChannelFrequencies, Fs, NumStations);
    
    % 可视化分选结果
    figure('Position', [200, 200, 1000, 600]);
    
    % 特征空间可视化
    subplot(2,3,[1,2,4,5]);
    colors = lines(max(source_clusters));
    for i = 1:max(source_clusters)
        cluster_points = source_clusters == i;
        scatter3(source_features.ChannelFeatures(cluster_points, 1), ...
                 source_features.ChannelFeatures(cluster_points, 2), ...
                 source_features.ChannelFeatures(cluster_points, 3), ...
                 50, colors(i,:), 'filled', 'DisplayName', sprintf('辐射源%d', i));
        hold on;
    end
    
    % 噪声点
    noise_points = source_clusters == -1;
    if any(noise_points)
        scatter3(source_features.ChannelFeatures(noise_points, 1), ...
                 source_features.ChannelFeatures(noise_points, 2), ...
                 source_features.ChannelFeatures(noise_points, 3), ...
                 30, [0.5, 0.5, 0.5], 'x', 'DisplayName', '噪声');
    end
    
    xlabel('功率特征'); ylabel('带宽特征'); zlabel('调制特征');
    title('多辐射源信号分选结果');
    legend; grid on;
    
    % 各辐射源的信道分布
    subplot(2,3,3);
    num_sources = max(source_clusters);
    source_channel_counts = zeros(1, num_sources);
    for i = 1:num_sources
        source_channel_counts(i) = sum(source_clusters == i);
    end
    bar(1:num_sources, source_channel_counts);
    xlabel('辐射源编号'); ylabel('信道数量');
    title('各辐射源的信道数');
    
    % 频率分布
    subplot(2,3,6);
    for i = 1:num_sources
        source_channels = source_features.AllActiveChannels(source_clusters == i);
        source_freqs = ChannelFrequencies(source_channels) / 1e6;
        plot(source_freqs, i*ones(size(source_freqs)), 'o', 'MarkerSize', 8, 'LineWidth', 2);
        hold on;
    end
    xlabel('频率 (MHz)'); ylabel('辐射源编号');
    title('各辐射源的频率分布');
    grid on;
    
    fprintf('信号分选完成!\n');
    

end

function [ChannelOutputs_Baseband, DetectionStatistics,StationActiveChannels,ChannelFrequencies] = channelize_signal(signal, Fs, NumChannels)
    % 信道化处理
    ChannelBandwidth_Total = Fs;
    ChannelBandwidth_Individual = ChannelBandwidth_Total / NumChannels;
    
    % 设计低通原型滤波器
    NumTaps = 128;
    PrototypeFilter = fir1(NumTaps-1, 1/NumChannels);
    
    % 预先计算信道频率
    ChannelFrequencies = linspace(-Fs/2 + ChannelBandwidth_Individual/2, ...
                                  Fs/2 - ChannelBandwidth_Individual/2, NumChannels);
    
    ChannelOutputs_Baseband = cell(1, NumChannels);

    for ch = 1:NumChannels
        % 下变频到基带
        Fc = ChannelFrequencies(ch);
        t = (0:length(signal)-1) / Fs;
        Downconverter = exp(-1j * 2 * pi * Fc * t);
        DownconvertedSignal = signal .* Downconverter;
        
        % 低通滤波
        LowpassFiltered = filter(PrototypeFilter, 1, DownconvertedSignal);
        
        % 降采样
        DecimationFactor = NumChannels;
        DecimatedSignal = LowpassFiltered(1:DecimationFactor:end);
        
        ChannelOutputs_Baseband{ch} = DecimatedSignal;

    end

   % === 使用频域能量检测器进行自适应门限检测 ===
    DetectionThresholds = zeros(1, NumChannels);
    StationActiveChannels = false(1, NumChannels);
    DetectionStatistics = zeros(1, NumChannels);
    for ch = 1:NumChannels
        % 获取当前信道的基带IQ数据
        baseband_iq = ChannelOutputs_Baseband{ch};
        
        % 使用频域能量检测器进行检测
        [detection_result, signal_power, noise_power, detection_statistic,threshold] = ...
            frequency_domain_energy_detector(baseband_iq, Fs/NumChannels);
        
        DetectionThresholds(1,ch) = threshold;
        StationActiveChannels(1,ch) = detection_result;
        DetectionStatistics(1,ch) = detection_statistic;
        
        % 输出检测详情
        if detection_result
            fprintf('    信道%d: 检测到信号! 信号功率=%.6f, 噪声=%.6f, 统计量=%.1fdB\n', ...
                ch, signal_power, noise_power, detection_statistic);
        else
            if mod(ch, 5) == 1  % 减少输出，每5个信道输出一个
                fprintf('    信道%d: 无信号, 信号功率=%.6f, 噪声=%.6f, 统计量=%.1fdB\n', ...
                    ch, signal_power, noise_power, detection_statistic);
            end
        end
    end
    
end


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
        threshold  = 5;  % 频域信噪比
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




%% 7. 多辐射源信号分选模块
function [source_clusters, source_features] = multi_source_sorting(ChannelOutputs_Baseband, StationActiveChannels, ChannelFrequencies, Fs, NumStations)
% 多辐射源信号分选算法
% 输入参数：
%   ChannelOutputs_Baseband: 各站各信道的基带数据
%   StationActiveChannels: 各站的活跃信道检测结果
%   ChannelFrequencies: 信道中心频率
%   Fs: 采样率
%   NumStations: 观测站数量
% 输出参数：
%   source_clusters: 分选结果，每个辐射源的信道集合
%   source_features: 各信道的特征向量

    fprintf('开始多辐射源信号分选...\n');
    
    % 获取所有活跃信道
    AllActiveChannels = [];
    for l = 1:NumStations
        AllActiveChannels = union(AllActiveChannels, find(StationActiveChannels(l,:)));
    end
    
    fprintf('总活跃信道数: %d\n', length(AllActiveChannels));
    
    % 提取每个活跃信道的特征向量
    NumFeatures = 8; % 特征维度
    ChannelFeatures = zeros(length(AllActiveChannels), NumFeatures);
    
    for idx = 1:length(AllActiveChannels)
        ch = AllActiveChannels(idx);
        
        % 特征1-2: 平均功率和功率方差（各站统计）
        powers = zeros(1, NumStations);
        for l = 1:NumStations
            if StationActiveChannels(l, ch)
                signal_data = ChannelOutputs_Baseband{l, ch};
                powers(l) = mean(abs(signal_data).^2);
            else
                powers(l) = 0;
            end
        end
        ChannelFeatures(idx, 1) = mean(powers(powers > 0)); % 平均功率
        ChannelFeatures(idx, 2) = std(powers(powers > 0));  % 功率方差
        
        % 特征3-4: 信号带宽和中心频率稳定性
        bandwidths = zeros(1, NumStations);
        center_freqs = zeros(1, NumStations);
        valid_stations = 0;
        
        for l = 1:NumStations
            if StationActiveChannels(l, ch)
                signal_data = ChannelOutputs_Baseband{l, ch};
                [bw, cf] = estimate_signal_bandwidth(signal_data, Fs/size(StationActiveChannels,2));
                bandwidths(l) = bw;
                center_freqs(l) = cf;
                valid_stations = valid_stations + 1;
            end
        end
        
        ChannelFeatures(idx, 3) = mean(bandwidths(1:valid_stations)); % 平均带宽
        ChannelFeatures(idx, 4) = std(center_freqs(1:valid_stations)); % 频率稳定性
        
        % 特征5-6: 调制特征（循环平稳特征）
        cyclostationary_features = zeros(1, NumStations);
        for l = 1:NumStations
            if StationActiveChannels(l, ch)
                signal_data = ChannelOutputs_Baseband{l, ch};
                cyclostationary_features(l) = extract_cyclostationary(signal_data, Fs/size(StationActiveChannels,2));
            end
        end
        ChannelFeatures(idx, 5) = mean(cyclostationary_features);
        
        % 特征7: 时间到达模式（各站相对时延）  
          time_delays = estimate_relative_delays(ChannelOutputs_Baseband, ch, NumStations, Fs/size(StationActiveChannels,2));
          ChannelFeatures(idx, 6) = std(time_delays); % 时延一致性
        
        % 特征8: 信号空间特征
        spatial_features = extract_spatial_features(ChannelOutputs_Baseband, ch, NumStations);
        ChannelFeatures(idx, 7:8) = spatial_features(1:2);
    end
    
%     % 特征归一化
%     for i = 1:NumFeatures
%         ChannelFeatures(:, i) = (ChannelFeatures(:, i) - mean(ChannelFeatures(:, i))) / std(ChannelFeatures(:, i));
%     end
    
    [source_clusters, ~] = kmeans(ChannelFeatures(:,6), 3);

%     % 使用DBSCAN聚类进行信号分选
%     epsilon = 1.5; % 邻域半径
%     min_pts = 2;   % 最小点数
%     source_clusters = dbscan_clustering(ChannelFeatures(:,3), epsilon, min_pts, AllActiveChannels);
    
    % 输出分选结果
    fprintf('\n=== 信号分选结果 ===\n');
    num_sources = max(source_clusters);
    for i = 1:num_sources
        source_channels = AllActiveChannels(source_clusters == i);
        fprintf('辐射源%d: 信道 %s\n', i, mat2str(source_channels));
        fprintf('       对应频率: ');
        for ch = source_channels
            fprintf('%.3fMHz ', ChannelFrequencies(ch)/1e6);
        end
        fprintf('\n');
    end
    
    % 噪声/异常信道
    noise_channels = AllActiveChannels(source_clusters == -1);
    if ~isempty(noise_channels)
        fprintf('噪声/异常信道: %s\n', mat2str(noise_channels));
    end
    
    source_features.ChannelFeatures = ChannelFeatures;
    source_features.AllActiveChannels = AllActiveChannels;
end

function [bw, center_freq] = estimate_signal_bandwidth(signal, fs_decimated)
% 估计信号带宽和中心频率
    upsampling_factor = 20; % 上采样倍数
    fs_upsampled = fs_decimated * upsampling_factor; % 上采样后的采样率
    % 上采样参考信号
    data_upsampled = upsample(signal, upsampling_factor);
    
    % 设计低通滤波器（抗镜像滤波）
    cutoff_freq = fs_decimated / 2; % 原奈奎斯特频率
    lpf = designfilt('lowpassfir', ...
        'FilterOrder', 64, ...
        'CutoffFrequency', cutoff_freq, ...
        'SampleRate', fs_upsampled);
    
    % 滤波
    data_upsampled = filter(lpf, data_upsampled);



    [psd, freq] = pwelch(data_upsampled, hamming(256), 128, 1024, fs_upsampled, 'centered');
    psd = psd / max(psd);
%     figure
%     plot(freq,psd)
    % 找到峰值
    [~, peak_idx] = max(psd);
    center_freq = freq(peak_idx);
    
    % 计算-3dB带宽
    peak_power = psd(peak_idx);
    threshold = peak_power / 2; % -3dB
    
    % 向左找边界
    left_idx = peak_idx;
    while left_idx > 1 && psd(left_idx) > threshold
        left_idx = left_idx - 1;
    end
    
    % 向右找边界
    right_idx = peak_idx;
    while right_idx < length(psd) && psd(right_idx) > threshold
        right_idx = right_idx + 1;
    end
    
    bw = freq(right_idx) - freq(left_idx);
end

function cyc_feature = extract_cyclostationary(signal, fs)
% 提取循环平稳特征（简化版本）
    % 计算信号的循环自相关
    N = length(signal);
    max_lag = min(100, floor(N/4));
    
    % 计算不同循环频率的相关性
    alpha_values = [0, fs/N, 2*fs/N];
    cyc_values = zeros(1, length(alpha_values));
    
    for i = 1:length(alpha_values)
        alpha = alpha_values(i);
        % 简化的循环自相关计算
        R_alpha = mean(signal .* conj(circshift(signal, floor(alpha*N/fs))) .* ...
                      exp(-1j*2*pi*alpha*(0:N-1)/fs));
        cyc_values(i) = abs(R_alpha);
    end
    
    cyc_feature = max(cyc_values);
end

function time_delays = estimate_relative_delays(ChannelOutputs_Baseband, ch, NumStations, fs_decimated)
% 估计各站相对时延（上采样版本）
% 输入参数：
%   ChannelOutputs_Baseband: 信道化后的基带信号
%   ch: 当前信道编号
%   NumStations: 观测站数量
%   fs_decimated: 降采样后的采样率

    upsampling_factor = 20; % 上采样倍数
    fs_upsampled = fs_decimated * upsampling_factor; % 上采样后的采样率
    
    time_delays = zeros(1, NumStations);
    
    % 获取参考站（站1）的信号
    ref_data_original = ChannelOutputs_Baseband{1, ch};
    
    % 上采样参考信号
    ref_data_upsampled = upsample(ref_data_original, upsampling_factor);
    
    % 设计低通滤波器（抗镜像滤波）
    cutoff_freq = fs_decimated / 2; % 原奈奎斯特频率
    lpf = designfilt('lowpassfir', ...
        'FilterOrder', 64, ...
        'CutoffFrequency', cutoff_freq, ...
        'SampleRate', fs_upsampled);
    
    % 滤波
    ref_data_upsampled = filter(lpf, ref_data_upsampled);
    
    for l = 1:NumStations
        if l == 1
            time_delays(l) = 0; % 参考站时延为0
        else
            % 获取当前站的信号
            test_data_original = ChannelOutputs_Baseband{l, ch};
            
            if isempty(test_data_original) || length(test_data_original) < 10
                time_delays(l) = 0;
                continue;
            end
            
            % 上采样测试信号
            test_data_upsampled = upsample(test_data_original, upsampling_factor);
            test_data_upsampled = filter(lpf, test_data_upsampled);
            
            % 方法1：广义互相关（GCC-PHAT） - 对多径有较好鲁棒性
%             delay_samples = gcc_phat_delay_estimation(ref_data_upsampled, test_data_upsampled, fs_upsampled,ch);
            
            % 方法2：传统互相关（备选）
            delay_samples = xcorr_delay_estimation(ref_data_upsampled, test_data_upsampled, fs_upsampled,ch);
            
            time_delays(l) = delay_samples;
        end
    end
    
    fprintf('信道%d时延估计: ', ch);
    for l = 1:NumStations
        fprintf('站%d: %.3fμs ', l, time_delays(l)*1e6);
    end
    fprintf('\n');
end

function delay = gcc_phat_delay_estimation(ref_signal, test_signal, fs,ch)
% 广义互相关-相位变换（GCC-PHAT）时延估计
% 对多径环境和低信噪比有较好鲁棒性

    % 确保信号长度一致
    min_len = min(length(ref_signal), length(test_signal));
    ref_signal = ref_signal(1:min_len);
    test_signal = test_signal(1:min_len);
    
    % 计算互功率谱密度
    NFFT = 2^nextpow2(2*min_len-1);
    
    R1 = fft(ref_signal, NFFT);
    R2 = fft(test_signal, NFFT);
    
    % 计算互功率谱
    G12 = R1 .* conj(R2);
    
    % PHAT加权
    G12_phat = G12 ./ (abs(G12) + eps); % 加eps避免除零
    
    % 逆FFT得到广义互相关函数
    gcc = ifft(G12_phat);
    gcc = fftshift(gcc);
    
    % 生成时延轴
    lags = (-NFFT/2 : NFFT/2-1) / fs;
    
    % 找到峰值位置
    [~, max_idx] = max(abs(gcc));
    delay = lags(max_idx);
    
    if ch == 1
        figure
        plot(lags*fs,abs(gcc))
    end
    % 二次插值提高精度
    if max_idx > 1 && max_idx < length(gcc)
        y1 = abs(gcc(max_idx-1));
        y2 = abs(gcc(max_idx));
        y3 = abs(gcc(max_idx+1));
        
        % 抛物线插值
        delta = (y3 - y1) / (2*(2*y2 - y1 - y3));
        delay = delay + delta/fs;
    end
end

function delay = xcorr_delay_estimation(ref_signal, test_signal, fs,ch)
% 传统互相关时延估计（备选方法）

    % 确保信号长度一致
    min_len = min(length(ref_signal), length(test_signal));
    ref_signal = ref_signal(1:min_len);
    test_signal = test_signal(1:min_len);
    
    % 计算互相关
    [corr_seq, lags] = xcorr(test_signal, ref_signal,'biased');
    
    % 找到峰值
    [~, max_idx] = max(abs(corr_seq));
    delay_samples = lags(max_idx);
    
    % 转换为时间
    delay = delay_samples / fs;
    
    if ch == 4
        figure
        plot(lags,abs(corr_seq))
    end



    % 二次插值提高精度
    if max_idx > 1 && max_idx < length(corr_seq)
        y1 = abs(corr_seq(max_idx-1));
        y2 = abs(corr_seq(max_idx));
        y3 = abs(corr_seq(max_idx+1));
        
        delta = (y3 - y1) / (2*(2*y2 - y1 - y3));
        delay = delay + delta/fs;
    end
end

function delay = mle_delay_estimation(ref_signal, test_signal, fs)
% 最大似然时延估计（高精度版本）
% 适用于高信噪比场景

    min_len = min(length(ref_signal), length(test_signal));
    ref_signal = ref_signal(1:min_len);
    test_signal = test_signal(1:min_len);
    
    % 粗估计（使用GCC-PHAT）
    delay_coarse = gcc_phat_delay_estimation(ref_signal, test_signal, fs);
    
    % 精细搜索
    search_range = 2/fs; % ±2个采样点
    num_points = 41;
    tau_grid = linspace(delay_coarse - search_range, delay_coarse + search_range, num_points);
    
    likelihood = zeros(size(tau_grid));
    
    for i = 1:length(tau_grid)
        tau = tau_grid(i);
        
        % 时延补偿（分数时延插值）
        delayed_ref = fractional_delay(ref_signal, tau, fs);
        
        % 计算似然函数（相关系数）
        correlation = corrcoef(delayed_ref(100:end-100), test_signal(100:end-100));
        likelihood(i) = abs(correlation(1,2));
    end
    
    % 找到最大似然点
    [~, max_idx] = max(likelihood);
    delay = tau_grid(max_idx);
    
    % 抛物线插值
    if max_idx > 1 && max_idx < length(tau_grid)
        y1 = likelihood(max_idx-1);
        y2 = likelihood(max_idx);
        y3 = likelihood(max_idx+1);
        
        delta = (y3 - y1) / (2*(2*y2 - y1 - y3));
        delay = delay + delta * (tau_grid(2)-tau_grid(1));
    end
end

function delayed_signal = fractional_delay(signal, delay, fs)
% 分数时延补偿（频域实现）
    N = length(signal);
    f = (0:N-1) * fs / N;
    f(f > fs/2) = f(f > fs/2) - fs;
    
    % 频域相位旋转
    phase_shift = exp(-1j * 2 * pi * f * delay);
    
    % 应用相位旋转
    spectrum = fft(signal);
    delayed_spectrum = spectrum .* phase_shift;
    delayed_signal = real(ifft(delayed_spectrum));
end

function spatial_features = extract_spatial_features(ChannelOutputs_Baseband, ch, NumStations)
% 提取空间特征
    spatial_features = zeros(1, 3);
    
    % 特征1: 信号在各站的一致性
    powers = zeros(1, NumStations);
    for l = 1:NumStations
        signal_data = ChannelOutputs_Baseband{l, ch};
        powers(l) = mean(abs(signal_data).^2);
    end
    spatial_features(1) = std(powers) / mean(powers); % 功率变化系数
    
    % 特征2: 信号的空间分布均匀性
    spatial_features(2) = entropy(powers / sum(powers));
end

function clusters = dbscan_clustering(features, epsilon, min_pts, channel_indices)
% DBSCAN聚类算法
    n = size(features, 1);
    clusters = zeros(n, 1);
    visited = false(n, 1);
    cluster_id = 0;
    
    for i = 1:n
        if ~visited(i)
            visited(i) = true;
            
            % 找到epsilon邻域内的点
            neighbors = region_query(features, i, epsilon);
            
            if length(neighbors) < min_pts
                % 标记为噪声
                clusters(i) = -1;
            else
                cluster_id = cluster_id + 1;
                clusters = expand_cluster(features, i, neighbors, cluster_id, ...
                                        epsilon, min_pts, clusters, visited);
            end
        end
    end
    
    % 输出聚类信息
    fprintf('聚类完成: 发现 %d 个辐射源\n', cluster_id);
    for i = 1:cluster_id
        cluster_channels = channel_indices(clusters == i);
        fprintf('  聚类%d: %d个信道\n', i, length(cluster_channels));
    end
end

function neighbors = region_query(features, point_idx, epsilon)
% 查询epsilon邻域内的点
    n = size(features, 1);
    distances = zeros(n, 1);
    
    for i = 1:n
        distances(i) = norm(features(point_idx, :) - features(i, :));
    end
    
    neighbors = find(distances <= epsilon);
end

function clusters = expand_cluster(features, point_idx, neighbors, cluster_id, ...
                                 epsilon, min_pts, clusters, visited)
% 扩展聚类
    clusters(point_idx) = cluster_id;
    
    i = 1;
    while i <= length(neighbors)
        neighbor_idx = neighbors(i);
        
        if ~visited(neighbor_idx)
            visited(neighbor_idx) = true;
            
            new_neighbors = region_query(features, neighbor_idx, epsilon);
            
            if length(new_neighbors) >= min_pts
                neighbors = union(neighbors, new_neighbors);
            end
        end
        
        if clusters(neighbor_idx) == 0
            clusters(neighbor_idx) = cluster_id;
        end
        
        i = i + 1;
    end
end








