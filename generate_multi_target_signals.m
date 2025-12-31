function RxSignals_MultiTarget = generate_multi_target_signals(TargetPositions, TargetPowers, HopPatterns, StationPos, Fs, T_hop, NumHops,SymbolRates)
    % 生成多目标接收信号
    % 输入：
    %   TargetPositions - 目标位置矩阵 [x1,y1; x2,y2; ...]
    %   TargetPowers - 每个目标的信号功率 [p1, p2, ...]
    %   HopPatterns - 跳频图案细胞数组 {pattern1, pattern2, ...}
    %   StationPos - 观测站位置
    %   Fs - 采样率
    %   T_hop - 每跳持续时间
    %   NumHops - 跳数
    % 输出：
    %   RxSignals_MultiTarget - 多站接收信号矩阵
    
    NumTargets = size(TargetPositions, 1);
    NumStations = size(StationPos, 1);
    TotalSamples = NumHops * T_hop * Fs;
    
    % 基础参数（使用您之前定义的）
    Modulation = 'QPSK';
    FixedNoisePower = 1;  % 固定噪声功率为1（可根据需要调整）
    IF_CenterFreqs_MHz = [48.16, 48.48, 48.80, 49.12, 49.44, 49.76, 50.08, 50.40, ...
                          50.72, 51.04, 51.36, 51.68, 52.00, 52.32, 52.64, 52.96, ...
                          53.28, 53.60, 53.92, 54.24];
    IF_Center = 51.20e6;
    Baseband_Freqs = (IF_CenterFreqs_MHz * 1e6) - IF_Center;
    
    fprintf('生成 %d 个目标的跳频信号...\n', NumTargets);
    
    % 初始化接收信号
    RxSignals_MultiTarget = zeros(NumStations, TotalSamples);
    
    % 为每个目标生成信号
    for target_id = 1:NumTargets
        fprintf('  生成目标%d的信号...\n', target_id);
        SymbolRate = SymbolRates(target_id);
        
        % 生成单个目标的跳频信号
        target_signal = generate_single_target_signal(...
            HopPatterns{target_id}, Baseband_Freqs, ...
            SymbolRate, Fs, T_hop, NumHops, TargetPowers(target_id));
        
        % 计算该目标到各站的时延
        target_pos = TargetPositions(target_id, :);
        TrueTDOA = zeros(NumStations, 1);
        c = 3e8;
        
        for l = 1:NumStations
            TrueTDOA(l) = (norm(StationPos(l,:) - target_pos) - norm(StationPos(1,:) - target_pos)) / c;
        end
        
        % 模拟信号传播到各观测站
        for l = 1:NumStations
            % 应用时延
            DelaySamples = round(TrueTDOA(l) * Fs);
            DelayedSignal = circshift(target_signal, DelaySamples);
            
            % 添加噪声（使用固定SNR）
            SNR_dB = 20;

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
            
            if target_id ==  NumTargets
               RxSignals_MultiTarget(l, :) = RxSignals_MultiTarget(l, :) + ScaledSignal + Noise;
            else
               RxSignals_MultiTarget(l, :) = RxSignals_MultiTarget(l, :) + ScaledSignal;
            end
            
        
        end
    end
    
    fprintf('多目标信号生成完成\n');
end

function target_signal = generate_single_target_signal(hop_pattern, baseband_freqs, symbol_rate, Fs, T_hop, num_hops, power_scale)
    % 生成单个目标的跳频信号
    
    SamplesPerSymbol = floor(Fs / symbol_rate);
    SymbolsPerHop = floor(T_hop * symbol_rate);
    TotalSamples = num_hops * T_hop * Fs;
    SamplesPerHop = round(T_hop * Fs);


    TimeVector = (0:TotalSamples-1) / Fs;
    target_signal = zeros(1, TotalSamples);
    
    for h = 1:num_hops
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
        current_channel = hop_pattern(mod(h-1, length(hop_pattern)) + 1);
        CurrentFreq = baseband_freqs(current_channel);
        Carrier = exp(1j * 2 * pi * CurrentFreq * TimeVector(HopSamples));
        
        target_signal(HopSamples) = SymbolWaveform .* Carrier * sqrt(power_scale);
    end
end