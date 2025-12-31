function features = extract_pulse_features(PulseList, ChannelOutputs, StationPos)
    % 提取每个脉冲的多维特征
    
    NumPulses = length(PulseList);
    features = zeros(NumPulses, 6);  % [时间, 频率, 幅度, DOA, 带宽, 稳定性]
    
    for i = 1:NumPulses
        pulse = PulseList(i);
        
        % 基础特征
        features(i, 1) = pulse.timestamp;      % 到达时间
        features(i, 2) = pulse.frequency;      % 载频
        features(i, 3) = pulse.amplitude;      % 幅度
        
        % DOA估计（基于多站相位差）
        features(i, 4) = estimate_doa(pulse, ChannelOutputs, StationPos);
        
        % 信号带宽特征
        features(i, 5) = estimate_bandwidth(pulse, ChannelOutputs);
        
        % 幅度稳定性（多个站的一致性）
        features(i, 6) = amplitude_consistency(pulse);
    end
end

function doa = estimate_doa(pulse, ChannelOutputs, StationPos)
    % 基于多站相位差估计DOA
    
    NumStations = size(ChannelOutputs, 1);
    phases = zeros(1, NumStations);
    
    for l = 1:NumStations
        signal_data = ChannelOutputs{l, pulse.channel};
        % 提取脉冲对应时间段的相位信息
        phase = angle(mean(signal_data(pulse.time_indices)));
        phases(l) = phase;
    end
    
    % 基于相位差和站间几何关系估计DOA
    % 简化版：使用第一个站作为参考
    phase_diffs = phases - phases(1);
    doa = unwrap(phase_diffs) * 180/pi;  % 转换为度
    
    doa = mean(doa(2:end));  % 返回平均DOA估计
end