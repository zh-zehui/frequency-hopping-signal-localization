function [record_data, record_time]= extract_data_2(filePath,SampleRate,data_length,LowPassFilterStop_Freq,freq_offset,MonteCarlo,index,IQ_exchange,index_data)
% index  1  signal 做低通滤波
%        2  noise  不做低通滤波


% IQ_exchange  1 不对调
%              >1  对调
   

File = dir(fullfile(filePath));
File = File(~ismember({File.name}, {'.', '..'}));

record_data = [];
record_time = [];

   for i = index_data

        filename = [filePath,char(File(i).name)];
        % 从文件中读取内容
        fileContent = fileread(filename);
    
        % 使用正则表达式匹配数据
        pattern = '(\d+日) (\d+:\d+:\d+) (\d+) (\/) ([\d.])+M|(-?\d+),(-?\d+)';
        matches = regexp(fileContent, pattern, 'tokens');

        % 如果找到匹配的数据
        if ~isempty(matches)


            match = matches{1};
            record_time = [record_time,str2num(match{3})/204.8e6];
            
            signal = readmatrix(filename);
            
            if  IQ_exchange>1
                signal = signal(1:end,2) + 1j*signal(1:end,1);   %中频输入接收到IQ数据
            else
                signal = signal(1:end,1)+1j*signal(1:end,2);   %中频输入接收到IQ数据
            end
%             signal = (signal -mean(signal))/std(signal);

            signal = signal -mean(signal);
           
            if index == 1
               signal = Low_fileter(signal,LowPassFilterStop_Freq/2,SampleRate);
            end
%            %% 频率偏移校正
%              freq_offset = 2.05312e3; % 信号发生器发;
%            freq_offset  = 1.90933e3 ;     %靶标发1
%              freq_offset = 2.048e3;  % 假设有1kHz的频率偏移，需要根据实际情况调整          
             t_decimated = (0:length(signal)-1) / SampleRate;
             signal = signal .* exp(-1i*2*pi*freq_offset*t_decimated.');
             

             record_data = [record_data, single(signal(1:data_length))];
             
         else
             disp('未找到匹配的数据');
        end
    end

end