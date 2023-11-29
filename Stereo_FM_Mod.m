function [num_of_samples] = Stereo_FM_Mod(audio1, Fs1, audio2, Fs2, verbose)
% Function that modualtes two audios: audio1 and audio2 with Fs1 and Fs2 respectively
% generates one audio which is combined from the two audios.
% verbose: if 1 then the function plots figures

% Setting the audios to the same length if needed
if (length(audio1) ~= length(audio2))
    n = min(length(audio1), length(audio2));
    audio1 = audio1(1:n);
    audio2 = audio2(1:n);
end


% Upsamling to 200khz
F_up = 200e3;
y_upsampled_1 = resample(audio1, F_up, Fs1); 
y_upsampled_2 = resample(audio2, F_up, Fs2); 

 
% LPF filter with cutoff frequency 95kHz
Fc1 = 15e3; 
LPF = design(fdesign.lowpass('N,Fc', 200, Fc1, F_up));
y1_filtered = filter(LPF, y_upsampled_1);
y2_filtered = filter(LPF, y_upsampled_2); 


% Setting the audios to the same length if needed
if (length(y1_filtered) ~= length(y2_filtered))
    n = min(length(y1_filtered), length(y2_filtered));
    y1_filtered = y1_filtered(1:n);
    y2_filtered = y2_filtered(1:n);
end

% DSB-SC modulation
F_pb = 19e3; % Pilot
t = ((1:length(y1_filtered))/F_up)';
DSB_SC = (0.9*((y1_filtered+y2_filtered)/2 +...
        (y1_filtered-y2_filtered)/2.*sin(4*pi*F_pb.*t))...
        + 0.1*sin(2*pi*F_pb.*t))*75000; % Formula from the reading material


% FM modulation with delta_f = 75kHz
delta_f = 75e3;
FM_BB = exp(1j*2*pi*((delta_f/max(abs(DSB_SC)).*cumsum(DSB_SC)./F_up)));

 
% Upsampling the BB signal to 1MHz 
F_up2 = 1e6;
FM_BB_upsampled = resample(FM_BB, F_up2, F_up);

% LPF filter with cutoff frequency 128kHz
Fc2 = 128e3;
LPF = design(fdesign.lowpass('N,Fc',200,Fc2,F_up2));
FM_BB_filtered = filter(LPF, FM_BB_upsampled);

% Saving the data to bin file
WR_bin_file('data_tx_stereo.bin',FM_BB_filtered);
num_of_samples = length(FM_BB_filtered);

if(verbose)
    l1 = length (audio1);
    l2 = length (audio2);

    figure
    subplot(2,2,1);
    plot((1:l1)/ Fs1,audio1)
    title('Audio1 signal in time domain');
    xlabel('t [sec]');
    ylabel('s_1(t)');

    subplot(2,2,3);
    plot(((1:l1)-l1/2)/l1* Fs1/1000,10*log10(abs(fftshift(fft(audio1)))));
    title('Spectrum of audio1');
    xlabel('f [kHz]');
    ylabel('S_1(f) [dB]');
    
    subplot(2,2,2);
    plot((1:l2)/ Fs2,audio2)
    title('Audio2 signal in time domain');
    xlabel('t [sec]');
    ylabel('s_2(t)');

    subplot(2,2,4);
    plot(((1:l2)-l2/2)/l2* Fs2/1000,10*log10(abs(fftshift(fft(audio2)))));
    title('Spectrum of audio2');
    xlabel('f [kHz]');
    ylabel('S_2(f) [dB]');
    
    figure
    subplot(2,2,1);
    l1 = length (y1_filtered);
    plot((1:l1)/ F_up,y1_filtered)
    title('Audio1 filtered signal in time domain');
    xlabel('t [sec]');
    ylabel('s_1(t)');

    subplot(2,2,3);
    l1 = length (y1_filtered);
    plot(((1:l1)-l1/2)/l1*F_up/1000,10*log10(abs(fftshift(fft(y1_filtered)))));
    title('Audio1 - spectrum of the filtered signal');
    xlabel('f [kHz]');
    ylabel('S_1(f) [dB]');
    
    subplot(2,2,2);
    l2 = length (y2_filtered);
    plot((1:l2)/ F_up,y2_filtered)
    title('Audio2 filtered signal in time domain');
    xlabel('t [sec]');
    ylabel('s_2(t)');

    subplot(2,2,4);
    l2 = length (y2_filtered);
    plot(((1:l2)-l2/2)/l2*F_up/1000,10*log10(abs(fftshift(fft(y2_filtered)))));
    xlabel('f [kHz]');ylabel('S_2(f) [dB]');
    title('Audio2 - spectrum of the filtered signal');
    
    figure
    subplot(2,2,1);
    l2 = length (DSB_SC);
    plot((1:l2)/ F_up,DSB_SC)
    title('Generated audio in time domain');
    xlabel('t [sec]');
    ylabel('s(t)');

    subplot(2,2,3);
    plot(((1:l1)-l1/2)/l1*F_up/1000,10*log10(abs(fftshift(fft(DSB_SC)))));
    title('Spectrum of the generated audio');
    xlabel('f [kHz]');
    ylabel('S(f) [dB]');
    
    subplot(2,2,2);
    l2 = length (FM_BB);
    plot((1:l2)/ F_up,FM_BB)
    title('FM modulated signal in time domain');
    xlabel('t [sec]');
    ylabel('s_{FM}(t)');

    subplot(2,2,4)
    FM_BB_fft = fft(FM_BB);
    l = length (FM_BB_fft);
    plot(((1:l)-l/2)/l*F_up/1000,10*log10(abs(fftshift(FM_BB_fft))))
    title('Spectrum of the FM modulated signal');
    xlabel('f [kHz]');
    ylabel('S_{FM}(f) [dB]');

    figure
    subplot(2,1,1)
    l2 = length (FM_BB_filtered);
    plot((1:l2)/ Fc2,FM_BB_filtered)
    title('FM modulated filtered signal in time domain');
    xlabel('t [sec]');
    ylabel('s(t)');

    subplot(2,1,2)
    FM_BB_filtered_fft = fft(FM_BB_filtered);
    l = length (FM_BB_filtered_fft);
    plot(((1:l)-l/2)/l*F_up/1000,10*log10(abs(fftshift(FM_BB_filtered_fft))))
    title('Spectrum of the FM modulated filtered signal');
    xlabel('f [kHz]');
    ylabel('S_{FM}(f) [dB]');

% end
end

