function [] = Mono_Fm_DeMod(file_name, num_of_samples, verbose, play)
% num_of_samples: 6306500 for our modulated signal file 'data_tx.bin', 
% 4e6 for the provided file 'rx_data_4MSamples1MsPs_978MHz.bin'.
% Gets two additional parameters:
% - verbose: the function plots figures if 1
% - play: the function plays the audio if 1

Fs = 1e6;

% Loading data 
data = RD_bin_file(file_name, num_of_samples);

% Discarding the first 10000 of the initial samples for a clear audio reception 
data(1:10e4) = 0;

% LPF filter with cutoff frequency 90kHz
Fc2 = 90e3;
LPF = design(fdesign.lowpass('N,fc', 200, Fc2, Fs));
data_filtered = filter(LPF, data);


% Downsampling to 200kHz
F_dn = 200e3;
data_dws = resample(data_filtered, F_dn, Fs);

% Demodulation of the FM signal
delta_f = 75e3;
data_demodulated = diff(unwrap(angle(data_dws)))*F_dn./(2*pi*delta_f);
% We get the angle of each sample, then unwarp removes phase jumps greater
% than pi, then we apply diff to get phase changes 

% LPF filter with cutoff frequency 15kHz
Fc = 15e3;
LPF = design(fdesign.lowpass('N,Fc',150,Fc,F_dn));
y_filtered = filter(LPF, data_demodulated);

% Downsampling to the original sampling rate of the piano audio 48kHz
F_dn2 = 48e3;
audio = resample(y_filtered,F_dn2,F_dn);
audio_fft = fft(audio);

if (verbose)

    figure
    subplot(2,1,1)
    l = length(data);
    plot((1:l)/Fs,data)
    title('Recieved signal in time domain');
    xlabel('t [sec]');
    ylabel('y(t)');

    subplot(2,1,2)
    l = length (data);
    plot(((1:l)-l/2)/l*Fs/1000, 10*log10(abs(fftshift(fft(data)))));
    title('Spectrum of the recieved signal');
    xlabel('f [kHz]');
    ylabel('Y(f) [dB]');

    figure
    subplot(2,1,1)
    l = length(data_dws);
    plot((1:l)/F_dn,data_dws)
    title('Filtered and down sapmled recieved signal in time domain');
    xlabel('t [sec]');
    ylabel('y(t)');

    subplot(2,1,2)
    data_dws_fft = fft(data_dws);
    l = length (data_dws_fft);
    plot(((1:l)-l/2)/l*F_dn/1000, 10*log10(abs(fftshift(data_dws_fft))))
    title('Spectrum of the filtered and down sapmled recieved signal');
    xlabel('f [kHz]');
    ylabel('Y(f) [dB]');
    
    figure
    subplot(2,1,1)
    l = length (data_demodulated);
    plot((1:l)/F_dn,data_demodulated)
    title('Demodulated signal in time domain');
    xlabel('t [sec]');
    ylabel('y(t)');

    subplot(2,1,2)
    data_demodulated_fft = fft(data_demodulated);
    l = length (data_demodulated_fft);
    plot(((1:l)-l/2)/l*F_dn/1000, 10*log10(abs(fftshift(data_demodulated_fft))))
    title('Spectrum of the demodulated signal');
    xlabel('f [kHz]');
    ylabel('S_d(f) [dB]');
    
    figure
    subplot(2,1,1)
    l = length (audio_fft);
    plot((1:l)/F_dn2,audio)
    title('Audio signal in time domain');
    xlabel('t [Sec]');
    ylabel('x(t)');

    subplot(2,1,2)
    plot(((1:l)-l/2)/l*F_dn2/1000,10*log10(abs(fftshift(audio_fft))))
    title('Spectrum of the signal');
    xlabel('f [kHz]');
    ylabel('X(f) [dB]');

    if(file_name == 'rx_data_4MSamples1MsPs_978MHz.bin')
        figure
        l = length (data_demodulated_fft);
        plot(((1:l)-l/2)/l*F_dn/1000, 10*log10(abs(fftshift(data_demodulated_fft))))
        title('Spectrum of the demodulated signal');
        xlabel('f [kHz]');
        ylabel('S_d(f) [dB]');
    end
end

if(play)
    soundsc(audio,F_dn2) 
end
end
