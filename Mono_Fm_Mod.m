function [num_of_samples] = Mono_Fm_Mod(y, Fs, verbose)
% Function that applies frequency modulation on Mono signal, and generated
% bin file with the modulated signal
% if verbose is 1 then the function plots figures

% Upsamling to 200khz
F_up = 200e3;
y_upsampled = resample(y, F_up, Fs);  

% LPF filter with cutoff frequency 15kHz - Mono from 0-15kHz
Fc = 15e3; 
LPF = design(fdesign.lowpass('N,fc', 200, Fc, F_up));
y_filtered = filter(LPF, y_upsampled);

% FM modulation with delta_f = 75kHz
delta_f = 75e3;
FM_BB = exp(1j*2*pi*((delta_f/max(abs(y_filtered)).*cumsum(y_filtered)./F_up)));

% Upsampling the BB signal to 1MHz 
F_up2 = 1e6;
FM_BB_upsampled = resample(FM_BB, F_up2, F_up);  % Up samling to 200khz

% LPF filter with cutoff frequency 90kHz
Fc2 = 90e3; % According to Carsons rule Fc2 = fm + delta_f
LPF = design(fdesign.lowpass('N,fc', 200, Fc2, F_up2));
FM_BB_filtered = filter(LPF, FM_BB_upsampled);

% Saving the data to bin file
WR_bin_file('data_tx_Mono.bin',FM_BB_filtered);

if (verbose)
    figure
    subplot(2,1,1)
    l = length(y);
    plot((1:l)/Fs,y)
    title('Audio signal in time domain');
    xlabel('t [sec]');
    ylabel('s(t)');

    subplot(2,1,2)
    plot(((1:l)-l/2)/l*(Fs/1000),10*log10(abs(fftshift(fft(y)))));
    title('Spectrum of the signal');
    xlabel('f [kHz]');
    ylabel('S(f) [dB]');
    
    figure
    subplot(2,1,1)
    l = length(y_filtered);
    plot((1:l)/Fs,y_filtered);
    title('Filtered audio signal');
    xlabel('t [sec]');
    ylabel('s(t)');

    subplot(2,1,2)
    l = length(y_filtered);
    plot(((1:l)-l/2)/l*F_up/1000,10*log10(abs(fftshift(fft(y_filtered)))));
    title('Spectrum of the signal');
    xlabel('f [kHz]');
    ylabel('S(f) [dB]');

    figure
    subplot(2,1,1)
    l = length(FM_BB);
    plot((1:l)/Fs,FM_BB);
    title('Audio FM signal');
    xlabel('t [sec]');
    ylabel('s_{FM}(t)');

    subplot(2,1,2)
    l = length(FM_BB);
    plot(((1:l)-l/2)/l*F_up/1000,10*log10(abs(fftshift(fft(FM_BB)))))
    title('Audio FM signal in frequency domain');
    xlabel('f [kHz]');
    ylabel('S_{FM}(f) [dB]');

    figure
    subplot(2,1,1)
    l = length(FM_BB_filtered);
    plot((1:l)/Fs,FM_BB_filtered);
    title('Filtered audio FM signal');
    xlabel('t [sec]');
    ylabel('s(t)');

    subplot(2,1,2)
    l = length(FM_BB_filtered);
    plot(((1:l)-l/2)/l*F_up/1000,10*log10(abs(fftshift(fft(FM_BB_filtered)))));
    title('Spectrum of the signal');
    xlabel('f [kHz]');
    ylabel('S(f) [dB]');
  
end
num_of_samples = length(FM_BB_filtered);
end

% t = (1:length(y_filtered))/F_up;
% FM_BB = exp(1j*2*pi*(delta_f*t + (Fc*cumsum(y_filtered)./F_up)'));