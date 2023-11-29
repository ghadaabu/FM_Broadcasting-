function [] = Stereo_FM_DeMod(num_of_samples, verbose, play)
% Function that demodulates our generated bin file which contains data of
% an audio. Gets two parameters:
% - verbose: the function plots figures if 1
% - play: the function plays the audio if 1

data = RD_bin_file('data_rx.bin',num_of_samples);
Fs = 1e6;

% Discarding the first 10000 of the initial samples for a clear audio reception
data(1:10e3) = 0;

% LPF filter with cutoff frequency 128kHz
Fc1 = 128e3;
LPF = design(fdesign.lowpass('N,fc',200,Fc1,Fs));
data_filtered = filter(LPF, data);

% Downsampling to 200kHz
F_dn1 = 200e3;
data_dws = resample(data_filtered,F_dn1,Fs);

% Demodulation
delta_f = 75e3;
data_demodulated = diff(unwrap(angle(data_dws)))*F_dn1./(2*pi*delta_f);
% We get the angle of each sample, then unwarp removes phase jumps greater
% than pi, then we apply diff to get phase changes

% LPF filter with cutoff frequency 15kHz
Fc2 = 15e3;
LPF = design(fdesign.lowpass('N,Fc',500,Fc2,F_dn1));
Mono_audio = filter(LPF, data_demodulated);

fp = 19e3;
t = ((1:length(data_demodulated))/F_dn1);
I = data_demodulated.*cos(2*pi*fp*t);  
Q = data_demodulated.*sin(2*pi*fp*t); 

I = filter(LPF, I);
Q = filter(LPF, Q);

phase_err = atan(Q./I); 
stereo_audio = 2/0.9*data_demodulated.*sin(4*pi*fp*t + phase_err); % carrier frequency is 2*fp = 38kHz
stereo_audio = filter(LPF, stereo_audio);
stereo_audio(isnan(stereo_audio)) = 0;

%stereo_audio=(L-R)
%Mono_audio=(L+R)
left = (Mono_audio + stereo_audio)/2; %((L+R)+(L-R))/2=L
right = (Mono_audio - stereo_audio)/2;%((L+R)-(L-R))/2=R

% Downsampling to 44.1kHz
F_dn2 = 44.1e3;
left = resample(left,F_dn2,F_dn1);
right = resample(right,F_dn2,F_dn1);

% Composing the audio such that audio1 is heard on the right ear and
% audio2 is heard on the left ear
audio(:,1) = right;
audio(:,2) = left;
if (play)
    soundsc(audio(:,1),(44.1)*1e3);
%     soundsc(audio,(48)*1e3);

end

if (verbose)
    figure
    subplot(3,1,1)
    l = length (data);
    plot(((1:l)-l/2)/l*Fs/1000,10*log10(abs(fftshift(fft(data)))));
    title('Spectrum of the recieved signal');
    xlabel('f [kHz]');
    ylabel('S_{received}(f) [dB]');
    
    subplot(3,1,2)
    DATA_filt_res = fft(data_dws);
    l = length (DATA_filt_res);
    plot(((1:l)-l/2)/l*F_dn1/1000,10*log10(abs(fftshift(DATA_filt_res))))
    title('Spectrum of the filtered and downsampled signal');
    xlabel('f [kHz]');
    ylabel('S_{filtered}(f) [dB]');

    subplot(3,1,3)
    demod_DATA = fft(data_demodulated);
    l = length (demod_DATA);
    plot(((1:l)-l/2)/l*F_dn1/1000,10*log10(abs(fftshift(demod_DATA))))
    title('Spectrum of the demodulated signal');
    xlabel('f [kHz]');
    ylabel('S_{demod}(f) [dB]');

    figure
    subplot(2,1,1)
    plot(((1:l)-l/2)/l*F_dn1/1000,10*log10(abs(fftshift(fft(Mono_audio)))));
    title('Mono spectrum');
    xlabel('f [kHz]');
    ylabel('S_{Mono}(f) [dB]');

    subplot(2,1,2)
    plot(((1:l)-l/2)/l*F_dn1/1000,10*log10(abs(fftshift(fft(stereo_audio)))));
    title('Stereo spectrum');
    xlabel('f [kHz]');
    ylabel('S_{Stereo}(f) [dB]');

    figure
    subplot(2,2,1)
    Left = fft(left);
    l=length(Left);
    plot((1:l)/F_dn1,left)
    title('Signal of the left ear audio');
    xlabel('t [sec]');
    ylabel('s_{Left}(t)');

    subplot(2,2,3);
    plot(((1:l)-l/2)/l*F_dn1/1000,10*log10(abs(fftshift(Left))));
    title('Spectrum of the left ear audio');
    xlabel('f [kHz]');
    ylabel('S_{Left}(f) [dB]');

    subplot(2,2,2);
    plot((1:l)/F_dn1,right)
    title('Signal of the right ear audio');
    xlabel('t [sec]');
    ylabel('s_{Right}(t)');

    subplot(2,2,4);
    Right = fft(right);
    plot(((1:l)-l/2)/l*F_dn1/1000,10*log10(abs(fftshift(Right))));
    title('Spectrum of the right ear audio');
    xlabel('f [kHz]');
    ylabel('S_{Right}(f) [dB]');
end
end

