[Piano_audio, Fs_Piano] = audioread('piano2.wav');
Piano_audio = Piano_audio(:,1) + Piano_audio(:,2);
% soundsc(Piano_audio,Fs_Piano)

[music_audio, Fs_music] = audioread('music1.wav');
music_audio = music_audio(:,1) + music_audio(:,2);
% soundsc(music_audio,Fs_music)

%% PART 1 - Mono Audio Stream
% num_of_samples = Mono_Fm_Mod(Piano_audio, Fs_Piano, 0); % Mono_Fm_Mod(Piano_audio, Fs_Piano, verbose) verbose is a flag that enables the function to prints figures
% Mono_Fm_DeMod('data_rx.bin', num_of_samples, 1, 0) % Mono_Fm_DeMod(file_name, num_of_samples, verbose, play)
% % Mono_Fm_DeMod('rx_data_4MSamples1MsPs_978MHz.bin', 4e6, 0, 1) 
% % 
% % %% PART 2 - Stereo Audio Stream
num_of_samples = Stereo_FM_Mod(Piano_audio, Fs_Piano, music_audio, Fs_music, 1); %
% Stereo_FM_DeMod(num_of_samples, 1, 0)