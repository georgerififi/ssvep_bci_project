clc; clear all;
close all;
warning off;

%% Load data
% directory
PATH_CUR = pwd;
idcs = strfind(PATH_CUR, '/');
PATH_DATA = [PATH_CUR(1:idcs(end)-1) '/dataset'];
addpath(PATH_DATA);

% define dataset parameters
fs = 250; % sampling rate 
t_start = .5 ; % .5sec
t_end_minus = .5; % .5sec
idx_first = t_start*fs + 1; % 126
idx_last_minus = t_end_minus*fs; % 125
electr = [48, 54, 55, 56, 57, 58, 61, 62, 63]; % electrodes
load('Freq_Phase.mat'); %  loads frequencies of targets in proper order 'freqs' and the phases vector 'phases' also.
fr = freqs;
ph = phases;
n_fr = length(fr); % number of frequencies 
n_ph = length(phases); % number of phases
load('S1.mat'); % load the 1st subject to define the following..
[n_channels, full_time, n_fr, n_sessions] = size(data(electr, idx_first:end-idx_last_minus, :, :));
n_h = 5; % number of harmonics
time_step = fs/4;  % time lengths, 1s = 250 time points, .25sec = 62.5 time points
time_length = ceil(time_step:time_step:full_time); % 63, 125, 188, 250, ..., 1250 
n_tlength = length(time_length); % number of time lengths to try
dt = 1/fs;
TW = time_length*dt;
t = 0:dt:((time_length(end))-1)/fs;

% define params to plot sepctrogram
channel_id = 62; % Oz  
Hz_id = find(fr==15); % 15Hz

% load data into cell arrays, for each subject
n_subjects = 35;
all_subjects = cell(1,n_subjects);

for n = 1:n_subjects
    subject_name = strcat("S",num2str(n),".mat");
    load(subject_name);
    all_subjects{1,n} = data(electr, idx_first:end-idx_last_minus, :, :); % keep only specified electrodes' recordings and time after first 0.5s : (9,1250,40,6) <--------------------------
    
    % keep a temp to plot the spectrogram
    temp = squeeze(data(channel_id, idx_first:end-idx_last_minus, Hz_id,:)); %(1,1250,1,6) -> (1250,6)
    temp_mean(:,n) = mean(temp,2); % (1250,6) -> (1250,1)
end
temp_mean_all = mean(temp_mean, 2); % (1250,35) -> (1250,1)

% plot the spectrogram and the snr
filt_temp = designfilt('bandpassiir', ...
                        'FilterOrder',2,...
                        'PassbandFrequency1',7,...
                        'PassbandFrequency2',70,...
                        'PassbandRipple',1,...
                        'SampleRate',fs); % bandpass filter [7,70] 
fvtool(filt_temp);
x = filtfilt(filt_temp, temp_mean_all);
n = length(x);
X = (1/n) * fft(x);
X_abs = abs(X);
f = (0:n-1)*(fs/n);
%Xshift = (1/n)*fftshift(X);
%fshift = (-n/2:n/2-1)*(fs/n);
f_window = 5; % df is 0.2Hz, thus [f-1,f+1] is 5 X's on each side
filter_denom = 1/(2*f_window+1)*[ones(f_window,1); 0; ones(f_window,1)];
denom = conv(X_abs, filter_denom, 'same');
snr = 20*log10(X_abs ./ denom);
figure; 
subplot(311); 
plot(t,x); title('x'); xlabel('time (s)')
subplot(312);
plot(f, X_abs, 'DisplayName', sprintf('%.2fHz, channel: %d', fr(Hz_id), channel_id)); 
grid on; grid minor; set(gca,'xtick',[0:5:max(f)]); xlim([0,90]);
xlabel('f'); title('|X|, for averaged over all sessions and subjects')
legend show;
subplot(313); 
plot(f, snr); xlim([0,90]); grid on; title('SNR')


clear x X X_abs f filt_temp;


% construct reference signals for canoncor
for i = 1:n_fr
    y(i,:,:) = refsig(fr(i), fs, full_time, n_h); % 40,10,1250
end

% construct the filter bank of  fbcca
n_fb = 10; % number of filters in the bank
filters = cell(1,n_fb);
for freq_id = 1:n_fb
    bpFilt_temp = designfilt('bandpassiir','FilterOrder',2, ...
                            'PassbandFrequency1',freq_id*8,'PassbandFrequency2',88, ...
                            'PassbandRipple',1,'SampleRate',fs);
    filters{1,freq_id} = bpFilt_temp;
end

%% SSVEP recognition through cca, fbcca, lda-classification
methds = {'cca', 'fbcca'};
n_methds = length(methds);
n_correct = zeros(n_subjects, n_tlength, n_fr, n_methds);
for mthd_id = 2%1:n_methds
    for subj_id = 1:n_subjects % for all subjects
        fprintf('doing subject %d, methd: %d\n', subj_id, mthd_id);
        for sess_id = 1:n_sessions % for every session
            for tl_id = 1:n_tlength % for all time lengths
                indcs = 1:time_length(tl_id); % define current window 
                for fr_id = 1:n_fr % for all frequencies fr
                    x = all_subjects{1, subj_id}(:, indcs, fr_id, sess_id); % (9,window_length)
                    switch methds{mthd_id}
                        case 'cca'
                            pred = cca(x,y,n_fr,indcs);
                        case 'fbcca'
                            pred = fbcca(x,y,n_fr,indcs,filters);
                    end
                    if pred == fr_id % fr_id, is after all, the target id
                       n_correct(subj_id, tl_id, fr_id, mthd_id) = n_correct(subj_id, tl_id, fr_id, mthd_id) + 1;
                    end
                end
            end
        end
        acc = 100 * sum(n_correct(subj_id,:,:,mthd_id), 3) / (n_sessions*n_fr); % accuracy of each subject, for each time length
        disp(['accuracy: ' sprintf(' %.2f', acc)]);
    end
end
save('results', 'n_correct');

%% Plot results 
figure;
for mthd_id = 1:n_methds
    accuracy_per_subject{mthd_id} = squeeze(100 * sum(n_correct(:,:,:,mthd_id),3) / (n_sessions*n_fr)); % 35,20
    acc_mean = mean(accuracy_per_subject{mthd_id});
    acc_std = std(accuracy_per_subject{mthd_id});
    errorbar(TW, acc_mean, acc_std, 'DisplayName', sprintf('mthd: %d', mthd_id))
    hold on;
end
xlabel('Time window (s)');
ylabel('Accuracy (%)');
ylim([0 100]);
title('accuracies');
grid on;
legend show;