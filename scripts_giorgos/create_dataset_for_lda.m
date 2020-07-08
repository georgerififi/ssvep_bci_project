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
t_start = .5 ; %  .5sec
t_end_minus = .5; % .5sec
idx_first = t_start*fs + 1;
idx_last_minus = t_end_minus*fs;
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

% load data into cell arrays, for each subject
n_subjects = 35;
all_subjects = cell(1,n_subjects);

for n = 1:n_subjects
    subject_name = strcat("S",num2str(n),".mat");
    load(subject_name);
    all_subjects{1,n} = data(electr, idx_first:end-idx_last_minus, :, :); % % keep only specified electrodes' recordings and time after first 0.5s: (9,1250,40,6) <--------------------------
end


% construct reference signals for canoncor
for i = 1:n_fr
    y(i,:,:) = refsig(fr(i), fs, full_time, n_h); % 40,10,1250
end

%% create the datasets 
x_final = cell(1,n_tlength);
y_final = cell(1,n_tlength);
for tl_id = 1:n_tlength % for all time lengths
	indcs = 1:time_length(tl_id); % define current window
	x_final{tl_id} = zeros(n_subjects*n_sessions*n_fr, n_fr); % init
	y_final{tl_id} = zeros(n_subjects*n_sessions*n_fr, 1); % init
	count = 1; % init
	fprintf('doing time length: %d\n', time_length(tl_id));
	for subj_id = 1:n_subjects % for all subjects
    	for sess_id = 1:n_sessions % for all session 
            for fr_id = 1:n_fr % for all frequencies fr
                x = all_subjects{1, subj_id}(:, indcs, fr_id, sess_id); % (9,window_length)
                for fk_id = 1:n_fr % for all reference frequencies fk
				    yk = squeeze(y(fk_id, :, indcs)); % (1,10,window) -> (10,window)
				    [~,~,r] = canoncorr(x', yk'); % (9,window_length), (10,window_length) -> r is (9,1), one value for each harmonic of fk
				    r_all(fk_id) = r(1); % keep only the max-harmonic response value
				end
				x_final{tl_id}(count,:) = r_all;
				y_final{tl_id}(count) = fr_id;
				count = count+1;
            end
        end
    end
end

save('dataset_lda', 'x_final', 'y_final'); 