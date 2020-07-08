clear; clc; close all;

% first time 
%create_dataset_for_lda 

% if you have done it again then just load it
load('dataset_lda');
n_tlength = length(x_final);
K = 5;


accuracy = zeros(n_tlength, K);
for i = 1:n_tlength
	
	fprintf('running for time length %d\n', i);
	[X,Y] = split_data_to_folds(x_final{i},y_final{i}, K);
	
	for j = 1:K  
		x_test = X{j};
		y_test = Y{j};
		
		X_temp = X;
		X_temp{j} = [];
		x_train = cell2mat(X_temp');
		
		Y_temp = Y;
		Y_temp{j} = [];
		y_train = cell2mat(Y_temp');

		mdl = fitcdiscr(x_train, y_train);

		y_pred = predict(mdl, x_test);
		accuracy(i,j) = sum(y_pred == y_test) / length(y_test);
		fprintf('accuracy for %d fold: %.2f\n', j, accuracy(i,j));
	end
end

accuracy = 100*accuracy;

% plot
figure;
load('dataset_specs');

%errorbar(1:n_tlength, mean(accuracy, 2), std(accuracy,0,2), 'DisplayName', 'lda'); hold on;
errorbar(TW, mean(accuracy, 2), std(accuracy,0,2), 'DisplayName', sprintf('lda, est at %d fold CV', K)); hold on;

% compare with ,,, 

load('cca')
acc_mean = mean(accuracy_per_subject);
acc_std = std(accuracy_per_subject);
errorbar(TW, acc_mean, acc_std, 'DisplayName', 'cca'); hold on;

load('fbcca')
acc_mean = mean(accuracy_per_subject);
acc_std = std(accuracy_per_subject);
errorbar(TW, acc_mean, acc_std, 'DisplayName', 'fbcca');

grid on; ylim([0 100]); xlabel('time window (seconds)'); title('accuracies'); legend show;