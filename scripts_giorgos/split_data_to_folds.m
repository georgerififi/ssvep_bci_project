function [folds_x, folds_y] = splitDataToFolds(X, y, k)
%
% STRATIFIED CROSS-VALIDATION SPLITING
%
% Outputs :
% - folds_x{1,i}, with i = 1..k contains the variables
% - folds_y{1,i}, ... contains the labels
%
% X, y are distributed among the folds, s.t. if there are 2 categories, 
% then their distribution is the same in every fold.

[m, n] = size(X);

folds = cell(1, k);

n_fold = m/k;
n_mod = mod(m, k);

cat = unique(y);
cat_sorted = sort(cat);
n_cat = length(cat);

X_sorted = [];
y_sorted = [];
n_samples = zeros(1, n_cat); % num_samples per category
for i = 1:n_cat
	idx = (y == cat_sorted(i));
	n_samples(i) = sum(idx ==1);
	X_sorted = [X_sorted ; X(idx, :)];
	y_sorted = [y_sorted ; y(idx, :)];
end

% --------------
% create k-folds s.t. each fold has almost the same distribution i.e. # of occurences.
% -------------

% initialize useful temps
used = zeros(n_cat, 1);
remaining = n_samples;
start_prev = ones(n_cat, 1);
for j = 2:n_cat
	start_prev(j) = start_prev(j-1) + n_samples(j-1);
end

% start distributing samples of each label to every fold, s.t. every label 
% is represented with the same distribution in each fold.
for i = 1:k % for k-folds
	temp_x = [];
	temp_y = [];
	for j = 1:n_cat % for n_cat categories
		start = start_prev(j) + used(j); %  j=1: 1,   j=2: 1+1854=1855. 
		used(j) = ceil(remaining(j) / (k-i+1));% j=1: 1854/(11-1+1) = 169 ,     j=2: 646/(11-1+1)=59.
		stop = start + (used(j) -1); % j=1: 169,   j=2: 1913.
		temp_x = [temp_x ; X_sorted(start:stop, :)]; % j=1: 1..169,   j=2: 1855..1913
		temp_y = [temp_y ; y_sorted(start:stop)];
		remaining(j) = remaining(j) - used(j);
		start_prev(j) = start;
	end
	folds_x{i} = temp_x;
	folds_y{i} = temp_y;
end

% validate that this is indeed the case, at least for the 1st and 2nd categories.
flag=0;
if flag==1
	prompt = 'do you want to validate that indeed all labels are represented with the same distribution in every fold?\n(no:0, yes:1) \n';
	flag = input(prompt); 
	if flag==1
		for i = 1:k
			ratio = sum(folds_y{1, i} == cat_sorted(1)) / sum(folds_y{1, i} == cat_sorted(2)); % ratio of the first to second label
			sprintf('ratio #(label(%d)) : #(label(%d)),  for %d-th fold is: %f', cat_sorted(1), cat_sorted(2), i, ratio)
		end
	end
end

end