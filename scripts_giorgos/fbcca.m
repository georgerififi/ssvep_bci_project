function pred = fbcca(x,y,n_fr,indcs,filters)

a = 1.25; b = 0.25;
n_filters = length(filters);

for fk_id = 1:n_fr % for all reference frequencies fk
    yk = squeeze(y(fk_id, :, indcs)); % (1,10,window) -> (10,window)
    for sb = 1:n_filters
    	x_temp = filtfilt(filters{1,sb},x); % would be faster to have already implemented filtfilt for the whole signal (not only the current time window)
    										% but would be unfair, as the filtering would have been applied in the whole 5sec recording (longer recording, better fourier...)
    										% but in the online application of the algorithm we will only have the 1sec or 2 sec and so on, recording.  
	    [~,~,r] = canoncorr(x_temp', yk'); % r is (9,1), one value for each harmonic of fk
	    r_sb(1,sb) = r(1); % keep only the max-harmonic response value
	    w_sb(1,sb) = sb^-a + b; % calc the respective weight
	end
	r_all(1,fk_id) = w_sb*((r_sb).^2)'; 
end

[~, pred] = max(r_all); % r_all is (1,40), one value for each frequency fk

end