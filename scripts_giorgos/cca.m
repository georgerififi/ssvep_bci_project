function pred = cca(x,y,n_fr,indc)

for fk_id = 1:n_fr % for all reference frequencies fk
    yk = squeeze(y(fk_id, :, indcs)); % (1,10,window) -> (10,window)
    [~,~,r] = canoncorr(x', yk'); % (9,window_length), (10,window_length) -> r is (9,1), one value for each harmonic of fk
    r_all(fk_id) = r(1); % keep only the max-harmonic response value
end
[~, pred] = max(r_all); % r_all is (1,40), one value for each frequency fk