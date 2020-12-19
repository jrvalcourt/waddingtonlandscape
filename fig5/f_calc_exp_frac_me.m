function exp_frac_me = f_calc_exp_frac_me(logratios, sigmoidA, sigmoidB)
nbins=50;
[counts, bincenters] = hist(logratios, nbins);
sigmoid = f_sigmoid(bincenters, sigmoidA, sigmoidB);
dist = counts ./ sum(counts);
exp_frac_me = sum(sigmoid .* dist);
end