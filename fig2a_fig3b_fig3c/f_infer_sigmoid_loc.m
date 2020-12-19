function bestShift = f_infer_sigmoid_loc(logratios, frac_me, sigmoidA)
    doScore = @(x)scoreShift(x, logratios, frac_me, sigmoidA);
    bestShift = fminbnd(doScore, min(logratios), max(logratios));
end

function sqerror = scoreShift(b, logratios, target_frac_me, sigmoidA)
    nbins=50;
    [counts, bincenters] = hist(logratios, nbins);
    sigmoid = f_sigmoid(bincenters, sigmoidA, b);
    dist = counts ./ sum(counts);
    calc_frac_me = sum(sigmoid .* dist);
    error = (calc_frac_me - target_frac_me);
    sqerror = error .* error;
end
