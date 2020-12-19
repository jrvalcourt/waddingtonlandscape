function bestShift = f_infer_sigmoid_loc(logratios, frac_me, sigmoidA)
    doScore = @(x)scoreShift(x, logratios, frac_me, sigmoidA);
    bestShift = fminbnd(doScore, min(logratios), max(logratios));
end

function sqerror = scoreShift(sigmoidB, logratios, target_frac_me, sigmoidA)
    calc_frac_me = f_calc_exp_frac_me(logratios, sigmoidA, sigmoidB);
    error = (calc_frac_me - target_frac_me);
    sqerror = error .* error;
end
