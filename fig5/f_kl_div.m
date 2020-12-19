function kl = f_kl_div(a, b, nbins, pseudocount)
    
    [counts_a, edges] = histcounts(a, nbins);
    counts_b = histcounts(b, edges);

    p_a = (counts_a + pseudocount) ./ sum(counts_a + pseudocount);
    p_b = (counts_b + pseudocount) ./ sum(counts_b + pseudocount);
    
    kl = sum(p_a .* log2(p_a ./ p_b));
end