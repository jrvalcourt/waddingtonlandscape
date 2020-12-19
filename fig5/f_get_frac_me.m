function frac = f_get_frac_me(data, coords_me, coords_ecto)
    in_me   = inpolygon(data.Red, data.Green, ...
                        coords_me(:,1), coords_me(:,2));
    in_ecto = inpolygon(data.Red, data.Green, ...
                        coords_ecto(:,1), coords_ecto(:,2));
    frac_me   = sum(in_me)   / length(in_me);               
    frac_ecto = sum(in_ecto) / length(in_ecto); 
    
    frac = frac_me / (frac_me + frac_ecto);
end