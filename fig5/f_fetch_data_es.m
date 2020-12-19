function data = f_fetch_data_es(n)

% Load 3-day data
Data = load(n);
Data = Data.Data;
% Data is a cell array of the following order:
% Cell 1: CFP_negative OCT4
% Cell 2: CFP_negative SOX2
% Cell 3: CFP_positive OCT4
% Cell 4: CFP_positive SOX2
% Turn this cell array into two tables
varNames = {'Red', 'Green'};
DATA_neg_full = table(Data{1}, Data{2}, 'VariableNames', varNames);

% Find cells with values above 0 in both columns (CFP-negative)
red_above_idx = find(DATA_neg_full.Red > 0);
green_above_idx = find(DATA_neg_full.Green > 0);
all_above_idx = intersect(red_above_idx, green_above_idx);
data = DATA_neg_full(all_above_idx,:);

end