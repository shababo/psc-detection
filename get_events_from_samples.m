function [events_map, events_peak] = get_events_from_samples(results)

sweeps_window = [3000 5000];
time_window = [10 Inf];
amps_window = [0 Inf];


results_trunc = arrayfun(@(y) ...
truncate_samples(y,sweeps_window,time_window,amps_window),results);

events_map = arrayfun(@(y) get_map_sample(y),results_trunc);

events_peak = [];
% samples_trials_combo = arrayfun(@(x) [x.times],results_trunc,'UniformOutput',0);
% 
% samples_psths = ...
%     cellfun(@(x) histcounts(x,1:20:2001)',samples_trials_combo,'UniformOutput',0);
% 
% % samples_psths = cellfun(@(x) cell2mat(x)',...
% %     samples_psths,'UniformOutput',0);
% 
% new_times = cellfun(@(x) detect_peaks(x',150,10,0,1,0,0,0,0)',...
%     samples_psths,'UniformOutput',0);
% %     assignin('base','new_times',new_times)
% %     return
% events_peak = events_map;
% for j = 1:numel(events_peak)
% 
%     events_peak{j}.times = new_times{j}{1};
% 
% end



