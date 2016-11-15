function [glm_out,xcorr_images,samples_psths] = glm_xcorr_all_pairs(resultsdir,dates,slice_nums,cell_nums,tags, trials_detection)
xcorr_images = 1;
samples_psths = 1;

num_experiments = length(dates);
glm_out = struct();
xcorr_images = cell(num_experiments,1);
samples_psths = cell(num_experiments,1);

sweeps_window = [3000 5000];
time_window = [150 450];
amps_window = [20 Inf];

for i = 1:num_experiments
    
    matchstr =  [dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
        num2str(cell_nums(i)) tags{i} '_ch1']
    rebuildmap_file = ...
        ['old-data/' dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
        num2str(cell_nums(i)) tags{i} '_ch1_trace_grid-rebuildmap.mat'];
    
    results_grid_ch1 = rebuild_detection_grid(...
    resultsdir,matchstr,...
    rebuildmap_file);
    size(results_grid_ch1{1,1})
    trials_detection{i}
    results_grid_trunc_ch1 = ...
    cellfun(@(x) arrayfun(@(y) ...
    truncate_samples(y,sweeps_window,time_window,amps_window),x(trials_detection{i})),...
    results_grid_ch1,'UniformOutput',0);

%     samples_trials_combo_ch1 = cellfun(@(y) arrayfun(@(x) [x.times],y(trials_detection{i}),...
%         'UniformOutput',0),results_grid_ch1,'UniformOutput',0);
%     samples_psths_ch1 = ...
%         cellfun(@(y) cellfun(@(x) histcounts(x,1:1:2000)',y,'UniformOutput',0),...
%         samples_trials_combo_ch1,'UniformOutput',0); 
%     samples_psths_ch1 = cellfun(@(x) smoothts(cell2mat(x)','g',100,20),...
%         samples_psths_ch1,'UniformOutput',0);
%     
%     num_event_array_grid = cellfun(@(y) cell2mat(arrayfun(@(x) mean([x.num_events]),y(trials_detection{i}),...
%         'UniformOutput',0)),results_grid_ch1,'UniformOutput',0);
    glm_out(i).ch1 = run_glm_num_events_spatial(results_grid_trunc_ch1);
%     glm_out(i).ch1 = run_glm_rate_spatial(results_grid_trunc_ch1);


    matchstr =  [dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
        num2str(cell_nums(i)) tags{i} '_ch2']
    rebuildmap_file = ...
        ['old-data/' dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
        num2str(cell_nums(i)) tags{i} '_ch2_trace_grid-rebuildmap.mat'];
    
    results_grid_ch2 = rebuild_detection_grid(...
    resultsdir,matchstr,...
    rebuildmap_file);
    
    results_grid_trunc_ch2 = ...
    cellfun(@(x) arrayfun(@(y) ...
    truncate_samples(y,sweeps_window,time_window,amps_window),x(trials_detection{i})),...
    results_grid_ch2,'UniformOutput',0);

%     samples_trials_combo_ch2 = cellfun(@(y) arrayfun(@(x) [x.times],y(trials_detection{i}),...
%         'UniformOutput',0),results_grid_ch2,'UniformOutput',0);
%     samples_psths_ch2 = ...
%         cellfun(@(y) cellfun(@(x) histcounts(x,1:1:2000)',y,'UniformOutput',0),...
%         samples_trials_combo_ch2,'UniformOutput',0); 
%     samples_psths_ch2 = cellfun(@(x) smoothts(cell2mat(x)','g',100,20),...
%         samples_psths_ch2,'UniformOutput',0);
    
%     assignin('base','samples_psths_ch2',samples_psths_ch2)

%     glm_out(i).ch2 = run_glm_rate_spatial(results_grid_trunc_ch2);
    glm_out(i).ch2 = run_glm_num_events_spatial(results_grid_trunc_ch2);
% 
%     xcorr_images{i} = cellfun(@(x,y) xcorr_peak_trials(x,y,[150 450],0),...
%         samples_psths_ch1,samples_psths_ch2);
% 
%     samples_psths{i} = cell(2,1);
%     samples_psths{i}{1} = samples_psths_ch1;
%     samples_psths{i}{2} = samples_psths_ch2;

end

