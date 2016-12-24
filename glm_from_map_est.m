function [glm_out,xcorr_images,samples_psths] = glm_from_map_est(map_estimates,glm_dummy)
xcorr_images = 1;
samples_psths = 1;
glm_out = 1;

num_experiments = length(map_estimates);
glm_out = struct();
xcorr_images = cell(num_experiments,1);
samples_psths = cell(num_experiments,1);

sweeps_window = [3000 5000];
time_window = [150 800];
amps_window = [10 Inf];

time_windows = {[100 800]};
baseline_window = [1500 1999];
bg_trunc = 4/3;

for i = 1:num_experiments
    
    disp('i in glm_xcorr:')
    i
    this_exp_map_ch1 = map_estimates{i}{1}{1};
    this_exp_map_ch2 = map_estimates{i}{1}{2};
%     
%     matchstr =  [dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
%         num2str(cell_nums(i)) tags{i} '_ch1']
%     rebuildmap_file = ...
%         ['old-data/' dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
%         num2str(cell_nums(i)) tags{i} '_ch1_trace_grid-rebuildmap.mat'];
%     
%     results_grid_ch1 = rebuild_detection_grid(...
%     resultsdir,matchstr,...
%     rebuildmap_file);
%     size(results_grid_ch1{1,1})
%     trials_detection{i}
%     results_grid_trunc_ch1 = ...
%     cellfun(@(x) arrayfun(@(y) ...
%     truncate_samples(y,sweeps_window,time_window,amps_window),x(trials_detection{i})),...
%     results_grid_ch1,'UniformOutput',0);
%     

%     samples_psths_ch1 = cellfun(@(x) arrayfun(@get_map_sample,x,'UniformOutput',0),...
%         results_grid_trunc_ch1,'UniformOutput',0);


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
    disp('ch1')
%     ch1(length(time_windows)) = glm_dummy;
    [~, baseline_map] = run_glm_num_events_spatial_mapest(this_exp_map_ch1,baseline_window,1);
    baseline_map = baseline_map/(diff(baseline_window)/20000);
    bg_data = sort(baseline_map(:));
    rate = mean(bg_data(1:ceil(length(bg_data)/bg_trunc)));
    for j = 1:length(time_windows)
        disp('j in for each time window:')
        disp(num2str(j))
        
%         results_grid_trunc_ch1 = ...
%             cellfun(@(x) arrayfun(@(y) ...
%             truncate_samples(y,sweeps_window,time_windows{j},amps_window),x(trials_detection{i})),...
%             results_grid_ch1,'UniformOutput',0);
        [~, resp_map] = run_glm_num_events_spatial_mapest(this_exp_map_ch1,time_windows{j},1);
        
        
        num_trials = length(this_exp_map_ch1{1,1});
        ch1(j).num_trials = num_trials;
        ch1(j).resp_map = resp_map;
        ch1(j).p_val = 1- poisscdf(ch1(j).resp_map - 1,rate*(diff(time_windows{j})/20000));
        ch1(j).rate = rate*(diff(time_windows{j})/20000);
        ch1(j).baseline_map = baseline_map*(diff(time_windows{j})/20000);
    end
    glm_out(i).ch1 = ch1;
    
%     glm_out(i).ch1 = run_glm_rate_spatial(results_grid_trunc_ch1);


%     matchstr =  [dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
%         num2str(cell_nums(i)) tags{i} '_ch2']
%     rebuildmap_file = ...
%         ['old-data/' dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
%         num2str(cell_nums(i)) tags{i} '_ch2_trace_grid-rebuildmap.mat'];
%     
%     results_grid_ch2 = rebuild_detection_grid(...
%     resultsdir,matchstr,...
%     rebuildmap_file);
%     
%     results_grid_trunc_ch2 = ...
%     cellfun(@(x) arrayfun(@(y) ...
%     truncate_samples(y,sweeps_window,time_window,amps_window),x(trials_detection{i})),...
%     results_grid_ch2,'UniformOutput',0);

%     samples_psths_ch2 = cellfun(@(x) arrayfun(@get_map_sample,x,'UniformOutput',0),...
%         results_grid_trunc_ch2,'UniformOutput',0);

%     samples_trials_combo_ch2 = cellfun(@(y) arrayfun(@(x) [x.times],y(trials_detection{i}),...
%         'UniformOutput',0),results_grid_ch2,'UniformOutput',0);
%     samples_psths_ch2 = ...
%         cellfun(@(y) cellfun(@(x) histcounts(x,1:1:2000)',y,'UniformOutput',0),...
%         samples_trials_combo_ch2,'UniformOutput',0); 
%     samples_psths_ch2 = cellfun(@(x) smoothts(cell2mat(x)','g',100,20),...
%         samples_psths_ch2,'UniformOutput',0);
    disp('ch2')
%     ch1(length(time_windows)) = glm_dummy;
    [~, baseline_map] = run_glm_num_events_spatial_mapest(this_exp_map_ch2,baseline_window,1);
    baseline_map = baseline_map/(diff(baseline_window)/20000);
    bg_data = sort(baseline_map(:));
    rate = mean(bg_data(1:ceil(length(bg_data)/bg_trunc)));
    for j = 1:length(time_windows)
        disp('j in for each time window:')
        disp(num2str(j))
        
%         results_grid_trunc_ch1 = ...
%             cellfun(@(x) arrayfun(@(y) ...
%             truncate_samples(y,sweeps_window,time_windows{j},amps_window),x(trials_detection{i})),...
%             results_grid_ch1,'UniformOutput',0);
        [~, resp_map] = run_glm_num_events_spatial_mapest(this_exp_map_ch2,time_windows{j},1);
        num_trials = length(this_exp_map_ch2{1,1});
        ch2(j).num_trials = num_trials;
        ch2(j).resp_map = resp_map;
        ch2(j).p_val = 1 - poisscdf(ch2(j).resp_map-1,rate*(diff(time_windows{j})/20000));
        ch2(j).rate = rate*(diff(time_windows{j})/20000);
        ch2(j).baseline_map = baseline_map*(diff(time_windows{j})/20000);
    end
    glm_out(i).ch2 = ch2;
%     glm_out(i).ch2 = run_glm_rate_spatial(results_grid_trunc_ch2);
%     glm_out(i).ch2 = run_glm_num_events_spatial(results_grid_trunc_ch2);
% 
%     xcorr_images{i} = cellfun(@(x,y) xcorr_peak_trials(x,y,[150 450],0),...
%         samples_psths_ch1,samples_psths_ch2);
% 
    samples_psths{i} = cell(2,1);
%     samples_psths{i}{1} = samples_psths_ch1;
%     samples_psths{i}{2} = samples_psths_ch2;

end
end
function map_sample = get_map_sample(posterior)
    [~,map_ind] = min(posterior.obj);
    map_sample = ...
            truncate_samples(posterior,[map_ind map_ind]);
end