function [detect_ests,results_trunc] = run_map_inference_multispot_grid(resultsdir,matchstr,location_inds,sweeps_window,time_window,amps_window)

results_stack = rebuild_detection_grid_multi(resultsdir,matchstr);
thresh = 150/20;
window = 10;


    results_stack_trunc = ...
    arrayfun(@(y) ...
    truncate_samples(y,sweeps_window,time_window,amps_window),[results_stack.trials]);

    map_ests = arrayfun(@(y) get_map_sample(y),results_stack_trunc);


    samples_trials_combo = arrayfun(@(x) [x.times],results_stack_trunc,...
        'UniformOutput',0);
    samples_psths = ...
        cellfun(@(x) histcounts(x,1:1:2001),samples_trials_combo,'UniformOutput',0); 
%     assignin('base','samples_psths',samples_psths)
%     samples_psths = cellfun(@(x) cell2mat(x')',...
%         samples_psths,'UniformOutput',0);

    event_times = cellfun(@(x) detect_peaks(x,thresh,window,0,1,0,0,0)',...
        samples_psths);
    % 
    assignin('base','event_times',event_times);
    assignin('base','map_ests',map_ests);
    % assignin('base','samples_trials_combo',samples_trials_combo);
    % assignin('base','samples_psths',samples_psths);


    detect_ests = map_ests;
    for j = 1:numel(map_ests)

%         detect_ests(j).times = event_times{j};
        detect_ests(j) = get_ests_from_fixed_times(results_stack_trunc(j),event_times{j},detect_ests(j));

    end
    results_trunc = results_stack_trunc;
if ~isempty(location_inds)    
    
% else

%     results_grid = unstack_results_multi(results_stack,location_inds);
% 
%     results_grid_trunc = ...
%     cellfun(@(x) arrayfun(@(y) ...
%     truncate_samples(y,sweeps_window,time_window,amps_window),x),...
%     results_grid,'UniformOutput',0);
% 
%     map_ests = cellfun(@(x) arrayfun(@(y) get_map_sample(y),x,'UniformOutput',0),...
%         results_grid_trunc,'UniformOutput',0);
% 
% 
%     samples_trials_combo = cellfun(@(y) arrayfun(@(x) [x.times],y,...
%         'UniformOutput',0),results_grid_trunc,'UniformOutput',0);
%     samples_psths = ...
%         cellfun(@(y) cellfun(@(x) histcounts(x,1:1:2001)',y,'UniformOutput',0),...
%         samples_trials_combo,'UniformOutput',0); 
%     samples_psths = cellfun(@(x) cell2mat(x')',...
%         samples_psths,'UniformOutput',0);
% 
%     event_times = cellfun(@(x) detect_peaks(x,thresh,window,0,1,0,0,0)',...
%         samples_psths,'UniformOutput',0);
%     % 
%     % assignin('base','event_times',event_times);
%     % assignin('base','map_ests',map_ests);
%     % assignin('base','samples_trials_combo',samples_trials_combo);
%     % assignin('base','samples_psths',samples_psths);
% 
% 
%     detect_ests = map_ests;
%     for j = 1:numel(map_ests)
%     %     size(detect_ests{j})
%     %     size(event_times{j})
%         for k = 1:length(event_times{j})
%             detect_ests{j}{k}.times = event_times{j}{k};
%             detect_ests{j}{k}.amp = get_event_amps(results_grid_trunc{j}(k),event_times{j}{k});
%         end
%     end
%     results_trunc = results_grid_trunc;
    results_trunc = unstack_results_multi(results_trunc,location_inds,0);
    detect_ests = unstack_results_multi(detect_ests,location_inds,0);

end




