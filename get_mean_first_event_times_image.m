function [mean_events_image, events_per_trial] = ...
    get_mean_first_event_times_image(results_grid, burn_in, do_plot, est_flag)

mean_events_image = zeros(size(results_grid));
events_per_trial = cell(size(results_grid));

for i = 1:size(results_grid,1)
    for j = 1:size(results_grid,2)
        
        num_events_here = 0;
        events_per_trial{i,j} = zeros(length(results_grid{i,j}),1);
        event_times = [];
        for k = 1:length(results_grid{i,j})
            if est_flag
                events_per_trial{i,j}(k) = length(results_grid{i,j}{k}.times);
                num_events_here = num_events_here + length(results_grid{i,j}{k}.times);
            else
            posterior = ...
                truncate_samples(results_grid{i,j}(k),...
                [burn_in length(results_grid{i,j}(k).num_events)], [0 300],[7 Inf]);
            map_samp = get_map_sample(posterior);
            
            est_mean_num_events = min(map_samp.times);
            if ~isempty(est_mean_num_events)
%                 num_events_here = num_events_here + est_mean_num_events;
%                 events_per_trial{i,j}(k) = est_mean_num_events;
                event_times = [event_times est_mean_num_events];
            end
            end
        end
        if length(event_times) > 3
            mean_events_image(i,j) = mean(event_times);
        else
            mean_events_image(i,j) = 0;
        end
    end
end

if do_plot
    figure;
    imagesc(mean_events_image)
    colormap hot
    colorbar
end