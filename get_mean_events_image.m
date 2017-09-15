function [mean_events_image, events_per_trial] = ...
    get_mean_events_image(results_grid, burn_in, do_plot, est_flag)

mean_events_image = zeros(size(results_grid));
events_per_trial = cell(size(results_grid));

for i = 1:size(results_grid,1)
    for j = 1:size(results_grid,2)
        
        num_events_here = 0;
        events_per_trial{i,j} = zeros(length(results_grid{i,j}),1);
        for k = 1:length(results_grid{i,j})
            if est_flag
                events_per_trial{i,j}(k) = length(results_grid{i,j}(k).times);
                num_events_here = num_events_here + length(results_grid{i,j}(k).times);
            else
            posterior = ...
                truncate_samples(results_grid{i,j}(k),...
                [burn_in length(results_grid{i,j}(k).num_events)], [0 300],[10 Inf]);
            est_mean_num_events = mean(posterior.num_events);
            num_events_here = num_events_here + est_mean_num_events;
            events_per_trial{i,j}(k) = est_mean_num_events;
            end
        end
        if length(results_grid{i,j}) > 4
            mean_events_image(i,j) = num_events_here/length(results_grid{i,j});
        end
    end
end

if do_plot
    cmap = parula;
    cmap = [0 0 0; cmap];
    figure;
    imagesc(mean_events_image)
    colormap(cmap)
    colorbar
end