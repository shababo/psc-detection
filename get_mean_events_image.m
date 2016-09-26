function mean_events_image = get_mean_events_image(results_grid, burn_in, do_plot)

mean_events_image = zeros(size(results_grid));

for i = 1:size(results_grid,1)
    for j = 1:size(results_grid,2)
        
        num_events_here = 0;
        for k = 1:length(results_grid{i,j})
            posterior = ...
                truncate_samples(results_grid{i,j}(k),...
                [burn_in length(results_grid{i,j}(k).num_events)]);
            k = mean(posterior.num_events);
            num_events_here = num_events_here + k;
        end
        mean_events_image(i,j) = num_events_here/length(results_grid{i,j});
    end
end

if do_plot
    figure;
    imagesc(mean_events_image)
    colormap hot
    colorbar
end