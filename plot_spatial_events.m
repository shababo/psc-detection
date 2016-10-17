function [events_grid, sync_score] = plot_spatial_events(posterior_grid,traces)

sync_score = zeros(size(posterior_grid,1),size(posterior_grid,2));
count = 1;
figure(1000)
subplot1(21,21)
plot(1:100,sin(1:100));
for i = 1:size(posterior_grid,1)
    for j = 1:size(posterior_grid,2)
        
        i
        j
        [this_location_events, events_by_trace] = kmeans_estimate(posterior_grid{i,j},traces{i,j},0);
        events_grid{i,j} = events_by_trace;
        figure(999)
        subplot(size(posterior_grid,1),size(posterior_grid,2),count)
        event_vec = zeros(1,100);
        event_vec(ceil(this_location_events(:,4)/20)) = 1;
        plot(smoothts(event_vec,'g',10,3))
        ylim([0 .75])
        axis off
        
        figure(1000)
        subplot1(count)
        scatter(this_location_events(:,4),this_location_events(:,1),'.')
        
        sync_score(i,j) = max(smoothts(event_vec,'g',10,3));
        
        count = count + 1;
    end
end