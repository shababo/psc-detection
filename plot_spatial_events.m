function sync_score = plot_spatial_events(posterior_grid)

sync_score = zeros(size(posterior_grid,1),size(posterior_grid,2));
count = 1;
for i = 1:size(posterior_grid,1)
    for j = 1:size(posterior_grid,2)
        
        i
        j
        this_location_events = kmeans_estimate(posterior_grid{i,j},[]);
        
        subplot(size(posterior_grid,1),size(posterior_grid,2),count)
        event_vec = zeros(1,100);
        event_vec(ceil(this_location_events(:,4)/20)) = 1;
        plot(smoothts(event_vec,'g',10,3))
        ylim([0 .75])
        axis off
        
        sync_score(i,j) = max(smoothts(event_vec,'g',10,3));
        
        count = count + 1;
    end
end