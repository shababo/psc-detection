function all_event_features = kmeans_estimate(posteriors,traces)

burn_in = 1;
num_traces = length(posteriors);
colors_groups = hsv(100);
colors_groups = colors_groups(randperm(100),:);
colors_lines =  lines(100);
event_count = 1;

all_event_features = zeros(0,4);
figure
for i = 1:num_traces
    k = ceil(mean(posteriors(i).num_events(burn_in:end)))
    
    if k > 0 && length(posteriors(i).amp) > 1
        samples_matrix = [posteriors(i).amp' posteriors(i).tau1' posteriors(i).tau2' posteriors(i).times'];
        size(samples_matrix)

        [labels, event_feature_means] = kmeans(samples_matrix,k);

        
        gscatter(posteriors(i).times,-posteriors(i).amp,labels,colors_lines(i,:),[],[],0)
        hold on
        gscatter(posteriors(i).times,posteriors(i).tau1,labels,colors_lines(i,:),[],[],0)
        hold on
        gscatter(posteriors(i).times,posteriors(i).tau2,labels,colors_lines(i,:),[],[],0)
        hold on
        scatter(event_feature_means(:,4), -event_feature_means(:,1),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        scatter(event_feature_means(:,4), event_feature_means(:,2),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        scatter(event_feature_means(:,4), event_feature_means(:,3),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        
       

        all_event_features = [all_event_features; event_feature_means];
    end
     plot(1:2000,traces(i,:) - traces(i,1) - 200,'color',colors_lines(i,:))
        xlim([1 2000])
        event_count = event_count + k;
    
end