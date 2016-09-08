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
%         
%         subsample_rate = 10;
%         subsample_i = randsample(size(samples_matrix,1),ceil(size(samples_matrix,1)/subsample_rate));
%         samples_matrix = samples_matrix(subsample_i,:);
%         dpmm_out = dpmm(samples_matrix,25);
%         assignin('base','dpmm_out',dpmm_out)
%             good_clusters = find(dpmm_out(end).counts > 25);
%             event_feature_means = bsxfun(@rdivide,dpmm_out(end).sums(good_clusters,:),dpmm_out(end).counts(good_clusters)');
%         labels = dpmm_out(end).classes;
        [labels, event_feature_means] = kmeans(samples_matrix,k);
        subsample_i = 1:length(posteriors(i).times);
        

        
        gscatter(posteriors(i).times(subsample_i),-posteriors(i).amp(subsample_i),labels,colors_lines(i,:),[],1,0)
        hold on
        gscatter(posteriors(i).times(subsample_i),posteriors(i).tau1(subsample_i),labels,colors_lines(i,:),[],1,0)
        hold on
        gscatter(posteriors(i).times(subsample_i),posteriors(i).tau2(subsample_i),labels,colors_lines(i,:),[],1,0)
        hold on
        
        
        scatter(event_feature_means(:,4), -event_feature_means(:,1),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        scatter(event_feature_means(:,4), event_feature_means(:,2),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        scatter(event_feature_means(:,4), event_feature_means(:,3),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        
       

        all_event_features = [all_event_features; event_feature_means];
    end
     plot(1:2000,traces(i,:) - traces(i,1) - 200 - 25*(i-1),'color',colors_lines(i,:))
     plot(1:2000, -1.0*build_curve(event_feature_means,0,.1,1/20000,2000) - 200 - 25*(i-1),'color',colors_lines(i,:),'Linewidth',2)

     xlim([1 2000])
     event_count = event_count + k;
    
end