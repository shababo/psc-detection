function [all_event_features, events_by_trace] = kmeans_estimate(posteriors,traces,do_plot)

burn_in = 2000;
num_traces = length(posteriors);
colors_groups = hsv(100);
colors_groups = colors_groups(randperm(100),:);
colors_lines =  lines(100);

events_by_trace = cell(length(posteriors),1);

offset = 100;

all_event_features = zeros(0,4);
if do_plot
    figure
end
for i = 1:num_traces
    posterior = truncate_samples(posteriors(i),[burn_in length(posteriors(i).num_events)]);
    k = ceil(mean(posterior.num_events));
    subsample_i = 1:length(posterior.times);
    
    k_opts = k-1:k+1;
    minerr = Inf;
    event_feature_means = [];
    for k_tmp = k_opts
    
        if k_tmp > 0 && length(posterior.amp) >= k_tmp

            samples_matrix = [posterior.amp' posterior.tau1' posterior.tau2' posterior.times'];
%             size(samples_matrix)
    %         
    %         subsample_rate = 10;
    %         subsample_i = randsample(size(samples_matrix,1),ceil(size(samples_matrix,1)/subsample_rate));
    %         samples_matrix = samples_matrix(subsample_i,:);
    %         dpmm_out = dpmm(samples_matrix,25);
    %         assignin('base','dpmm_out',dpmm_out)
    %             good_clusters = find(dpmm_out(end).counts > 25);
    %             event_feature_means = bsxfun(@rdivide,dpmm_out(end).sums(good_clusters,:),dpmm_out(end).counts(good_clusters)');
    %         labels = dpmm_out(end).classes;
            
            [labels_tmp, event_feature_means_tmp] = kmeans(samples_matrix(subsample_i,:)',k_tmp);
            event_feature_means_tmp = event_feature_means_tmp';
            denoised_curve = -1.0*build_curve(event_feature_means_tmp,0,.1,1/20000,2000);
            denoised_curve = denoised_curve - median(denoised_curve);
            this_error = sqrt(mean((denoised_curve - (traces(i,:) - median(traces(i,:)))).^2));
            if this_error < minerr
                minerr = this_error;
                event_feature_means = event_feature_means_tmp;
                labels = labels_tmp;
                k = k_tmp;
            end
        end
        
    end
    
    all_event_features = [all_event_features; event_feature_means];
    events_by_trace{i} = event_feature_means;
%     k
    if do_plot && k > 0 && length(posterior.amp) >= k    
        gscatter(posterior.times(subsample_i),-posterior.amp(subsample_i),labels,colors_lines(i,:),[],1,0)
        hold on
        gscatter(posterior.times(subsample_i),posterior.tau1(subsample_i),labels,colors_lines(i,:),[],1,0)
        hold on
        gscatter(posterior.times(subsample_i),posterior.tau2(subsample_i),labels,colors_lines(i,:),[],1,0)
        hold on
        
        
        scatter(event_feature_means(:,4), -event_feature_means(:,1),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        scatter(event_feature_means(:,4), event_feature_means(:,2),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        scatter(event_feature_means(:,4), event_feature_means(:,3),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        
       

        
        plot(1:2000,traces(i,:) - traces(i,end) - 200 - offset*(i-1),'color',colors_lines(i,:))
        plot(1:2000, -1.0*build_curve(event_feature_means,0,.1,1/20000,2000) - 200 - offset*(i-1),'color',colors_lines(i,:),'Linewidth',2)
        xlim([1 2000])
    elseif do_plot
        plot(1:2000,traces(i,:) - traces(i,end) - 200 - offset*(i-1),'color',colors_lines(i,:))
        xlim([1 2000])
    end

end