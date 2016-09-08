function dpcrp = crp_estimate(posteriors,traces)

burn_in = 1;
num_traces = length(posteriors);
colors_groups = hsv(100);
colors_groups = colors_groups(randperm(100),:);
colors_lines =  lines(100);
event_count = 1;

figure

for i = 1:num_traces
    init_k = ceil(mean(posteriors(i).num_events(burn_in:end)))
    
    if init_k > 0
        
        samples_matrix = [posteriors(i).amp' posteriors(i).tau1' posteriors(i).tau2' posteriors(i).times']';
        size(samples_matrix)
    
%         qq0      = estimateGaussianWishart(samples_matrix);
%         alphaa   = 1;
%         alphab   = 1;
%         
%         labels = kmeans(samples_matrix',init_k);
% 
%         
%         dp0      = dp_init(samples_matrix, alphaa, alphab, qq0, labels);
%         numiter  = 1000;
%         tracecrp.numclass  = zeros(1,numiter);
%         dpcrp    = dp0;
%         
%         for iter=1:numiter
% 
%             dp = dpcrp;
%             dp_crp;
%             dp_conparam;
%             dpcrp = dp;
%             dp.numclass
%             tracecrp.numclass(iter) = dp.numclass;
%             
%             cla; hold on;
%             plot(samples_matrix(4,:),samples_matrix(1,:),'.','markersize',10,'linewidth',10);
%             [mu sigma] = map(dp.classqq,dp.classnd);
%             for ii=1:dp.numclass
%                 plotellipse(mu([4 1],ii),sigma([4 1],[4 1],ii),'b','linewidth',3);
%             end
%             axis equal; 
%             hold off
%             drawnow
%         end

%         [y,model] = mixGaussGb(samples_matrix);
          dpmm_fit = dpmm(Y,100);
        
        
%         gscatter(posteriors(i).times,posteriors(i).amp,labels,colors_groups(event_count:event_count+k-1,:),[],[],0)
%         hold on
%         scatter(event_feature_means(:,4), -event_feature_means(:,1),100,colors_lines(i,:),'x','LineWidth',2)
%         hold on
%         scatter(event_feature_means(:,4), event_feature_means(:,2),100,colors_lines(i,:),'x','LineWidth',2)
%         hold on
%         scatter(event_feature_means(:,4), event_feature_means(:,3),100,colors_lines(i,:),'x','LineWidth',2)
%         hold on
%         
%         plot(1:2000,traces(i,:) - traces(i,1) - 200,'color',colors_lines(i,:))
%         event_count = event_count + k;

    end
end