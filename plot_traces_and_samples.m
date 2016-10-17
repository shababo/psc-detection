function plot_traces_and_samples(traces,posterior)

num_traces = size(traces,1);
colors = lines(num_traces);

burn_in = 1000;

for i = 1:num_traces
    scatter(posterior(i).times',-posterior(i).amp',2,colors(i,:),'.'); hold on;
    scatter(posterior(i).times',posterior(i).tau1',2,colors(i,:),'.'); hold on;
    scatter(posterior(i).times',posterior(i).tau2',2,colors(i,:),'.'); hold on;
    plot(traces(i,:) - traces(i,1) - 200,'Color',colors(i,:)); hold on
end
    
hold off