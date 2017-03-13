function est = get_ests_from_fixed_times(result,event_times,start_est)

est = start_est;
est.amp = zeros(size(event_times));
est.tau1 = zeros(size(event_times));
est.tau2 = zeros(size(event_times));
est.times = event_times;
est.num_events = length(event_times);
for i = 1:length(event_times)
    
    event_window = event_times(i) + [-20 20];
    trunc_samples = truncate_samples(result,[1 length(result.num_events)],event_window,[0 Inf]);
    est.amp(i) = mean(trunc_samples.amp);
    est.tau1(i) = mean(trunc_samples.tau1);
    est.tau2(i) = mean(trunc_samples.tau2);
end