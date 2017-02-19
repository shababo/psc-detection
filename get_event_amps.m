function amps = get_event_amps(result,event_times)

amps = zeros(size(event_times));

for i = 1:length(event_times)
    
    event_window = event_times(i) + [-20 20];
    trunc_samples = truncate_samples(result,[1 length(result.num_events)],event_window,[0 Inf]);
    amps(i) = mean(trunc_samples.amp);
    
end