function posterior_trunc = truncate_samples(posterior,sweep_bounds)


posterior_trunc = struct();
% num_new_events = sum(posterior.num_events(sweep_bounds(1):sweep_bounds(2)));
event_sample_start = max(sum(posterior.num_events(1:sweep_bounds(1)-1)),1);
event_sample_end = sum(posterior.num_events(1:sweep_bounds(2)));
if event_sample_end > 0
    posterior_trunc.amp = posterior.amp(event_sample_start:event_sample_end);
else
    posterior_trunc.amp = [];
end
posterior_trunc.base = posterior.base(sweep_bounds(1):sweep_bounds(2));
if event_sample_end > 0
    posterior_trunc.tau1 = posterior.tau1(event_sample_start:event_sample_end);
    posterior_trunc.tau2 = posterior.tau2(event_sample_start:event_sample_end);
else
    posterior_trunc.tau1 = [];
    posterior_trunc.tau2 = [];
end
posterior_trunc.num_events = posterior.num_events(sweep_bounds(1):sweep_bounds(2));
posterior_trunc.phi = posterior.phi(sweep_bounds(1):sweep_bounds(2),:);
posterior_trunc.noise = posterior.noise(sweep_bounds(1):sweep_bounds(2));
posterior_trunc.obj = posterior.obj(sweep_bounds(1):sweep_bounds(2));
if event_sample_end > 0
    posterior_trunc.times = posterior.times(event_sample_start:event_sample_end);
else
    posterior_trunc.times = [];
end