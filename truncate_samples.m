function posterior_trunc = truncate_samples(posterior,new_num_sweeps)


posterior_trunc = struct();
num_new_events = sum(posterior.num_events(1:new_num_sweeps));
posterior_trunc.amp = posterior.amp(1:num_new_events);
posterior_trunc.base = posterior.base(1:new_num_sweeps);
posterior_trunc.tau1 = posterior.tau1(1:num_new_events);
posterior_trunc.tau2 = posterior.tau2(1:num_new_events);
posterior_trunc.num_events = posterior.num_events(1:new_num_sweeps);
posterior_trunc.phi = posterior.phi(1:new_num_sweeps,:);
posterior_trunc.noise = posterior.noise(1:new_num_sweeps);
posterior_trunc.obj = posterior.obj(1:new_num_sweeps);
posterior_trunc.times = posterior.times(1:sum(posterior_trunc.num_events));
