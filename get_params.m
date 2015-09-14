function params = get_params

%% savename

params.source_path = '/vega/stats/users/bms2156/psc-detection';
params.savepath = '';
params.savename = 'tmp.mat';
params.full_save_string = [params.savepath '/' params.savename];

%% use an rng seed

params.rand = 1;
params.seed = 1234;

%% data params

params.dt = 1/20000;

% direction of events: upward is 1, downard is -1
params.event_sign = -1;

%% subtraces

% first sample, if you want to start at 1, omit
% params.start_ind = 1;

% if you want to go to the end of the traces, omit
% params.duration = 1000;

% if you want all traces, omit
% params.traces_ind = 1:10;

%% inference params

% event amplitude bounds
params.a_max = Inf;
params.a_min = 1;

% baseline bounds
params.b_min = -50;
params.b_max = 50;

% event kernel params.
% params.kernel = @kernel_function;
params.tau1_min = 1/20000;
params.tau1_max = 20/20000;
params.tau2_min = 20/20000;
params.tau2_max = 75/20000;
params.kernel_samples = 3*params.tau2_max;

% poisson/rate
params.p_spike = 2.8e-6;

% ar noise model
params.p = 4;
params.phi_0 = zeros(params.p,1);
params.Phi_0 = 10*eye(params.p); %inverse covariance 3

params.noise_var_init = 5;

%% sampling params

params.num_sweeps = 4000;
params.burn_in_sweeps = 0;

% sampling spike times
params.time_proposal_var = 10;

params.tau1_prop_std = 2/20000;
params.tau1_prop_std = 20/20000;

params.amp_prop_std = .2;

params.baseline_prop_std = 2;

params.add_drop_sweeps = 5;
params.spike_time_sweeps = 3;
params.amp_sweeps = 5;
params.baseline_sweeps = 1;
params.tau1_sweeps = 1;
params.tau2_sweeps = 1;

params.exclusion_bound
params.Dt
params.A
params.b
%% convolution/template-matching initialization method
params.init_method.tau = .002; % min seconds
params.init_method.amp_thresh = 10;
params.init_method.conv_thresh = 1;


