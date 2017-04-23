function params = build_trace_param_files(traces_full,params_base,tracedir,paramdir,varargin)

if ~isempty(tracedir)
    tracedir = [tracedir '/'];
end

if ~isempty(paramdir)
    paramdir = [paramdir '/'];
end

if ~exist(tracedir,'dir')
    mkdir(tracedir);
end
if ~exist(paramdir,'dir')
    mkdir(paramdir);
end

if ~isempty(varargin)
    params_base.traces_filename = varargin{1};
    params_base.savename = [params_base.traces_filename(1:end-4) '-0000.mat'];
    params_base.full_save_string = [params_base.savename];
end

[~,traces_basename] = fileparts(params_base.traces_filename); 

if iscell(traces_full)
    [traces_full, rebuild_map] = stack_traces(traces_full);
    rebuild_savename = [tracedir traces_basename '-rebuildmap.mat'];
    save(rebuild_savename,'rebuild_map')
end

max_traces_per_job = 160;

num_traces = size(traces_full,1);
num_jobs = ceil(num_traces/max_traces_per_job);

for i = 1:num_jobs
    
    start_i = (i-1)*max_traces_per_job + 1;
    end_i = min(i*max_traces_per_job,num_traces);
    
    traces = traces_full(start_i:end_i,:);
    
    params = params_base;
    params.traces_filename = [params.traces_filename(1:end-4) '-subjob-' sprintf('%03d',i) '.mat'];
    params.savename = [params.savename(1:end-8) '-subjob-' sprintf('%03d',i) '-' params.savename(end-7:end-4) '.mat'];
    params.full_save_string = params.savename;
    
    [~,basename] = fileparts(params.traces_filename);
    save([tracedir basename],'traces')
    
    
    paramfilename = [paramdir basename '-params.mat'];
    save(paramfilename,'params')
    
end