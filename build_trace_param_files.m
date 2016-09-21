function build_trace_param_files(traces_full,params_base,tracedir,paramdir)

if ~isempty(tracedir)
    tracedir = [tracedir '/'];
end

if ~isempty(paramdir)
    paramdir = [paramdir '/'];
end

if iscell(traces_full)
    [traces_full, rebuild_map] = stack_traces(traces_full);
    rebuild_savename = [tracedir params_base.traces_filename(1:end-4) '-rebuildmap.mat'];
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
    
    save([tracedir params.traces_filename],'traces')
    [~,basename] = fileparts(params.traces_filename);
    
    paramfilename = [paramdir basename '-params.mat'];
    save(paramfilename,'params')
    
end