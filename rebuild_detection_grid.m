function results_grid = rebuild_detection_grid(resultsdir, matchstr, rebuild_file,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    do_posterior = varargin{1};
else
    do_posterior = 1;
end

filestoload = regexpdir(resultsdir,matchstr);

for i = 1:length(filestoload)
    
    result_struct = load(filestoload{i});
    
    if isfield(result_struct,'results')
        if i == 1
            results_all = result_struct.results;
        else
            results_all = [results_all result_struct.results];
        end
    end
end

size(results_all)

load(rebuild_file)

results_grid = unstack_results(results_all, rebuild_map,do_posterior);