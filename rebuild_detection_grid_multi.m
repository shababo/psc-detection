function results_all = rebuild_detection_grid_multi(resultsdir, matchstr, varargin)

if length(varargin) > 1 && ~isempty(varargin{2})
    do_posterior = varargin{2};
else
    do_posterior = 1;
end

filestoload = regexpdir(resultsdir,matchstr);
length(filestoload)

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

if ~isempty(varargin) && ~isempty(varargin{1})
    stim_location_key = varargin{1};
else
    return
end

results_all = unstack_results(results_all, rebuild_map,do_posterior);