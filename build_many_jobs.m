function build_many_jobs(dates,slice_nums,cell_nums,tags,trials,params,map_index,tracedir,paramdir)

num_experiments = length(dates);

data_files = cell(num_experiments,1);

for i = 1:num_experiments
    
    data_files{i} = ...
        [dates{i} '_slice' num2str(slice_nums(i)) '_cell' num2str(cell_nums(i)) '' tags{i} '.mat'];
    data_files{i}
    load(data_files{i})
    [map_ch1,map_ch2] = see_grid(data,trials{i},map_index,0);
    [dates{i} '_slice' num2str(slice_nums(i)) '_cell' num2str(cell_nums(i)) '' tags{i} '_ch2_trace_grid.mat']
    params.traces_filename = ...
        ['data/' dates{i} '_slice' num2str(slice_nums(i)) '_cell' num2str(cell_nums(i)) '' tags{i} '_ch2_trace_grid.mat'];
    params.savename = ...
        ['data/' dates{i} '_slice' num2str(slice_nums(i)) '_cell' num2str(cell_nums(i)) '' tags{i} '_ch2_trace_grid-0000.mat'];
    params.full_save_string = params.savename;
    build_trace_param_files(map_ch2,params,...
        tracedir,...
        paramdir);
end
end


