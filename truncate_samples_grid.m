function posterior_grid_trunc = truncate_samples_grid(posterior_grid,new_num_sweeps)

posterior_grid_trunc = cell(size(posterior_grid));

for i = 1:size(posterior_grid,1)
    for j = 1:size(posterior_grid,2)
        posterior_grid_trunc{i,j} = posterior_grid{i,j};
        for k = 1:length(posterior_grid{i,j})
            posterior_grid_trunc{i,j}(k) = truncate_samples(posterior_grid{i,j}(k),new_num_sweeps);
        end
    end
end