function map_sample = get_map_sample(posterior)
    [~,map_ind] = min(posterior.obj);
    map_sample = ...
            truncate_samples(posterior,[map_ind map_ind]);
end