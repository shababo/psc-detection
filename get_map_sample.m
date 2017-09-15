function map_sample = get_map_sample(posterior)
    for i = 1:length(posterior)
        [~,map_ind] = max(posterior(i).obj);
        map_sample(i) = ...
                truncate_samples(posterior(i),[map_ind map_ind]);
    end
end