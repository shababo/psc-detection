exps_to_run = [1:12 14:30];
% exps_to_run = 12;
% delete(gcp('nocreate'))
% pool = parpool(3);
for i = 1:length(exps_to_run)
    this_exp = exps_to_run(i)
%     try
%         [glm_out_vb_good3{this_exp}] = ...
        [~,~,map_estimates_full{this_exp}] = ...
            glm_xcorr_all_pairs('/media/shababo/data/old-data',dates(this_exp),slice_nums(this_exp),...
            cell_nums(this_exp),tags(this_exp),trials_detection(this_exp),glm_dummy);
...%            glm_from_map_est(map_estimates(this_exp),glm_dummy);
%                 ROI_VB_som_data(map_estimates(this_exp),[1 pair_id(this_exp)])
%     catch e
%         disp([num2str(this_exp) ' fail'])
%     end
end
% delete(pool);



clc;
% close all
% exps_to_run = [1 4:5 8 11:12 14:17 19:23 27 29:30];
% exps_to_run = [16:17 19:25 27:30];
% exps_to_run = 12;
exps_to_run = [1:12 14:30];
% exps_to_run = [1:15 17:30];
% 
delete(gcp('nocreate'))
this_pool = parpool();
for ii = 1:length(exps_to_run)

    this_exp = exps_to_run(ii)
    
%     try

        events1 = map_estimates{this_exp}{1}{1};
        events2 = map_estimates{this_exp}{1}{2};

    %     events1= samples_psths_new{this_exp}{1};
    %     events2 = samples_psths_new{this_exp}{2};

        disp('got events')
    %     input_locs1 = find(glm_out_newer(this_exp).ch1.glmnet_fit.beta(2:end,find(glm_out_newer(this_exp).ch1.lambda == ...
    %             glm_out_newer(this_exp).ch1.lambda_1se)+2));%+5*(-pair_id(this_exp) + 1)
    %     input_locs2 = find(glm_out_newer(this_exp).ch2.glmnet_fit.beta(2:end,find(glm_out_newer(this_exp).ch2.lambda == ...
    %         glm_out_newer(this_exp).ch2.lambda_1se)+2+5*(-pair_id(this_exp)+1)));

%     ch1_glm_map = zeros(21,21);
%     ch2_glm_map = zeros(21,21);
%     for j = 1:2
%         
%         ch1_glm_map = ch1_glm_map + reshape(glm_to_plot{this_exp}.ch1(j).glmnet_fit.beta(2:end,...
%             find(glm_to_plot{this_exp}.ch1(j).lambda == ...
%             glm_to_plot{this_exp}.ch1(j).lambda_1se)),21,21);
%         
%         ch2_glm_map = ch2_glm_map + reshape(glm_to_plot{this_exp}.ch2(j).glmnet_fit.beta(2:end,...
%             find(glm_to_plot{this_exp}.ch2(j).lambda == ...
%             glm_to_plot{this_exp}.ch2(j).lambda_1se)),21,21);
%         
%     end
    input_locs1 = find(input_maps(this_exp).ch1(:) > .1);
    input_locs2 = find(input_maps(this_exp).ch2(:) > .1);
% 
%         num_nonzeros = unique(sum(glm_to_plot(this_exp).ch1.glmnet_fit.beta > 0));
%         target_num = num_nonzeros(end-1);
%         lambda_ind = find(sum(glm_to_plot(this_exp).ch1.glmnet_fit.beta > 0) == target_num,1,'last');
%         input_locs1 = find(...
%             glm_to_plot(this_exp).ch1.glmnet_fit.beta(2:end,lambda_ind) > .25);
% 
%         num_nonzeros = unique(sum(glm_to_plot(this_exp).ch2.glmnet_fit.beta > 0));
%         target_num = num_nonzeros(end-1);
%         lambda_ind = find(sum(glm_to_plot(this_exp).ch2.glmnet_fit.beta > 0) == target_num,1,'last');
%         input_locs2 = find(...
%             glm_to_plot(this_exp).ch2.glmnet_fit.beta(2:end,lambda_ind) > .25);
    % 
    %     all_locs = union(input_locs1,input_locs2);
        all_locs = intersect(input_locs1,input_locs2);
        num_locs(this_exp) = length(all_locs)
        num_locs1(this_exp) = length(input_locs1);
        num_locs2(this_exp) = length(input_locs2);
        num_locs_union(this_exp) = length(union(input_locs1,input_locs2));
%         continue
    %     all_locs = 1:441;
    %     all_locs = all_locs(26)
        [is,js] = ind2sub([21 21],all_locs);

        null = 1;
    %     if null
    %         iters = 10000;
    %         for j = 1:iters
    %             input_locs1 =  randsample(1:441,length(input_locs1));
    % 
    %             num_event_locs_null(this_exp,j) = length(union(input_locs1,input_locs2));
    %             num_shared_locs_null(this_exp,j) = length(intersect(input_locs1,input_locs2));
    %         end
    %     else
    %         num_event_locs = length(union(input_locs1,input_locs2));
    %         num_shared_locs = length(intersect(input_locs1,input_locs2));
    %     end

        disp('got locs')


        iters = 100000;
        if ~isempty(all_locs)



            length(all_locs)
    %         all_locs_save(this_exp) = length(all_locs);
            for i = 1:length(all_locs)
                if mod(i,10) == 0
                    i
                end
        %         psth1 = psth1 + sum(events1{is(i),js(i)},1)/norm_factor;
        %         psth2 = psth2 + sum(events2{is(i),js(i)},1)/norm_factor;
                null_dist = zeros(iters,1);
                exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).loc_ind = [is(i) js(i)];

    %             this_loc_events1 = zeros(size(events1{is(i),js(i)},1),130)';
    %             this_loc_events2 = zeros(size(events1{is(i),js(i)},1),130)';
                this_loc_events1 = zeros(length(events1{is(i),js(i)}),130)';
                this_loc_events2 = zeros(length(events1{is(i),js(i)}),130)';

                for j= 1:length(events1{is(i),js(i)})
        %             raw_jpsth = raw_jpsth + events1{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
        %             raw_psth1 = raw_psth1 + events1{is(i),js(i)}(j,:)'/norm_factor*events1{is(i),js(i)}(j,:)/norm_factor;
        %             raw_psth2 = raw_psth2 + events2{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
    %                 [~, event_times] = findpeaks(events1{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events1{is(i),js(i)}(j,20:end)),'MinPeakDistance',10);
                    event_times = events1{is(i),js(i)}{j}.times;
                    this_loc_events1(ceil((event_times-1)/10),j) = 1;

    %                 [~, event_times] = findpeaks(events2{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events2{is(i),js(i)}(j,20:end)),'MinPeakDistance',10);
                    event_times = events2{is(i),js(i)}{j}.times;
                    this_loc_events2(ceil((event_times-1)/10),j) = 1;
                end

                concat1_data = this_loc_events1(:);
                concat2_data = this_loc_events2(:);
                exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).zero_lag_corr = concat1_data(1:end-1)'*concat2_data(1:end-1)+ ...
                        concat1_data(2:end)'*concat2_data(1:end-1) + ...
                        concat1_data(1:end-1)'*concat2_data(2:end);

                tmp = repmat(1:ceil(length(concat1_data)/6),6,1);
                tmp = tmp(:);
                tmp = tmp(1:length(concat1_data));
                poisson_rate1 = accumarray(tmp(:),concat1_data);
                poisson_rate2 = accumarray(tmp(:),concat2_data);  
                disp('before jittering')
                parfor k = 1:iters
                    event_bins1 = find(poisson_rate1);
% 
%         num_nonzeros = unique(sum(glm_to_plot(this_exp).ch2.glmnet_fit.beta > 0));
%         target_num = num_nonzeros(end-1);
%         lambda_ind = find(sum(glm_to_plot(this_exp).ch2.glmnet_fit.beta > 0) == target_num,1,'last');
%         input_locs2 = find(...
%             glm_to_plot(this_exp).ch2.glmnet_fit.beta(2:end,lambda_ind) > .25);
                    event_bins2 = find(poisson_rate2);
                    shuffled_events1 = zeros(size(concat1_data));
                    shuffled_events2 = zeros(size(concat2_data));
                    all_bins = union(event_bins1,event_bins2);
                    for l = 1:length(all_bins)
                        this_bin = all_bins(l);
                        num_events1 = poisson_rate1(this_bin);
                        num_events2 = poisson_rate2(this_bin);

                        shuffled_events1(randsample(1:12,num_events1) + (this_bin-1)*12) = 1;
                        shuffled_events2(randsample(1:12,num_events2) + (this_bin-1)*12) = 1;
                    end
                    shuffled_events1 = shuffled_events1(1:length(concat1_data));
                    shuffled_events2 = shuffled_events2(1:length(concat2_data));
                    null_dist(k) = shuffled_events1(1:end-1)'*shuffled_events2(1:end-1) + ...
                        shuffled_events1(2:end)'*shuffled_events2(1:end-1) + ...
                        shuffled_events1(1:end-1)'*shuffled_events2(2:end);
                end
                disp('end jitter')
                exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).null_dist = null_dist;
                [f,x] = ecdf(exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).null_dist);
                idx = find(x < exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).zero_lag_corr,1,'last');
                if ~isempty(idx)
                    exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).p_val = 1 - f(idx);
                else
                    exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).p_val = 1;
                end
                if i == 65435435435
                    break
                end
            end

            exp_shuffle_stats_map_est4_intersect(this_exp).input_map = ones(21,21);
            for i = 1:length(exp_shuffle_stats_map_est4_intersect(this_exp).locs)
                inds = exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).loc_ind;
                if ~isempty(inds)
                    exp_shuffle_stats_map_est4_intersect(this_exp).input_map(inds(1),inds(2)) = exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).p_val;
                end
                if i == 100000
                    break
                end
            end

            figure;
            imagesc(exp_shuffle_stats_map_est4_intersect(this_exp).input_map)
            title(['Experiment: ' num2str(this_exp)])

        else
    %         common_input_score_norm_intersect(this_exp) = 0;
            exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).null_dist = [];
            exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).loc_ind = [];
            exp_shuffle_stats_map_est4_intersect(this_exp).locs(i).zero_lag_corr = 0;
            exp_shuffle_stats_map_est4_intersect(this_exp).input_map = zeros(21,21);
        end
%     catch e
%         disp([num2str(this_exp) ' fail'])
%     end
    
end

delete(this_pool)