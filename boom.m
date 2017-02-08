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
%%
% 
% clc;
% close all
% exps_to_run = [1 4:5 8 11:12 14:17 19:23 27 29:30];
% exps_to_run = [16:17 19:25 27:30];
% exps_to_run = 12;
% exps_to_run = [22:30];
exps_to_run = [20];
% 
% delete(gcp('nocreate'))
% this_pool = parpool(6);
for ii = 1:length(exps_to_run)

    this_exp = exps_to_run(ii)
%     this_exp = 20;
%     try

        events1 = map_estimates_full{this_exp}{1}{1};
        events2 = map_estimates_full{this_exp}{1}{2};

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
        input_locs1 = find(input_maps_tmp(this_exp).ch1(:) > .1);
        input_locs2 = find(input_maps_tmp(this_exp).ch2(:) > .1);
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
    %     all_locs_ww_q10 = union(input_locs1,input_locs2);
        all_locs_ww_q10 = intersect(input_locs1,input_locs2);
        num_locs_full_ww_q10(this_exp) = length(all_locs_ww_q10);
        num_locs1_full_ww_q10(this_exp) = length(input_locs1);
        num_locs2_full_ww_q10(this_exp) = length(input_locs2);
        num_locs_union_full_ww_q10(this_exp) = length(union(input_locs1,input_locs2));
%         continue
    %     all_locs_ww_q10 = 1:441;
    %     all_locs_ww_q10 = all_locs_ww_q10(26)
        [is,js] = ind2sub([21 21],all_locs_ww_q10);

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


        iters = 10000;
        if ~isempty(all_locs_ww_q10)



            length(all_locs_ww_q10)
    %         all_locs_save(this_exp) = length(all_locs_ww_q10);
            for i = 1:length(all_locs_ww_q10)
                if mod(i,10) == 0
                    i
                end
        %         psth1 = psth1 + sum(events1{is(i),js(i)},1)/norm_factor;
        %         psth2 = psth2 + sum(events2{is(i),js(i)},1)/norm_factor;
                null_dist = zeros(iters,1);
                exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).loc_ind = [is(i) js(i)];
                bin_size = 40;
    %             this_loc_events1 = zeros(size(events1{is(i),js(i)},1),130)';
    %             this_loc_events2 = zeros(size(events1{is(i),js(i)},1),130)';
                this_loc_events1 = zeros(length(events1{is(i),js(i)}),2000/bin_size)';
                this_loc_events2 = zeros(length(events1{is(i),js(i)}),2000/bin_size)';
                
                for j= 1:length(events1{is(i),js(i)})
        %             raw_jpsth = raw_jpsth + events1{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
        %             raw_psth1 = raw_psth1 + events1{is(i),js(i)}(j,:)'/norm_factor*events1{is(i),js(i)}(j,:)/norm_factor;
        %             raw_psth2 = raw_psth2 + events2{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
    %                 [~, event_times] = findpeaks(events1{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events1{is(i),js(i)}(j,20:end)),'MinPeakDistance',10);
                    [event_times, sort_inds] = sort(events1{is(i),js(i)}{j}.times);
%                     this_map_est.amp = this_map_est.amp(sort_inds);
%                     this_map_est.amp(find(diff(this_map_est.times) < 60) + 1) = [];
                    event_times(find(diff(event_times) < 60) + 1) = [];
                    event_times(event_times < 100 | event_times > 800) = [];
%                     event_times = events1{is(i),js(i)}{j}.times;
                    this_loc_events1(ceil((event_times-1)/bin_size),j) = 1;

    %                 [~, event_times] = findpeaks(events2{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events2{is(i),js(i)}(j,20:end)),'MinPeakDistance',10);
%                     event_times = events2{is(i),js(i)}{j}.times;
                    [event_times, sort_inds] = sort(events2{is(i),js(i)}{j}.times);
                    event_times(find(diff(event_times) < 60) + 1) = [];
                    event_times(event_times < 100 | event_times > 800) = [];
                    this_loc_events2(ceil((event_times-1)/bin_size),j) = 1;
                end

                concat1_data = this_loc_events1(:);
                concat2_data = this_loc_events2(:);
                exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).zero_lag_corr = concat1_data(1:end)'*concat2_data(1:end);%+ ...
%                         concat1_data(2:end)'*concat2_data(1:end-1) + ...
%                         concat1_data(1:end-1)'*concat2_data(2:end);
                jitter_window = 120;
                jitter_window_bins = jitter_window/bin_size;
                tmp = repmat(1:ceil(length(concat1_data)/jitter_window_bins),jitter_window_bins,1);
                tmp = tmp(:);
                tmp = tmp(1:length(concat1_data));
                poisson_rate1 = accumarray(tmp(:),concat1_data);
                poisson_rate2 = accumarray(tmp(:),concat2_data);  
                disp('before jittering')
                parfor k = 1:iters
                    event_bins1 = find(poisson_rate1);

                    event_bins2 = find(poisson_rate2);
                    shuffled_events1 = zeros(size(concat1_data));
                    shuffled_events2 = zeros(size(concat2_data));
                    all_bins = union(event_bins1,event_bins2);
                    for l = 1:length(all_bins)
                        this_bin = all_bins(l);
                        num_events1 = poisson_rate1(this_bin);
                        num_events2 = poisson_rate2(this_bin);

                        shuffled_events1(randsample(1:jitter_window_bins,num_events1) + (this_bin-1)*jitter_window_bins) = 1;
                        shuffled_events2(randsample(1:jitter_window_bins,num_events2) + (this_bin-1)*jitter_window_bins) = 1;
                    end
                    shuffled_events1 = shuffled_events1(1:length(concat1_data));
                    shuffled_events2 = shuffled_events2(1:length(concat2_data));
                    null_dist(k) = shuffled_events1(1:end)'*shuffled_events2(1:end);% + ...
%                         shuffled_events1(2:end)'*shuffled_events2(1:end-1) + ...
%                         shuffled_events1(1:end-1)'*shuffled_events2(2:end);
                end
                disp('end jitter')
                exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).null_dist = null_dist;
                [f,x] = ecdf(exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).null_dist);
                idx = find(x < exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).zero_lag_corr,1,'last');
                if ~isempty(idx)
                    if f(idx) == 1
                        exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).p_val = 1/iters;
                    else
                        exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).p_val = 1 - f(idx);
                    end
                else
                    exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).p_val = 1;
                end
                if i == 65435435435
                    break
                end
            end

            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).input_map = ones(21,21);
            for i = 1:length(exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs)
                inds = exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).loc_ind;
                if ~isempty(inds)
                    exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).input_map(inds(1),inds(2)) = exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).p_val;
                end
                if i == 100000
                    break
                end
            end

%             figure;
%             imagesc(exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).input_map)
%             title(['Experiment: ' num2str(this_exp)])

        else
    %         common_input_score_norm_intersect(this_exp) = 0;
            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).null_dist = [];
            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).loc_ind = [];
            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).zero_lag_corr = 0;
            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).input_map = zeros(21,21);
        end
%     catch e
%         disp([num2str(this_exp) ' fail'])
%     end
    
end

% delete(this_pool)



%%

exp_shuffle_stats_event_detection_2017_ww_q10 = struct();
% exps_to_run = [1:12 14:30];
exps_to_run = setdiff(1:31,exclude);
% exps_to_run = 13;
% 
% delete(gcp('nocreate'))
% this_pool = parpool(6);
for ii = 1:length(exps_to_run)

    this_exp = exps_to_run(ii)
%     this_exp = 20;
%     try

        events1 = new_event_detection_2017{this_exp}{1}{1};
        events2 = new_event_detection_2017{this_exp}{1}{2};
%         events1{9,9} = events1{9,9}(1);
%         events2{9,9} = events2{9,9}(1);
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
        input_locs1 = find(input_maps_2017_ww_q10(this_exp).ch1(:) > .1);
        input_locs2 = find(input_maps_2017_ww_q10(this_exp).ch2(:) > .1);
%         input_locs1 = 177;
%         input_locs2 = 177;
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
    %     all_locs_ww_q10 = union(input_locs1,input_locs2);
        all_locs_ww_q10 = intersect(input_locs1,input_locs2);
        num_locs_full_ww_q10(this_exp) = length(all_locs_ww_q10);
        num_locs1_full_ww_q10(this_exp) = length(input_locs1);
        num_locs2_full_ww_q10(this_exp) = length(input_locs2);
        num_locs_union_full_ww_q10(this_exp) = length(union(input_locs1,input_locs2));
%         continue
    %     all_locs_ww_q10 = 1:441;
    %     all_locs_ww_q10 = all_locs_ww_q10(26)
        [is,js] = ind2sub([21 21],all_locs_ww_q10);

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
        dt = 1;
        ts_length = 100/dt;
        iters = 10000;
        if length(exp_shuffle_stats_event_detection_2017_ww_q10) >= this_exp && ~isempty(exp_shuffle_stats_event_detection_2017_ww_q10(this_exp))
         exp_shuffle_stats_event_detection_2017_ww_q10(this_exp) = [];
        end
        if ~isempty(all_locs_ww_q10)



            length(all_locs_ww_q10)
    %         all_locs_save(this_exp) = length(all_locs_ww_q10);
            for i = 1:length(all_locs_ww_q10)
                if mod(i,10) == 0
                    i
                end
        %         psth1 = psth1 + sum(events1{is(i),js(i)},1)/norm_factor;
        %         psth2 = psth2 + sum(events2{is(i),js(i)},1)/norm_factor;
                null_dist = zeros(iters,1);
                exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).loc_ind = [is(i) js(i)];
                bin_size = 1/dt;
    %             this_loc_events1 = zeros(size(events1{is(i),js(i)},1),130)';
    %             this_loc_events2 = zeros(size(events1{is(i),js(i)},1),130)';
                this_loc_events1 = zeros(length(events1{is(i),js(i)}),ts_length/bin_size)';
                this_loc_events2 = zeros(length(events1{is(i),js(i)}),ts_length/bin_size)';
                
                for j= 1:length(events1{is(i),js(i)})
        %             raw_jpsth = raw_jpsth + events1{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
        %             raw_psth1 = raw_psth1 + events1{is(i),js(i)}(j,:)'/norm_factor*events1{is(i),js(i)}(j,:)/norm_factor;
        %             raw_psth2 = raw_psth2 + events2{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
    %                 [~, event_times] = findpeaks(events1{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events1{is(i),js(i)}(j,20:end)),'MinPeakDistance',10);
                    [event_times, sort_inds] = sort(events1{is(i),js(i)}{j}.times);
%                     this_map_est.amp = this_map_est.amp(sort_inds);
%                     this_map_est.amp(find(diff(this_map_est.times) < 60) + 1) = [];
%                     event_times(find(diff(event_times) < 3*dt) + 1) = [];
                    event_times(event_times < 10/dt | event_times > 40/dt) = [];
%                     event_times = events1{is(i),js(i)}{j}.times;
                    this_loc_events1(ceil((event_times-1)/bin_size),j) = 1;

    %                 [~, event_times] = findpeaks(events2{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events2{is(i),js(i)}(j,20:end)),'MinPeakDistance',10);
%                     event_times = events2{is(i),js(i)}{j}.times;
                    [event_times, sort_inds] = sort(events2{is(i),js(i)}{j}.times);
%                     event_times(find(diff(event_times) < 60) + 1) = [];
                    event_times(event_times < 10/dt | event_times > 40/dt) = [];
                    this_loc_events2(ceil((event_times-1)/bin_size),j) = 1;
                end

                concat1_data = this_loc_events1(:);
                concat2_data = this_loc_events2(:);
%                 figure;
%                 ms = 25;
%                 plot(find(concat1_data),zeros(length(find(concat1_data)),1),'.','color',[.4 .4 .4],'markersize',ms); hold on
%                  plot([find(concat1_data)-1 find(concat1_data)+1]',zeros(length(find(concat1_data)),2)','-','color',[.4 .4 .4])
%                 plot(find(concat2_data),zeros(length(find(concat2_data)),1) + .4,'.','color',[0 0 0],'markersize',ms); hold on
%                 plot([find(concat2_data)-1 find(concat2_data)+1]',zeros(length(find(concat2_data)),2)' + .4,'-','color',[0 0 0])
                
%                 exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).zero_lag_corr = concat1_data(1:end)'*concat2_data(1:end)+ ...
%                         concat1_data(2:end)'*concat2_data(1:end-1) + ...
%                         concat1_data(1:end-1)'*concat2_data(2:end) + ...;
%                         concat1_data(3:end)'*concat2_data(1:end-2) + ...
%                         concat1_data(1:end-2)'*concat2_data(3:end);
                cch = smoothts(xcorr(concat1_data,concat2_data)','b',3)*3;
                exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).zero_lag_corr = max(cch(length(concat1_data) - 0:length(concat1_data) + 0));
                jitter_window = 10/dt;
                jitter_window_bins = jitter_window/bin_size;
                tmp = repmat(1:ceil(length(concat1_data)/jitter_window_bins),jitter_window_bins,1);
                tmp = tmp(:);
                tmp = tmp(1:length(concat1_data));
                poisson_rate1 = accumarray(tmp(:),concat1_data);
                poisson_rate2 = accumarray(tmp(:),concat2_data);  
                disp('before jittering')
                parfor k = 1:iters
                    event_bins1 = find(poisson_rate1);

                    event_bins2 = find(poisson_rate2);
                    shuffled_events1 = zeros(size(concat1_data));
                    shuffled_events2 = zeros(size(concat2_data));
                    all_bins = union(event_bins1,event_bins2);
                    for l = 1:length(all_bins)
                        this_bin = all_bins(l);
                        num_events1 = poisson_rate1(this_bin);
                        num_events2 = poisson_rate2(this_bin);

                        shuffled_events1(randsample(1:jitter_window_bins,num_events1) + (this_bin-1)*jitter_window_bins) = 1;
                        shuffled_events2(randsample(1:jitter_window_bins,num_events2) + (this_bin-1)*jitter_window_bins) = 1;
                    end
                    shuffled_events1 = shuffled_events1(1:length(concat1_data));
                    shuffled_events2 = shuffled_events2(1:length(concat2_data));
%                     plot(find(shuffled_events1),zeros(length(find(shuffled_events1)),1) - 2*k,'.','color',[.4 .4 .4],'markersize',ms); hold on
%                     plot([find(shuffled_events1)-1 find(shuffled_events1)+1]',zeros(length(find(shuffled_events1)),2)' - 2*k,'-','color',[.4 .4 .4]); hold on
%                plot(find(shuffled_events2),zeros(length(find(shuffled_events2)),1) + .4 - 2*k,'.','color',[0 0 0],'markersize',ms); hold on
%                 plot([find(shuffled_events2)-1 find(shuffled_events2)+1]',zeros(length(find(shuffled_events2)),2)' + .4 - 2*k,'-','color',[0 0 0])
                    
%                     null_dist(k) = shuffled_events1(1:end)'*shuffled_events2(1:end) + ...
%                         shuffled_events1(2:end)'*shuffled_events2(1:end-1) + ...
%                         shuffled_events1(1:end-1)'*shuffled_events2(2:end) + ...;
%                         shuffled_events1(3:end)'*shuffled_events2(1:end-2) + ...
%                         shuffled_events1(1:end-2)'*shuffled_events2(3:end);
                    cch = smoothts(xcorr(shuffled_events1,shuffled_events2)','b',3)*3;
                    null_dist(k) = max(cch(length(shuffled_events1) - 0:length(shuffled_events1) + 0));
                end
                disp('end jitter')
                exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).null_dist = null_dist;
                [f,x] = ecdf(exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).null_dist);
                idx = find(x < exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).zero_lag_corr,1,'last');
                if ~isempty(idx)
                    if f(idx) == 1
                        exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).p_val = 1/iters;
                    else
                        exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).p_val = 1 - f(idx);
                    end
                else
                    exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).p_val = 1;
                end
                if i == 65435435435
                    break
                end
            end

            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).input_map = ones(21,21);
            for i = 1:length(exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs)
                inds = exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).loc_ind;
                if ~isempty(inds)
                    exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).input_map(inds(1),inds(2)) = exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).p_val;
                end
                if i == 100000
                    break
                end
            end

%             figure;
%             imagesc(exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).input_map)
%             title(['Experiment: ' num2str(this_exp)])

        else
    %         common_input_score_norm_intersect(this_exp) = 0;
            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).null_dist = [];
            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).loc_ind = [];
            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).locs(i).zero_lag_corr = 0;
            exp_shuffle_stats_event_detection_2017_ww_q10(this_exp).input_map = zeros(21,21);
        end
%     catch e
%         disp([num2str(this_exp) ' fail'])
%     end
    
end

% delete(this_pool)