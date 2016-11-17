% function output = scratch(linescan)
zero_ind = 6000;
duration = 4000;
fs = 20000;

dirname = '~/projects/mapping/data/0810/';

files = {'0810_cell02.mat','0810_cell04.mat','0810_cell06.mat','0810_cell08.mat','0810_cellXX.mat'};
%,'0810_cell07.mat'
num_cells = length(files);

condition_labels = {'950','900','850','800','750'};

num_conditions = length(condition_labels);

condition_inds = {{[12:16 22:26 47:51],[17:21 43:46],[27:31],[32:36],[37:41]},...
              {[18:24],[25:29],[35:39],[40:44],[45:48]},...
              {[14:18],[19:23],[33:37],[38:42],[43:48]},...
              {[14:18],[19:23],[24:28],[29:34],[35:39]},...
              {[28:35],[36:41],[42:46],[47:52],[53:59]}};

                        %{[7:15],[16:20],[21:25],[26:30],[]},...
dont_plot = 0;
data_ch = 1;

% get cell averages
cell_means = cell(num_cells,num_conditions);
cell_stds = cell(num_cells,num_conditions);
cell_current_max = zeros(num_cells,num_conditions);
cell_charge = zeros(num_cells,num_conditions);
condition_means = cell(num_conditions,1);

for i = 1:num_cells
    
    this_cell_file = [dirname files{i}];
    
    for j = 1:num_conditions
        

        traces = get_sweeps(this_cell_file,data_ch,condition_inds{i}{j},dont_plot);
        cell_means{i,j} = mean(traces);
        cell_means{i,j} = cell_means{i,j} - cell_means{i,j}(zero_ind);
        cell_stds{i,j} = std(traces);
        cell_current_max(i,j) = max(-cell_means{i,j}(zero_ind:zero_ind + duration));
        cell_charge(i,j) = sum(-cell_means{i,j}(zero_ind:zero_ind + duration));

        if i == 1
            condition_means{j} = cell_means{i,j}/num_cells;
        else
            condition_means{j} = condition_means{j} + cell_means{i,j}/num_cells;
        end
        
    end
end


condition_stds = cell(num_condtions,1);
for j = 1:num_conditions
    temp = [];
    for i = 1:num_cells
        temp(i,:) = cell_means{i,j};
    end
    condition_stds{j} = std(temp);
end
%%
normalized_cell_current_max = bsxfun(@ldivide,max(cell_current_max,[],2),cell_current_max);
normalized_cell_charge = bsxfun(@ldivide,max(cell_charge,[],2),cell_charge);



%% plot condition means


time = (0:length(condition_means{1})-1)/fs - zero_ind/fs;

figure;
set(gcf,'defaultAxesColorOrder',jet(num_conditions))
colors = {'r','m','g','c','b'};
for j = num_conditions:-1:1
    shadedErrorBar(time,condition_means{j} - condition_means{j}(zero_ind),condition_stds{j}/sqrt(num_cells),colors{j},1)
    hold on
end
hold off
% legend(fliplr(condition_labels))

ylim([-140 20])
xlim([-.02 .1])

%% plot summary

figure
plot([950 900 850 800 750],mean(normalized_cell_current_max))

figure
plot([950 900 850 800 750],mean(normalized_cell_charge))

%%
[targe_feature_mat_spike3, hot_grouping_spike3] = plot_event_features('data/evoked-pscs-strong-results-0000.mat',get_plot_params);
%%
[targe_feature_mat_spike1, hot_grouping_spike1] = plot_event_features('data/evoked-pscs-strong-results-0000.mat',get_plot_params);
%%
[targe_feature_mat_spike2, hot_grouping_spike2] = plot_event_features('data/evoked-pscs-strong-results-0000.mat',get_plot_params);
%%
amp1 = targe_feature_mat_spike1(hot_grouping_spike1 == 2,1);
amp2 = targe_feature_mat_spike2(hot_grouping_spike2 == 2,1);
amp3 = targe_feature_mat_spike3(hot_grouping_spike3 == 2,1);

bin_centers = 0:5:60;
num_bins = 8;
figure;
% histogram(amp3,bin_centers);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
% hold on
histogram([amp3; amp2],bin_centers);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b','EdgeColor','w','facealpha',0.75)
hold on
histogram(amp1,bin_centers);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','g','EdgeColor','w','facealpha',0.75)

figure;
histogram([amp1; amp2; amp3],bin_centers)


figure;
histogram([targe_feature_mat_spike1(:,1)],bin_centers);

%%

figure;
histogram(targe_feature_mat_spike1(hot_grouping_spike1 == 1 & hot_grouping_spike2 == 1 & hot_grouping_spike3 == 1,1),bin_centers);
hold on
histogram([amp1; amp2; amp3],bin_centers)

%%

amp1 = targe_feature_mat(hot_grouping == 2,1);


bin_centers = 0:2.5:40;
num_bins = 15;
figure;
% histogram(amp3,bin_centers);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
% hold on
histogram(amp1,num_bins);
%%

figure;
histogram([amp1; amp2; amp3],bin_centers)


figure;
histogram([targe_feature_mat_spike1(:,1); targe_feature_mat_spike2(:,1); targe_feature_mat_spike3(:,1)],bin_centers);



figure;
histogram([targe_feature_mat_spike1(hot_grouping_spike1 == 1,1); targe_feature_mat_spike2(hot_grouping_spike2 == 1,1); targe_feature_mat_spike3(hot_grouping_spike3 == 1,1)],bin_centers);
hold on
histogram([amp1; amp2; amp3],bin_centers)

%% save parameter files for cluster

p_spike = [1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2];
a_min = [0 1e-2 1e-1 1e0 1e1 1e2 5e-2 5e-1 5e0 5e1 5e2];
num_sweeps = [1e1 1e2 1e3 1e4];

savefile_basename = 'data/cluster-param-files/params-pspike-%0.0e-amin-%0.0e-num_sweeps-%0.0e.mat';

for i = 1:length(p_spike)
    for j = 1:length(a_min)
        for k = 1:length(num_sweeps) 
            params.cluster = 1;
            params.p_spike = p_spike(i);
            params.a_min = a_min(j);
            params.num_sweeps = num_sweeps(k);
            params = get_params(params);
            savefile = sprintf(savefile_basename,params.p_spike,params.a_min,params.num_sweeps);
            save(strrep(savefile,'+',''),'params')
            clear params
        end
    end
end
    
%% plot ROC

figure
for i = 1:size(roc,1)
    plot(roc(i,:,end,1),roc(i,:,end,2),'-o'); hold on
end
hold off
legend
title('each line is diff a_{min}, within varies p_{spike}')
%%
figure

for i = 1:size(roc,2)
    plot(roc(:,i,1),roc(:,i,2),'-o'); hold on
end
hold off
title('deconv')
% figure
% for i = 1:size(roc_crit,2)
%     plot(roc_crit(:,i,1),roc_crit(:,i,2),'-o'); hold on
% end
% hold off
% title('crit')





%%
    figure
for i = 2:size(roc,1)
    plot(squeeze(roc(i,end,:,1)),squeeze(roc(i,end,:,2)),'-o'); hold on
end
hold off
legend
title('each line is diff p_spike, within varies a_min')

%% plot ROC

figure
for i = 1:size(roc,2)
    plot(roc(3:end,i,1),roc(3:end,i,2),'-o'); hold on
end



hold off
legend
%% all rocs
figure
plot(roc_bayes(3:end,1,1),roc_bayes(3:end,1,2),'.-b','linewidth',2,'markersize',25); hold on

plot(roc_cb(3:end,1,1),roc_cb(3:end,1,2),'.-g','linewidth',2,'markersize',25); hold on
% plot(roc_deconv(:,1,1),roc_deconv(:,1,2),'.-m','linewidth',2,'markersize',25); hold on
plot(roc_wiener(3:end,1,1),roc_wiener(3:end,1,2),'.-r','linewidth',2,'markersize',25); hold on
plot(roc_rand(:,1,1),roc_rand(:,1,2),'--k','linewidth',1,'markersize',5);

hold off

%% false pos

figure
plot(threshold(3:end),roc_cb(3:end,1,1),'.-g','linewidth',2,'markersize',25)
hold on
% plot(threshold,roc_deconv(:,1,1),'.-m','linewidth',2,'markersize',25)
% hold on
plot(threshold(3:end),roc_wiener(3:end,1,1),'.-r','linewidth',2,'markersize',25)
hold on
plot(threshold(3:end),roc_bayes(3:end,1,1),'.-b','linewidth',2,'markersize',25)

%% true pos

figure
plot(threshold(3:end),roc_cb(3:end,1,2),'.-g','linewidth',2,'markersize',25)
hold on
% plot(threshold,roc_deconv(:,1,2),'.-m','linewidth',2,'markersize',25)
% hold on
plot(threshold(3:end),roc_wiener(3:end,1,2),'.-r','linewidth',2,'markersize',25)
hold on
plot(threshold(3:end),roc_bayes(3:end,1,2),'.-b','linewidth',2,'markersize',25)

%% ratio tp/fp


figure
plot(threshold,roc_cb(:,1,2)./(roc_cb(:,1,1) + 1),'.-g','linewidth',2,'markersize',25)
hold on
plot(threshold,roc_deconv(:,1,2)./(roc_deconv(:,1,1) + 1),'.-m','linewidth',2,'markersize',25)
hold on
plot(threshold,roc_wiener(:,1,2)./(roc_wiener(:,1,1) + 1),'.-r','linewidth',2,'markersize',25)
hold on
plot(threshold,roc_bayes(:,1,2)./(roc_bayes(:,1,1) + 1),'.-b','linewidth',2,'markersize',25)


%% ratio fp/tp


figure
plot(threshold,1./(roc_cb(:,1,2)./roc_cb(:,1,1) ),'.-g','linewidth',2,'markersize',25)
hold on
plot(threshold,1./(roc_deconv(:,1,2)./roc_deconv(:,1,1) ),'.-m','linewidth',2,'markersize',25)
hold on
plot(threshold,1./(roc_wiener(:,1,2)./roc_wiener(:,1,1) ),'.-r','linewidth',2,'markersize',25)
hold on
plot(threshold,1./(roc_bayes(:,1,2)./roc_bayes(:,1,1) ),'.-b','linewidth',2,'markersize',25)

%% temporal accuracy histogram

bayes_errs = [];
cb_errs = [];
wiener_errs = [];

for i = 1:10
    
    bayes_errs = [bayes_errs double(squeeze(timing_score_bayes(6,1,i).correct_err))];
    cb_errs = [cb_errs double(squeeze(timing_score_cb(6,1,i).correct_err))];
    wiener_errs = [wiener_errs double(squeeze(timing_score_wiener(7,1,i).correct_err))];
    
end

figure
subplot(131)
histogram(cb_errs,'FaceColor','g','Normalization','pdf')
hold on
plot(ones(2,1)*median(cb_errs),[0 .25],'--k'); hold off
title(['Template Matching, mean error: ' num2str(median(cb_errs)/20) 'msec'])
ylim([0 .25])
xlim([0 20])
set(gca,'xticklabel',{'0','.25','.50','.75','1.0'})
subplot(132)
histogram(wiener_errs,'FaceColor','r','Normalization','pdf')
hold on
plot(ones(2,1)*median(wiener_errs),[0 .25],'--k'); hold off
title(['Wiener Filter, mean error: ' num2str(median(wiener_errs)/20) 'msec'])
ylim([0 .25])
xlim([0 20])
set(gca,'xticklabel',{'0','.25','.50','.75','1.0'})
subplot(133)
histogram(bayes_errs,'FaceColor','b','Normalization','pdf')
hold on
plot(ones(2,1)*median(bayes_errs),[0 .25],'--k'); hold off
title(['Bayesian, mean error: ' num2str(median(bayes_errs)/20) 'msec'])
ylim([0 .25])
xlim([0 20])
set(gca,'xticklabel',{'0','.25','.50','.75','1.0'})

%% simulate random detection

rates = [0 .1 .5 1 2 4 8 16 32];

event_times = cell(length(rates),1);

for i = 1:length(rates)
    event_times{i} = cell(1,10);
    for j = 1:10
        num_events = poissrnd(rates(i));
        event_times{i}{j} = ceil(rand(1,num_events)*20000);
    end
end

save('data/random-detection-results-harder.mat','event_times');

%% plot noise examples noise model figure

% load('data/work/example_noise_traces_work.mat')

figure;
subplot(221)
plot_trace_stack(traces,20,zeros(3,3),'-',[.010 10])
title('Voltage Clamp Recordings')

subplot(222)
plot_trace_stack(ar2_sim_data,20,zeros(3,3),'-',[.010 10])
title('Simulated Noise From Fits - AR(2)')

subplot(223)
plot_trace_stack(ar6_sim_data,20,zeros(3,3),'-',[.010 10])
title('Simulated Noise From Fits - AR(6)')
edi
subplot(224)
plot_trace_stack(ar10_sim_data,20,zeros(3,3),'-',[.010 10])
title('Simulated Noise From Fits - AR(10)')

%% plot individual traces from real data example - fig0

load('/home/shababo/projects/mapping/code/psc-detection/stimfit/real-fit-events.mat')
figure;
legend_names = cell(1,3);

% for i = 1:3
    t = 0:1:.015*20000;
    tau_decay = event1_fit_params.tau1*20; decay = exp(-t/tau_decay);
    tau_rise = event1_fit_params.tau2*20; rise = -exp(-t/tau_rise);
    % plot(decay); hold on; plot(rise)
    plot(-(decay + rise)/max(decay+rise)*14,'linewidth',2)
    hold on; 
    legend_names{1} = ['event 1'];
    
    tau_decay = event2_fit_params.tau1*20; decay = exp(-t/tau_decay);
    tau_rise = event2_fit_params.tau2*20; rise = -exp(-t/tau_rise);
    % plot(decay); hold on; plot(rise)
    plot(-(decay + rise)/max(decay+rise)*19,'linewidth',2)
    hold on; 
    legend_names{2} = ['event 2'];
    
    tau_decay = event3_fit_params.tau1*20; decay = exp(-t/tau_decay);
    tau_rise = event3_fit_params.tau2*20; rise = -exp(-t/tau_rise);
    % plot(decay); hold on; plot(rise)
    plot(-(decay + rise)/max(decay+rise)*16,'linewidth',2)
    hold on; 
    legend_names{3} = ['event 3'];
% end

hold on
legend(legend_names)
axis off

%%

filenames = {'12_1_slice1_cell1.mat',...
             '12_1_slice1_cell2.mat',...
             '12_1_slice2_cell1.mat',...
             '12_1_slice3_cell1.mat',...
             '12_1_slice3_cell2.mat',...
             '12_1_slice4_cell1.mat',...
             '12_1_slice4_cell2.mat',...
             '12_1_slice5_cell1.mat',...
             '12_1_slice5_cell2.mat',...
             '12_2_slice1_cell1.mat',...
             '12_2_slice2_cell1.mat',...
             '12_2_slice2_cell2.mat',...
             '12_3_slice1_cell1.mat',...
             '12_3_slice1_cell2.mat',...
             '12_3_slice1_cell3.mat',...
             '12_3_slice2_cell1.mat',...
             '12_3_slice2_cell2.mat',...
             '12_3_slice3_cell1.mat',...
             '12_3_slice4_cell1.mat',...
             '12_3_slice4_cell2.mat',...
             '12_3_slice5_cell1.mat',...
             '12_3_slice6_cell1.mat',...
             '12_3_slice6_cell2.mat'};
         %%
filenames = {'12_17_slice1chrims_cell1.mat',...
    '12_17_slice1chrims_cell2.mat',...
    '12_17_slice2chrims_cell1.mat',...
    '12_17_slice2chrims_cell2.mat',...
    '12_17_slice2chrims_cell3.mat',...
    '12_17_slice3chrims_cell1.mat',...
    '12_17_slice3chrims_cell2.mat',...
    '12_17_slice4chrims_cell1.mat',...
    '12_17_slice4chrims_cell2.mat',...
    '12_17_slice5chrims_cell1.mat',...
    '12_17_slice5chrims_cell2.mat'};
         
         %%
         
for i = 1:length(filenames)
    compute_relative_obj_position(['data/' filenames{i}],[]);
end

%% Chrimson Good Currents
clear all
% 
filenames = {...%'12_1_slice1_cell1.mat',...
'12_1_slice1_cell2.mat',...
'12_1_slice3_cell1.mat',...
...%'12_1_slice3_cell2.mat',...
'12_1_slice4_cell2.mat',...
'12_1_slice5_cell1.mat',...
'12_1_slice5_cell2.mat',...
'12_2_slice1_cell1.mat',...
'12_3_slice1_cell2.mat',...
...%'12_3_slice2_cell1.mat',...
...% '12_3_slice3_cell1.mat',...
'12_17_slice1chrims_cell1.mat',...
'12_17_slice1chrims_cell2.mat',...
'12_17_slice2chrims_cell1.mat',...
'12_17_slice3chrims_cell1.mat',...
'12_17_slice3chrims_cell2.mat'};


run_count_id = {...%9
2,...
10,...
...%3
13,...
15,...
2,...
6,...
6,...
...%13,...
...%7,...
{3,4},{5},{4,5},{5,6},{4,5}};

%% Chrimson w/ TF

filenames = {'12_19_slice1_cell1.mat',...
    '12_19_slice1_cell2.mat',...
    '12_19_slice1_cell3.mat',...
    '12_20_slice2_cell3.mat'};

run_count_id = {3, 3, 3 ...
    3};

trial_ids1 = [-60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
             0 0 0 0 0 0 0 0 0  -60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -60 -45 -30 -15 0 15 30 45 60]';
         
%% Chrimson w/o TF XYZ

filenames = {'12_22_slice2_cell1.mat',...
    '12_22_slice4_cell1.mat',...
    '12_22_slice4_cell2.mat',...
    '12_22_slice4_cell4.mat'};

run_count_id = {3, {3,4}, 7 ...
    2};

trial_ids1 = [-60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
             0 0 0 0 0 0 0 0 0  -60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -60 -45 -30 -15 0 15 30 45 60]';
         
%% soma-chr2 good currents
clear all
filenames = {'12_15_slice1_cell2.mat','12_15_slice1_cell3.mat','12_15_slice2_cell1.mat','12_15_slice3_cell1.mat',...
    '12_17_slice1_cell1.mat','12_17_slice1_cell3.mat','12_17_slice2_cell1.mat','12_17_slice3_cell3.mat','12_17_slice3_cell4.mat',}

  run_count_id = {4, 3, {4, 5}, {3, 4},{5,6,7},{3,4},{3,4,5},{3,4},4}; 



  
%   run_count_id = 8;
% filenames = {'12_3_slice1_cell2.mat'};
% filenames = {...
% '12_3_slice2_cell1.mat',...
% '12_3_slice3_cell1.mat'};
% 
% run_count_id = [
% 13
% 7];
% 

%% trial ids

trial_ids1 = [-60 -45 -30 -15 0 15 30 45 60 0 0 0 0 0 0 0 0 0 
             0 0 0 0 0 0 0 0 0  -60 -45 -30 -15 0 15 30 45 60
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

trial_ids2 = [
  -64.0000         0         0
  -48.0000         0         0
  -32.0000         0         0
  -16.0000         0         0
   0 0 0
   15.0000         0         0
   30.0000         0         0
   45.0000         0         0
   60.0000         0         0
   0 -64 0
         0  -48.0000         0
         0  -32.0000         0
         0  -16.0000         0
                  0         0         0
         0   16.0000         0
         0   32.0000         0
                  0   48.0000         0
         0   64.0000         0
         ];
     
trial_ids3 = [
  -64.0000         0         0
  -56.0000         0         0
  -48.0000         0         0
  -40.0000         0         0
  -32.0000         0         0
  -24.0000         0         0
  -16.0000         0         0
   -8.0000         0         0
   0 0 0
    7.5000         0         0
   15.0000         0         0
   22.5000         0         0
   30.0000         0         0
   37.5000         0         0
   45.0000         0         0
   52.5000         0         0
   60.0000         0         0
   0 -64 0
            0  -56.0000         0
         0  -48.0000         0
         0  -40.0000         0
         0  -32.0000         0
         0  -24.0000         0
         0  -16.0000         0
         0   -8.0000         0
                  0         0         0
         0    8.0000         0
         0   16.0000         0
         0   24.0000         0
         0   32.0000         0
         0   40.0000         0
                  0   48.0000         0
         0   56.0000         0
         0   64.0000         0
         ];

%%

peak_currents_cells_by_trial = cell(length(filenames),size(trial_ids1,1));
peak_currents_cells_by_trial_mean = zeros(length(filenames),size(trial_ids1,1));
peak_currents_cells_by_trial_std = zeros(length(filenames),size(trial_ids1,1));
peak_currents_trial_mean = zeros(size(trial_ids1,1),1);
peak_currents_trial_std = zeros(size(trial_ids1,1),1);

baseline_window = 20000*[.299 .300]; measure_window = 20000*[.301 .330];

for i = 1:length(filenames)
    
    load(['data/' filenames{i}])
    
    [traces, traces_metadata] = get_sweeps_dir('data',filenames{i},0,1,0,Inf,'run_count',run_count_id{i});
    traces = traces{1};
    
    params1.run_count = run_count_id{i};
    match_inds = match_trials(params1, traces_metadata{1});
    traces = traces(match_inds,:);
    temp = traces_metadata{1};
    traces_metadata = temp(match_inds);
    
    trial_types = zeros(size(traces,1),1);
    

    for j = 1:size(trial_ids1,1)
        
        params2.relative_position = trial_ids1(j,:);
%         params3.relative_position = trial_ids2(j,:);
%         match_inds = unique([match_trials(params2, traces_metadata) match_trials(params3, traces_metadata)]);
        match_inds = unique([match_trials(params2, traces_metadata)]);
        size(match_inds)
        if isempty(match_inds)
            ['data/' filenames{i}]
%             assignin('base','whatthe',traces_metadata)
        end
        trial_types(match_inds) = j;
        
        upper_limit = 500;
        
        peak_currents_cells_by_trial{i,j} = get_current_amp(traces(match_inds,:),baseline_window,measure_window);
        peak_currents_cells_by_trial{i,j}(peak_currents_cells_by_trial{i,j} > upper_limit) = [];
        peak_currents_cells_by_trial_mean(i,j) = mean(peak_currents_cells_by_trial{i,j});
        peak_currents_cells_by_trial_std(i,j) = std(peak_currents_cells_by_trial{i,j});
        
    end


end


for j = 1:length(trial_ids1)
    
    peak_currents_trial_mean(j) = mean(peak_currents_cells_by_trial_mean(:,j));
    peak_currents_trial_std(j) = mean(peak_currents_cells_by_trial_std(:,j));
    
end

%% X Y ONLY

positions = -60:15:60;
colors = lines(length(filenames));
switch_ind = length(peak_currents_trial_std)/2;
figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,1:switch_ind)');
hold on;
for i = 1:length(filenames)
    for j = 1:switch_ind
        plot(positions,peak_currents_cells_by_trial_mean(i,1:switch_ind),'color',colors(i,:))
        scatter(positions(j)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end

end
% legend(filenames)
title('x')
figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,switch_ind+1:end)');
hold on;
for i = 1:length(filenames)
    for j = switch_ind+1:size(trial_ids1,1)
        plot(positions,peak_currents_cells_by_trial_mean(i,switch_ind+1:end),'color',colors(i,:))
        hold on
        scatter(positions(j-switch_ind)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end
end
% legend(filenames)
title('y')

%%

positions = -60:15:60;

%% X Y Z


colors = lines(length(filenames));
switch_ind = length(peak_currents_trial_std)/3;
figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,1:switch_ind)');
hold on;
for i = 1:length(filenames)
    for j = 1:switch_ind
        plot(positions,peak_currents_cells_by_trial_mean(i,1:switch_ind),'color',colors(i,:))
        scatter(positions(j)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end

end
% legend(filenames)
title('x')
figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,switch_ind+1:end)');
hold on;
for i = 1:length(filenames)
    for j = switch_ind+1:switch_ind*2
        plot(positions,peak_currents_cells_by_trial_mean(i,switch_ind+1:switch_ind*2),'color',colors(i,:))
        hold on
        scatter(positions(j-switch_ind)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end
end
% legend(filenames)
title('y')

figure; %h = plot(positions,peak_currents_cells_by_trial_mean(:,switch_ind+1:end)');
hold on;
for i = 1:length(filenames)
    for j = switch_ind*2+1:size(trial_ids1,1)
        plot(positions,peak_currents_cells_by_trial_mean(i,switch_ind*2+1:end),'color',colors(i,:))
        hold on
        scatter(positions(j-switch_ind*2)*ones(length(peak_currents_cells_by_trial{i,j}),1),peak_currents_cells_by_trial{i,j},[],repmat(colors(i,:),length(peak_currents_cells_by_trial{i,j}),1),'filled');
        hold on;
    end
end
% legend(filenames)
title('z')
%% normalization

peak_currents_cells_by_trial_norm = cell(length(filenames),size(trial_ids1,1));
peak_currents_cells_by_trial_mean_norm = zeros(length(filenames),size(trial_ids1,1));
peak_currents_cells_by_trial_std_norm = zeros(length(filenames),size(trial_ids1,1));
peak_currents_trial_mean_norm = zeros(size(trial_ids1,1),1);
peak_currents_trial_std_norm = zeros(size(trial_ids1,1),1);


for i = 1:length(filenames)
    
    for j = 1:size(trial_ids1,1)
        
        start = floor((j-1)/9)*9 + 1
        stop = start + 8
        peak_currents_cells_by_trial_norm{i,j} = ...
            peak_currents_cells_by_trial{i,j}/max(peak_currents_cells_by_trial_mean(i,start:stop));
        peak_currents_cells_by_trial_mean_norm(i,j) = mean(peak_currents_cells_by_trial_norm{i,j});
        peak_currents_cells_by_trial_std_norm(i,j) = std(peak_currents_cells_by_trial_norm{i,j});
        
    end
end
        
for j = 1:length(trial_ids1)
    
    peak_currents_trial_mean_norm(j) = mean(peak_currents_cells_by_trial_mean_norm(:,j));
    peak_currents_trial_std_norm(j) = mean(peak_currents_cells_by_trial_std_norm(:,j));
    
end

%%


switch_ind = length(peak_currents_trial_std)/2;

figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,1:switch_ind)',peak_currents_cells_by_trial_std_norm(:,1:switch_ind)');
title('x')
figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,switch_ind+1:end)',peak_currents_cells_by_trial_std_norm(:,switch_ind+1:end)');
title('y')

%%


switch_ind = length(peak_currents_trial_std)/3;

figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,1:switch_ind)',peak_currents_cells_by_trial_std_norm(:,1:switch_ind)');
title('x')
figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,switch_ind+1:switch_ind*2)',peak_currents_cells_by_trial_std_norm(:,switch_ind+1:switch_ind*2)');
title('y')
figure; errorbar(peak_currents_cells_by_trial_mean_norm(:,switch_ind*2+1:end)',peak_currents_cells_by_trial_std_norm(:,switch_ind*2+1:end)');
title('z')

%%

figure; errorbar(positions,peak_currents_trial_mean(1:9)',peak_currents_trial_std(1:9)');
title('x')
figure; errorbar(positions,peak_currents_trial_mean(10:end)',peak_currents_trial_std(10:end)');
title('y')

figure; errorbar(positions,peak_currents_trial_mean_norm(1:9)',peak_currents_trial_std_norm(1:9)');
title('x norm')
figure; errorbar(positions,peak_currents_trial_mean_norm(10:end)',peak_currents_trial_std_norm(10:end)');
title('y norm')

%%

figure; errorbar(positions,peak_currents_trial_mean(1:9)',peak_currents_trial_std(1:9)');
title('x')
figure; errorbar(positions,peak_currents_trial_mean(10:end)',peak_currents_trial_std(10:end)');
title('y')

%%
figure; errorbar(positions,peak_currents_trial_mean_norm(1:9)',peak_currents_trial_std_norm(1:9)');
title('x norm')
figure; errorbar(positions,peak_currents_trial_mean_norm(10:18)',peak_currents_trial_std_norm(10:18)');
title('y norm')
figure; errorbar(positions,peak_currents_trial_mean_norm(19:end)',peak_currents_trial_std_norm(19:end)');
title('z norm')

%%

clear all12_3_slice1_cell2

load(['data/12_3_slice6_cell2.mat'])
    
[traces, traces_metadata] = get_sweeps_dir('data','12_3_slice6_cell2.mat',0,1,0,Inf,'run_count',8);
traces = traces{1};

params1.run_count = 8;
match_inds = match_trials(params1, traces_metadata{1});
traces = traces(match_inds,:);
temp = traces_metadata{1};
traces_metadata = temp(match_inds);

unique_values = get_unique_metadata_vals(traces_metadata,'relative_position');

%%
traces_by_location = cell(11,11);

start_ind = 20000*.280; end_ind = 20000*.380;

for i = 1:size(traces,1)
    
    ind1 = traces_metadata(i).relative_position(1)/10 + 6;
    ind2 = traces_metadata(i).relative_position(2)/10 + 6;

    traces_by_location{ind1,ind2} = [traces_by_location{ind1,ind2} traces(i,start_ind:end_ind)'];
    
end



%%
%%

trace_count = 0;
for i = 1:size(traces_by_location,1)
    for j = 1:size(traces_by_location,2)
        
        trace_sum = sum(traces_by_location{i,j}',1);
        
        trace_count = trace_count + size(traces_by_location{i,j}',1);
        
    end
end

traces_avg = trace_sum/trace_count;

figure; plot(traces_avg)
     



%%

    
           filenames = {'12_3_slice1_cell1.mat',...
             '12_3_slice1_cell2.mat',...
             '12_3_slice1_cell3.mat',...
             '12_3_slice2_cell1.mat',...
             '12_3_slice2_cell2.mat',...
             '12_3_slice3_cell1.mat',...
             '12_3_slice4_cell1.mat',...
             '12_3_slice4_cell2.mat',...
             '12_3_slice5_cell1.mat',...
             '12_3_slice6_cell1.mat',...
             '12_3_slice6_cell2.mat'};
         
         
         
%% testing wiener filter

trace = traces{1,1}(1,:);
nfft = length(trace) + length(template);
dt = 1/20000;
threshold = 2.0;
min_window = 20;

ar_noise_params.sigma_sq = 2.948727926352792;
ar_noise_params.phi = [1.000000000000000, -0.982949319747574, 0.407063852831604];
gamma = 1e6;

figure;
subplot(1,2,1)
[filtered_trace, event_times_init] = wiener_filter(trace,template,ar_noise_params,...
            nfft, 1/20000, threshold, min_window);
        
subplot(1,2,2)
plot(filtered_trace/max(filtered_trace)*5)
hold on
plot(trace + 40 - trace(1))
hold on
scatter(event_times_init, 8*ones(1,length(event_times_init)))

%% spike detection on grid
 
trace_grids_3_31_s2c2_r2_3 = {traces_by_location_3_31_s2c2_r3_5mw, traces_by_location_3_31_s2c2_r3_10mw, traces_by_location_3_31_s2c2_r3_15mw,...
    traces_by_location_3_31_s2c2_r3_25mw, traces_by_location_3_31_s2c2_r3_50mw, traces_by_location_3_31_s2c2_r3_100mw};

trace_grids_3_29_s1c2_r2 = {traces_by_location_3_29_s1c2_r2_25mw, traces_by_location_3_29_s1c2_r2_50mw, traces_by_location_3_29_s1c2_r2_100mw};

trace_grids_3_29_s1c4_r2 = {traces_by_location_3_29_s1c4_r2_25mw, traces_by_location_3_29_s1c4_r2_50mw, traces_by_location_3_29_s1c4_r2_100mw};

trace_grids_3_31_s1c1_r4_5 = {traces_by_location_3_31_s1c1_r4_5_25mw, traces_by_location_3_31_s1c1_r4_5_50mw, traces_by_location_3_31_s1c1_r4_5_100mw};

trace_grids_3_31_s1c2_r4_5 = {traces_by_location_3_31_s1c2_r5_5mw, traces_by_location_3_31_s1c2_r5_10mw, traces_by_location_3_31_s1c2_r5_15mw,...
    traces_by_location_3_31_s1c2_r4_25mw, traces_by_location_3_31_s1c2_r4_50mw, traces_by_location_3_31_s1c2_r4_100mw};

trace_grids_4_5_s2c1_r5 = {traces_by_location_4_5_s2c1_r5_25mw, traces_by_location_4_5_s2c1_r5_50mw, traces_by_location_4_5_s2c1_r5_100mw};

trace_grids_4_6_s3c2_r1 = {traces_by_location_4_6_s3c2_r1_25mw, traces_by_location_4_6_s3c2_r1_50mw, traces_by_location_4_6_s3c2_r1_100mw};

trace_grids_4_6_s3c5_r1 = {traces_by_location_4_6_s3c5_r1_25mw, traces_by_location_4_6_s3c5_r1_50mw, traces_by_location_4_6_s3c5_r1_100mw};

trace_grids_4_6_s3c7_r2 = {traces_by_location_4_6_s3c7_r2_25mw, traces_by_location_4_6_s3c7_r2_50mw, traces_by_location_4_6_s3c7_r2_100mw};

trace_grids_4_6_s3c8_r3 = {traces_by_location_4_6_s3c8_r3_25mw, traces_by_location_4_6_s3c8_r3_50mw, traces_by_location_4_6_s3c8_r3_100mw};

%%

trace_grids = trace_grids_4_6_s3c8_r3;

detection_grids = cell(size(trace_grids));

for i = 1:length(trace_grids)
    
    trace_grid_tmp = trace_grids{i};
    [traces_tmp, rebuild_map] = stack_traces(trace_grid_tmp);

    detection_results = detect_peaks(-1.0*bsxfun(@minus,traces_tmp,median(traces_tmp,2)),0.1,20,1,1,0)*70;
    detection_grids{i} = unstack_traces(detection_results,rebuild_map);
    
end
    
detection_results_4_6_s3c8_r3 = detection_results;
detection_grids_4_6_s3c8_r3 = detection_grids;


figure; compare_trace_stack_grid({trace_grids{:},detection_grids_4_6_s3c8_r3{:}},...
    5,1,[],0,{'25 mW', '50 mW', '100 mW'},2)

%% count spikes and get means

spike_counts = zeros([size(detection_grids{1}) length(detection_grids)]);
max_val = 0;

axs = [];

figure
colormap hot
for i = 1:length(detection_grids)
    
    this_grid = detection_grids{i};
    for j = 1:size(this_grid,1)
        for k = 1:size(this_grid,2)
            
            spike_counts(j,k,i) = length(find(this_grid{j,k}(:)))/size(this_grid{j,k},1);
            if spike_counts(j,k,i) > max_val
                max_val = spike_counts(j,k,i);
            end
            
        end
    end
    subplot(1,length(detection_grids),i)
%     subplot(2,ceil(length(detection_grids)/2),i)
    imagesc(spike_counts(:,:,i));
    axis square
    axis off
end

for i = 1:length(detection_grids)
    subplot(1,length(detection_grids),i)
%     subplot(2,ceil(length(detection_grids)/2),i)
    caxis([0 max_val])
end

spike_counts_4_6_s3c8_r3 = spike_counts;

%% delay times and get means

delays = zeros([size(detection_grids{1}) length(detection_grids)]);
max_val = 0;

axs = [];

figure

colormap hot
for i = 1:length(detection_grids)
    
    this_grid = detection_grids{i};
    for j = 1:size(this_grid,1)
        for k = 1:size(this_grid,2)
            
            these_delays = [];
            
            for m = 1:size(this_grid{j,k},1)
                these_delays = [these_delays find(this_grid{j,k}(m,:),1,'first')/20000 - .005];
            end
            
            delays(j,k,i) = mean(these_delays);
            
            if delays(j,k,i) > max_val
                max_val = delays(j,k,i);
            end
            
        end
    end
    subplot(2,ceil(length(detection_grids)/2),i)
    pcolor(delays(:,:,i));
    axis ij
    axis square
    axis off
end

delays_4_6_s3c8_r3 = delays;

for i = 1:length(detection_grids)
    subplot(2,ceil(length(detection_grids)/2),i)
    caxis([0 max_val])
end

%%

figure 

for i = 1:3
    
    subplot(3,3,i)
    pcolor(flipud(delays_3_31_s2c2_r2_3(:,:,i+3)));
    caxis([0 .05])
    axis square
    axis off
end

for i = 1:3
    
    subplot(3,3,i+3)
    pcolor(flipud(delays_3_29_s1c2_r2(:,:,i)));
    caxis([0 .05])
    axis square
    axis off
end


for i = 1:3
    
    subplot(3,3,i+6)
    pcolor(flipud(delays_3_31_s1c2_r4_5(:,:,i+3)));
    caxis([0 .05])
    axis square
    axis off
end


colormap hot

%%

figure 

for i = 1:3
    
    subplot(3,3,i)
    imagesc(spike_counts_4_6_s3c5_r1(:,:,i));
    caxis([0 2])
    axis square
    axis off
end

for i = 1:3
    
    subplot(3,3,i+3)
    imagesc(spike_counts_4_6_s3c2_r1(:,:,i));
    caxis([0 2])
    axis square
    axis off
end


for i = 1:3
    
    subplot(3,3,i+6)
    imagesc(spike_counts_4_6_s3c7_r2(:,:,i));
    caxis([0 2])
    axis square
    axis off
end
colormap hot

%%
figure
for i = 1:6
    
    subplot(2,6,i)
    imagesc(spike_counts_3_31_s2c2_r2_3(:,:,i));
    caxis([0 3])
%     axis square
    axis off
end

for i = 1:6
    
    subplot(2,6,i + 6)
    imagesc(spike_counts_3_31_s1c2_r4_5(:,:,i));
    caxis([0 3])
%     axis square
    axis off
end



colormap hot

%%
baseline_window = [1800 2000];
measure_window = [2000 3000];

amps_50 = get_current_amp(traces_s1c1_50mw,baseline_window,measure_window);
amps_75 = get_current_amp(traces_s1c1_75mw,baseline_window,measure_window);
amps_100 = get_current_amp(traces_s1c1_100mw,baseline_window,measure_window);
amps_125 = get_current_amp(traces_s1c1_125mw,baseline_window,measure_window);
amps_150 = get_current_amp(traces_s1c1_150mw,baseline_window,measure_window);
amps_175 = get_current_amp(traces_s1c1_175mw,baseline_window,measure_window);

baseline_window = [3800 4000];
measure_window = [4000 5000];

amps2_50 = get_current_amp(traces_s2c2_50mw,baseline_window,measure_window);
% amps2_75 = get_current_amp(traces_s2c2_75mw,baseline_window,measure_window);
amps2_100 = get_current_amp(traces_s2c2_100mw,baseline_window,measure_window);
amps2_125 = get_current_amp(traces_s2c2_125mw,baseline_window,measure_window);
amps2_150 = get_current_amp(traces_s2c2_150mw,baseline_window,measure_window);
amps2_175 = get_current_amp(traces_s2c2_175mw,baseline_window,measure_window);

figure; plot([50 75 100 125 150 175], [mean(amps_50) mean(amps_75) mean(amps_100) mean(amps_125) mean(amps_150) mean(amps_175)])
hold on; plot([50 100 125 150 175], [mean(amps2_50) mean2(amps2_100) mean(amps2_125) mean(amps2_150) mean(amps2_175)])

%%
trace_grid_ch1 = traces_by_location_5_12_s2c1_2_r4{1};
trace_grid_ch2= traces_by_location_5_12_s2c1_2_r4{2};

across_ch_corr_image = zeros(size(trace_grid_ch1));

for i = 1:size(trace_grid_ch1,1)
    for j = 1:size(trace_grid_ch1,2)
        for k = 1:size(trace_grid_ch1{i,j},1)
            corr_mat = corr([trace_grid_ch1{i,j}(k,:); trace_grid_ch2{i,j}(k,:)]');
            across_ch_corr_image(i,j) = across_ch_corr_image(i,j) + corr_mat(1,2);
        end
        across_ch_corr_image(i,j) = across_ch_corr_image(i,j)/size(trace_grid_ch1{i,j},1);
    end
end

figure; 
imagesc(across_ch_corr_image)
colormap hot
colorbar

%%

time_posteriors = zeros(length(results),2000);


for i = 1:length(results)
    time_posteriors(i,:) = histcounts(results(i).trials.times,0:2000);
end

time_posteriors2 = zeros(length(results),2000);


for i = 1:length(results2)
    time_posteriors2(i,:) = histcounts(results2(i).trials.times,0:2000);
end
%%

event_timeseries2 = get_event_times_init(results2,2000,1,10);
event_timeseries1 = get_event_times_init(results,2000,1,5);
event_timeseries1_smooth = smoothts(time_posteriors,'g',100,20);
event_timeseries2_smooth = smoothts(time_posteriors2,'g',100,20);
% time_posteriors(:,1:50) = 0;
events_ts_grid1_smooth = unstack_traces(event_timeseries1_smooth/5,params.rebuild_map);
events_ts_grid2_smooth = unstack_traces(event_timeseries2_smooth/5,params.rebuild_map);
events_ts_grid1_smooth = unstack_traces(event_timeseries1_smooth*300,params.rebuild_map);

events_ts_grid2_smooth = unstack_traces(event_timeseries2_smooth*300,params.rebuild_map);
figure; compare_trace_stack_grid_overlap({events_ts_grid1_smooth,events_ts_grid2_smooth},3,1,[],0,{'L4','L5'},1)

%%
event_timeseries2 = get_event_times_init(results_5_12_s2c2_r4_tracegrid,2000,1,10);
event_timeseries1 = get_event_times_init(results_5_12_s2c1_r4_tracegrid,2000,1,10);
event_timeseries1_smooth = smoothts(event_timeseries1,'g',100,20);
event_timeseries2_smooth = smoothts(event_timeseries2,'g',100,20);
events_ts_grid1_smooth = unstack_traces(event_timeseries1*50,params.rebuild_map);
events_ts_grid2_smooth = unstack_traces(event_timeseries2_smooth*1000,params.rebuild_map);
figure; compare_trace_stack_grid([traces_by_location_5_12_s2c1_2_r4(:); {events_ts_grid1_smooth}; {events_ts_grid2_smooth}],3,1,[],0,{'L4','L5'},1)
figure; compare_trace_stack_grid_overlap({events_ts_grid1_smooth,traces},3,1,[],0,{'L4','L5'},1)

% figure; compare_trace_stack_grid_overlap({events_ts_grid1_smooth,events_ts_grid2_smooth},3,1,[],0,{'L4','L5'},1)

%%
event_timeseries2 = zeros(length(results_5_12_s2c1_r4_tracegrid),2000);
event_timeseries1 = zeros(length(results_5_12_s2c2_r4_tracegrid),2000);

for i = 1:length(results_5_12_s2c1_r4_tracegrid)
    event_timeseries1(i,:) = histcounts(results_5_12_s2c1_r4_tracegrid(i).trials.times,0:1:2000);
    event_timeseries2(i,:) = histcounts(results_5_12_s2c2_r4_tracegrid(i).trials.times,0:1:2000);
    event_timeseries1(i,1:20) = 0;
    event_timeseries2(i,1:20) = 0;
end

events_ts_grid1_smooth = unstack_traces(event_timeseries1/10,params.rebuild_map);
events_ts_grid2_smooth = unstack_traces(event_timeseries2/10,params.rebuild_map);
% figure; compare_trace_stack_grid({events_ts_grid1_smooth,events_ts_grid2_smooth},3,1,[],0,{'L4','L5'},1)



%%


max_xcorr = cell(size(detection_grid_5_12_s2c1_r4_2000));
mad_xcorr_lag = cell(size(detection_grid_5_12_s2c1_r4_2000));

max_xcorr_img = zeros(size(detection_grid_5_12_s2c1_r4_2000));
mad_xcorr_lag_img = zeros(size(detection_grid_5_12_s2c1_r4_2000));

for i = 1:size(detection_grid_5_12_s2c1_r4_2000,1)
    for j = 1:size(detection_grid_5_12_s2c1_r4_2000,2)
        max_xcorr{i,j} = zeros(size(detection_grid_5_12_s2c1_r4_2000{i,j},1),1);
        mad_xcorr_lag{i,j} = zeros(size(detection_grid_5_12_s2c1_r4_2000{i,j},1),1);
        trace1_tmp = [];
        trace2_tmp = [];
        for k = 1:size(detection_grid_5_12_s2c1_r4_2000{i,j},1)
            trace1_tmp = [trace1_tmp smoothts(detection_grid_5_12_s2c1_r4_2000{i,j}(k,:),'g',100,10)];
            trace2_tmp = [trace2_tmp smoothts(detection_grid_5_12_s2c2_r4_2000{i,j}(k,:),'g',100,10)];
        end
        full_xcorr = xcorr(trace1_tmp,trace2_tmp,'coeff');
        [max_xcorr_img(i,j), mad_xcorr_lag_img(i,j)] = max(full_xcorr(9900:10100));
    end
end

%%

max_xcorr = cell(size(detection_grid_5_12_s2c1_r4_2000));
mad_xcorr_lag = cell(size(detection_grid_5_12_s2c1_r4_2000));

max_xcorr_img_shuff = zeros(size(detection_grid_5_12_s2c1_r4_2000));
mad_xcorr_lag_img = zeros(size(detection_grid_5_12_s2c1_r4_2000));

shuffle_rows = randperm(size(detection_grid_5_12_s2c1_r4_2000,1));
shuffle_cols = randperm(size(detection_grid_5_12_s2c1_r4_2000,2));
shuffle_trials = randperm(size(detection_grid_5_12_s2c1_r4_2000{1,1},1));
while isequal(1:size(detection_grid_5_12_s2c1_r4_2000{1,1},1),shuffle_trials)
    shuffle_trials = randperm(size(detection_grid_5_12_s2c1_r4_2000{1,1},1));
end

for i = 1:size(detection_grid_5_12_s2c1_r4_2000,1)
    for j = 1:size(detection_grid_5_12_s2c1_r4_2000,2)
        max_xcorr{i,j} = zeros(size(detection_grid_5_12_s2c1_r4_2000{i,j},1),1);
        mad_xcorr_lag{i,j} = zeros(size(detection_grid_5_12_s2c1_r4_2000{i,j},1),1);
        trace1_tmp = [];
        trace2_tmp = [];
        for k = 1:size(detection_grid_5_12_s2c1_r4_2000{i,j},1)
            trace1_tmp = [trace1_tmp smoothts(detection_grid_5_12_s2c1_r4_2000{i,j}(k,:),'g',100,10)];
            trace2_tmp = [trace2_tmp smoothts(detection_grid_5_12_s2c2_r4_2000{i,j}(shuffle_trials(k),:),'g',100,10)];
        end
        full_xcorr = xcorr(trace1_tmp,trace2_tmp,'coeff');
        [max_xcorr_img_shuff(i,j), mad_xcorr_lag_img(i,j)] = max(full_xcorr(9900:10100));
    end
end

%%

max_xcorr_img = zeros(size(events_ts_grid1_smooth));
mad_xcorr_lag_img = zeros(size(events_ts_grid1_smooth));

for i = 1:size(events_ts_grid1_smooth,1)
    for j = 1:size(events_ts_grid1_smooth,2)
        [max_xcorr_img(i,j), max_i] = max(max_xcorr{i,j});
        mad_xcorr_lag_img(i,j) = mad_xcorr_lag{i,j}(max_i);

    end
end

%%


figure;
subplot(121)
imagesc(max_xcorr_img)
colorbar
subplot(122)
% mad_xcorr_lag_img(mad_xcorr_lag_img >= 5060) = 0;
imagesc(mad_xcorr_lag_img)
% caxis([4940 5060])
colorbar

colormap hot

%%

figure; histogram(max_xcorr_img_shuff(:))

%%

max_xcorr_img_sig = max_xcorr_img;
mad_xcorr_lag_img_sig = mad_xcorr_lag_img;
max_xcorr_img_sig(max_xcorr_img_sig < quantile(max_xcorr_img_shuff(:),.0)) = 0;
mad_xcorr_lag_img_sig(max_xcorr_img_sig < quantile(max_xcorr_img_shuff(:),.0)) = 0;

figure;
subplot(121)
imagesc(max_xcorr_img_sig)
axis off
title('Correlations')
colorbar
subplot(122)
% mad_xcorr_lag_img(mad_xcorr_lag_img >= 5060) = 0;
imagesc(mad_xcorr_lag_img_sig - 100)
axis off
title('Max Correlation Lag (samples)')
% caxis([4940 5060])
colorbar

colormap hot


%%


xcorrs = cell(size(events_ts_grid1_smooth));

for i = 1:size(events_ts_grid1_smooth,1)
    for j = 1:size(events_ts_grid1_smooth,2)
        xcorrs{i,j} = zeros(size(events_ts_grid1_smooth{i,j}));
        for k = 1:size(events_ts_grid1_smooth{i,j},1)
            xcorrs{i,j}(k,:) = max(xcorr(events_ts_grid1_smooth{i,j}(k,:),events_ts_grid1_smooth{i,j}(k,:)));
        end
    end
end

%%
results_5_12_s2c1_r4_tracegrid


%%

load('/media/shababo/Layover/fly/High Res MaxProj of few z sections (4x more frames and acquisition rate ~3x faster)/Fly2d_Day1_FemaleV_MB247lexA_LexAopGC6s_20160712_40x_zoomed_200TPsift_snmu_results_r50.mat')
good_comps = [4 6 7 9 10 12 13 14 16 17 19 20 22 24 25 26];
edgesize = sqrt(size(spatial,1));
figure
set(gcf,'position',[66 1 1855 1121])

for j = 1:size(temporal,2)
    count = 1;
    for i = good_comps
        subplot(4,4,count)
        imagesc(reshape(spatial(:,i),[edgesize edgesize])*temporal(i,j)); colormap gray
        axis off
        axis tight
        title(['Component ' num2str(i)])
        caxis([0 2e3])
        count = count + 1;
    end

    components_movie(j) = getframe(gcf);

end

myVideo = VideoWriter('components_movie_onlygood.avi');
uncompressedVideo = VideoWriter('components_movie.avi', 'Uncompressed AVI');
myVideo.FrameRate = 5;  % Default 30

open(myVideo);
writeVideo(myVideo,components_movie);
close(myVideo);
%%

figure
imagesc(corr(spatial))
axis off
colorbar
title('Correlations between spatial components')

%%

figure
for i = 1:50
    
    subplot(50,1,i)
    plot(temporal(i,:))
    axis off
    
end

%%
figure
imagesc(temporal)
colormap gray
set(gca,'xticklabels',{})
title('Fluor Traces For Each Component')

%%
num_results = length(results)/5;

figure;
for i = 1:num_results
    
    subplot(num_results,1,i)
    all_points = [];
    for j = 1:5
        all_points = [all_points results((i-1)*5+j).trials.times];
    end
%     scatter(results(i).trials.times,results(i).trials.tau2)
    histogram(all_points,[0:1:2000])
    xlim([0 2000])
%     ylim([0 2000])
end


%%
latency_image1(10+8,4+11) = 0;
figure; subplot(131); 
imagesc(latency_image2(9:21,12:20)); axis off; 

subplot(132); imagesc(latency_image1(9:21,12:20)); axis off; 

subplot(133); axis off; imagesc(min(latency_image1(9:21,12:20),latency_image2(9:21,12:20)))
axis off
colormap hot

%%
sum_image1(10+8,4+11) = 0;
figure; subplot(131); 
imagesc(sum_image2(9:21,12:20)); axis off; 
caxis([0 1.5])
subplot(132); imagesc(sum_image1(9:21,12:20)); axis off; 
caxis([0 1.5])
subplot(133); axis off; imagesc(sum_image1(9:21,12:20) + sum_image2(9:21,12:20))
axis off
caxis([0 1.5])
colormap hot


%%


dates = {'10_7', '10_12', '10_12','10_7','10_7','10_7','10_8','10_8','10_8' , '10_8' , '10_12'}; %

slice_nums = [1 4 4 3 3 3 3 3 5 3 2];

cell_nums =[1 1 1 1 1 1 1 1 1 1 1]; 

tags = {'', '', '', '', '', '', 'cont', 'cont', '', '', ''};% 

trials = {[6],[1 2 3 4],[3 4 5 6 7 8],[1 4 5 6],[1 2 3],[4 5] ,[1 2]};% third from last
tracedir = '/media/shababo/Layover/projects/mapping/data/som-mapping-st/som-big-job-01-dist';

paramdir = '/media/shababo/Layover/projects/mapping/data/som-mapping-st/som-big-job-01-distparams';

%%
% mkdir(tracedir)
% mkdir(paramdir)

% load('/media/shababo/Layover/projects/mapping/data/som-mapping-st/map_index.mat')

% build_many_jobs(dates(6),slice_nums(6),cell_nums(6),tags(6),trials(6),params,map_index,tracedir,paramdir)

%% som mapping psc pairs

dates = {'10_8', '10_8', '10_12', '10_12', '10_12', '10_8', '10_7', '10_7', '10_7', '10_7',...
    '10_31','10_31','10_31','10_31','10_31','11_1','11_1','11_1','11_1','11_1','11_1','11_1','11_1',...
    '11_2','11_2','11_2','11_2','11_2','11_2','10_12'}; %,'10_7','10_7','10_7','10_8','10_8','10_8' , '10_8' , '10_12'}; %

slice_nums = [5 3 4 4 2 3 1 3 3 3 ...
    4 5 5 6 7 1 1 3 4 5 6 7 8 1 1 2 3 3 4 3]; %3 3 3 3 3 5 3 2];

cell_nums = [1 1 1 1 1 1 1 1 1 1 ...
    1 1 2 1 1 1 2 1 1 2 1 1 1 2 2 1 1 1 1 1]; % 1 1 1 1 1 1 1 1]; 

tags = {'', '', '', '', '', 'cont', '', '', '', '',...
    '','','','','','','','','','','','','','','','','next','next','', ''}; %, '', '', 'cont', 'cont', '', '', ''};% 

trials = {[1 2 3],[4 5],[1 2],[3 4],[1 2],[4 5 6],[6],[3 4],[5 6],[7 8],...
    [15 16],[9 10],[1 2 3],[6 7 8],[6 7],[6 7],[5 6],[9 10 11],[7 8],[5 6],[5 6],[5 6],[5 6],...
    [5 6],[ 7 8],[5 6],[5 6],[ 7 8],[9 10],[1 2]}; %,[1 4 5 6],[1 2 3],[4 5] ,[1 2]};% third from last
trials_detection = {[1:9],[1:6],[1:6],[7:12],1:6,4:12,1:3,1:6,7:12,13:18,...
    1:6,1:6,1:9,1:9,1:6,1:6,1:6,1:9,1:6,1:6,1:6,1:6,1:6,...
    1:6,7:12,1:6,1:6,7:12,1:6,1:6}; %,[1:6], [7:12], [12:18],[1:3], [4:12],[1:9],[1:6] ,[1:6]};

pair_id = logical([0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 1 0 0 0 1 1 1 1 1 1 0]);

map_index_id = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];

id_names = {'L4-L5','L5-L5'};

skip_exps = [2 3 7 9 10 24 28];
pair_id_tmp = pair_id;
pair_id_tmp(skip_exps) = [];

exps_to_run = setdiff(9,skip_exps); %  13 26

%%

% exps_to_run = [14:25 27:30];
for i = 1:length(exps_to_run)
    this_exp = exps_to_run(i);
    [glm_out_map_based(this_exp),~,~] = ...
        glm_xcorr_all_pairs('/media/shababo/data/old-data',dates(this_exp),slice_nums(this_exp),...
        cell_nums(this_exp),tags(this_exp),trials_detection(this_exp));
end


%%

% exps_to_run = [1:5 7:12 14:17 19:25 27:29];
% exps_to_run = [16:23];
exps_to_run = 14;
% 
close all

glm_to_plot = glm_out_map_based;

num_experiments = length(exps_to_run);
for ii = 1:length(exps_to_run)
    
    this_exp = exps_to_run(ii);
    
%     this_exp = 7;

       
    figure;
    subplot(121)
%     imagesc(reshape(...
%         glm_out_new(this_exp).ch1.glmnet_fit.beta(2:end,...
%         ceil((find(glm_out_new(this_exp).ch1.lambda == ...
%         glm_out_new(this_exp).ch1.lambda_min) + ...
%         find(glm_out_new(this_exp).ch1.lambda == ...
%         glm_out_new(this_exp).ch1.lambda_1se))/2)),21,21)'>0); colormap gray
%     lambda_ind = ceil(find(glm_to_plot(this_exp).ch1.glmnet_fit.beta(:) > 0,1,'first')/...
%         size(glm_to_plot(this_exp).ch1.glmnet_fit.beta,1)+4);
%     lambda_ind = max(find(glm_to_plot(this_exp).ch1.lambda == ...
%         glm_to_plot(this_exp).ch1.lambda_1se),lambda_ind);
%     lambda_ind = ceil(find(glm_to_plot(this_exp).ch1.lambda == ...
%         glm_to_plot(this_exp).ch1.lambda_min) + find(glm_to_plot(this_exp).ch1.lambda == ...
%         glm_to_plot(this_exp).ch1.lambda_1se)/2);
%     min_1seup = glm_to_plot(this_exp).ch1.cvup(glm_to_plot(this_exp).ch1.lambda == ...
%         glm_to_plot(this_exp).ch1.lambda_min);
%     lambda_ind = find(glm_to_plot(this_exp).ch1.cvm < min_1seup,1,'last');
    num_nonzeros = unique(sum(glm_to_plot(this_exp).ch1.glmnet_fit.beta > 0));
    target_num = num_nonzeros(end-1);
    lambda_ind = find(sum(glm_to_plot(this_exp).ch1.glmnet_fit.beta > 0) == target_num,1,'last');
    ch1_glm_map = reshape(...
        glm_to_plot(this_exp).ch1.glmnet_fit.beta(2:end,lambda_ind),21,21)';
    ch1_glm_map(ch1_glm_map < .25) = 0;
    imagesc(ch1_glm_map); colormap gray
    title(['Experiment ' num2str(this_exp) ': CH1'])
    axis off
    subplot(122)
%     imagesc(reshape(...
%         glm_out_new(this_exp).ch2.glmnet_fit.beta(2:end,...
%         ceil((find(glm_out_new(this_exp).ch2.lambda == ...
%         glm_out_new(this_exp).ch2.lambda_min) + ...
%         find(glm_out_new(this_exp).ch2.lambda == ...
%         glm_out_new(this_exp).ch2.lambda_1se))/2)),21,21)'>0); colormap gray
%     lambda_ind = ceil(find(glm_to_plot(this_exp).ch2.glmnet_fit.beta(:) > 0,1,'first')/...
%         size(glm_to_plot(this_exp).ch2.glmnet_fit.beta,1)+4);
%     lambda_ind = max(find(glm_to_plot(this_exp).ch2.lambda == ...
%         glm_to_plot(this_exp).ch2.lambda_1se),lambda_ind);
%     lambda_ind = ceil(find(glm_to_plot(this_exp).ch2.lambda == ...
%         glm_to_plot(this_exp).ch2.lambda_min) + find(glm_to_plot(this_exp).ch2.lambda == ...
%         glm_to_plot(this_exp).ch2.lambda_1se)/2);
%     min_1seup = glm_to_plot(this_exp).ch2.cvup(glm_to_plot(this_exp).ch2.lambda == ...
%         glm_to_plot(this_exp).ch2.lambda_min);
%     lambda_ind = find(glm_to_plot(this_exp).ch2.cvm < min_1seup,1,'last');
    num_nonzeros = unique(sum(glm_to_plot(this_exp).ch2.glmnet_fit.beta > 0));
    target_num = num_nonzeros(end-1);
    lambda_ind = find(sum(glm_to_plot(this_exp).ch2.glmnet_fit.beta > 0) == target_num,1,'last');
    ch2_glm_map = reshape(...
        glm_to_plot(this_exp).ch2.glmnet_fit.beta(2:end,lambda_ind),21,21)';
    ch2_glm_map(ch2_glm_map < .25) = 0;
    imagesc(ch2_glm_map); colormap gray
    title(['Experiment ' num2str(this_exp) ': CH2'])
    axis off

%     subplot(133)
%     imagesc(xcorr_images_new{this_exp})
%     title(['Experiment ' num2str(this_exp)])
    
     data_file = ...
        [dates{this_exp} '_slice' num2str(slice_nums(this_exp)) '_cell' ...
        num2str(cell_nums(this_exp)) '' tags{this_exp} '.mat'];
    load(data_file)
    [map_ch1,map_ch2] = see_grid(data,trials{this_exp},map_indices{map_index_id(this_exp)},1);
    title(['Experiment ' num2str(this_exp)])
%     figure;
%     compare_trace_stack_grid_overlap({map_ch1,samples_psths_new{this_exp}{1}},Inf,1,[],0,{'L4','L5'},1)
%     title(['Experiment ' num2str(this_exp) ': Detection CH1'])
% 
%     figure;
%     %     plot_trace_stack_grid(samples_psths{this_exp}{1},Inf,1,1);
%     imagesc(cellfun(@(x) max(mean(x(:,150:450),1)),samples_psths_new{this_exp}{1}))
%     title(['Experiment ' num2str(this_exp) ': Max(mean(PSTH)) CH1'])
    
%     figure;
%     compare_trace_stack_grid_overlap({map_ch2,samples_psths_new{this_exp}{2}},Inf,1,[],0,{'L4','L5'},1)
%     title(['Experiment ' num2str(this_exp) ': Detection CH2'])
    
%     figure;
%     %     plot_trace_stack_grid(samples_psths{this_exp}{2},Inf,1,1);
%     imagesc(cellfun(@(x) max(mean(x(:,150:450),1)),samples_psths_new{this_exp}{2}))
%     title(['Experiment ' num2str(this_exp) ': Max(mean(PSTH)) CH2'])
    
%     figure;
%     compare_trace_stack_grid_overlap({samples_psths_new{this_exp}{1},samples_psths_new{this_exp}{2}},Inf,1,[],0,{'L4','L5'},1)
%     title(['Experiment ' num2str(this_exp) ': Detection CH1 vs. CH2'])
    
%     xcorr_tighter{this_exp} = cellfun(@(x,y) xcorr_peak_trials(x,y,[150 450],0),...
%         samples_psths_new{this_exp}{1},samples_psths_new{this_exp}{2});
%     xcorr_images_null{this_exp} = cellfun(@(x,y) xcorr_peak_trials(x,y,[150 450],1),...
%         samples_psths_new{this_exp}{1},samples_psths_new{this_exp}{2});
%     
%     figure; subplot(121); imagesc(xcorr_tighter{this_exp}); subplot(122); imagesc(xcorr_images_null{this_exp})
    
end

%%

figure;
subplot(121)
imagesc(reshape(...
    glm_out_l5.glmnet_fit.beta(2:end,find(glm_out_l5.lambda == ...
    glm_out_l5.lambda_1se)),21,21)'>0 ); colormap gray
title(['Experiment ' num2str(this_exp) ': CH1'])
axis off

subplot(122)
imagesc(reshape(...
    glm_out_l4.glmnet_fit.beta(2:end,find(glm_out_l4.lambda == ...
    glm_out_l4.lambda_min)+5),21,21)'>0); colormap gray
title(['Experiment ' num2str(this_exp) ': CH2'])
axis off

%%

[glm_out_test,xcorr_images_test,samples_psths_test] = ...
        glm_xcorr_all_pairs('/media/shababo/data/old-data',dates(this_exp),slice_nums(this_exp),...
        cell_nums(this_exp),tags(this_exp),trials_detection(this_exp));

%%

exps_to_run = [1:5 7:8 11:12 14:17 19:25 27:30];
% exps_to_run = 24;

for ii = 1:length(exps_to_run)

    this_exp = exps_to_run(ii)
    
    events1 = samples_psths_new{this_exp}{1};
    events2 = samples_psths_new{this_exp}{2};

    disp('got events')
    input_locs1 = find(glm_out_newer(this_exp).ch1.glmnet_fit.beta(2:end,find(glm_out_newer(this_exp).ch1.lambda == ...
            glm_out_newer(this_exp).ch1.lambda_1se)+2));%+5*(-pair_id(this_exp) + 1)
    input_locs2 = find(glm_out_newer(this_exp).ch2.glmnet_fit.beta(2:end,find(glm_out_newer(this_exp).ch2.lambda == ...
        glm_out_newer(this_exp).ch2.lambda_1se)+2+5*(-pair_id(this_exp)+1)));

%     all_locs = union(input_locs1,input_locs2);
%     all_locs = intersect(input_locs1,input_locs2);
    all_locs = 1:441;
    [js,is] = ind2sub([21 21],all_locs);
    
    null = 1;
    if null
        iters = 10000;
        for j = 1:iters
            input_locs1 =  randsample(1:441,length(input_locs1));

            num_event_locs_null(this_exp,j) = length(union(input_locs1,input_locs2));
            num_shared_locs_null(this_exp,j) = length(intersect(input_locs1,input_locs2));
        end
    else
        num_event_locs = length(union(input_locs1,input_locs2));
        num_shared_locs = length(intersect(input_locs1,input_locs2));
    end

    disp('got locs')

    if ~isempty(all_locs)
        psth1 = zeros(1,length(events1{1,1}(1,:)));
        psth2 = zeros(size(psth1));

        total_trials1 = 0;
        total_trials2 = 0;

        norm_factor = 1;
        raw_jpsth = zeros(length(events1{1,1}(1,:)));
        raw_psth1 = zeros(length(events1{1,1}(1,:)));
        raw_psth2 = zeros(length(events1{1,1}(1,:)));
        disp('getting psths')
        disp('for all locs...')
        length(all_locs)
%         all_locs_save(this_exp) = length(all_locs);
        for i = 1:length(all_locs)
    %         i
    %         psth1 = psth1 + sum(events1{is(i),js(i)},1)/norm_factor;
    %         psth2 = psth2 + sum(events2{is(i),js(i)},1)/norm_factor;

            for j= 1:size(events1{is(i),js(i)},1)
    %             raw_jpsth = raw_jpsth + events1{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
    %             raw_psth1 = raw_psth1 + events1{is(i),js(i)}(j,:)'/norm_factor*events1{is(i),js(i)}(j,:)/norm_factor;
    %             raw_psth2 = raw_psth2 + events2{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
                [~, event_times] = findpeaks(events1{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events1{is(i),js(i)}(j,:)),'MinPeakDistance',20);
                events1{is(i),js(i)}(j,:) = zeros(size(events2{is(i),js(i)}(j,:)));
                events1{is(i),js(i)}(j,event_times) = 1;

                [~, event_times] = findpeaks(events2{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events2{is(i),js(i)}(j,:)),'MinPeakDistance',20);
                events2{is(i),js(i)}(j,:) = zeros(size(events2{is(i),js(i)}(j,:)));
                events2{is(i),js(i)}(j,event_times) = 1;
            end
            events_tmp1 = events1{is(i),js(i)};
            events_tmp2 = events2{is(i),js(i)};


            psth1 = psth1 + sum(events1{is(i),js(i)},1)/norm_factor;
            psth2 = psth2 + sum(events2{is(i),js(i)},1)/norm_factor;

            raw_jpsth = raw_jpsth + events1{is(i),js(i)}'/norm_factor*events2{is(i),js(i)}/norm_factor;
            raw_psth1 = raw_psth1 + events1{is(i),js(i)}'/norm_factor*events1{is(i),js(i)}/norm_factor;
            raw_psth2 = raw_psth2 + events2{is(i),js(i)}'/norm_factor*events2{is(i),js(i)}/norm_factor;
            total_trials1 = total_trials1 + size(events1{is(i),js(i)},1);
            total_trials2 = total_trials2 + size(events2{is(i),js(i)},1);

        end

        psth1 = psth1/total_trials1;
        psth2 = psth2/total_trials2;

        jpsth = psth1'*psth2;
        jpsth1 = psth1'*psth1;
        jpsth2 = psth2'*psth2;
        raw_jpsth = raw_jpsth/total_trials1;
        raw_psth1 = raw_psth1/total_trials1;
        raw_psth2 = raw_psth2/total_trials1;

    %     eff_mat = sqrt(raw_psth1.*(1-raw_psth1))./sqrt(raw_psth2.*(1-raw_psth2));
    %     con_mat = sqrt(raw_psth2.*(1-raw_psth2))./sqrt(raw_psth1.*(1-raw_psth1));

        disp('got psths')
        disp('plotting...')
    %     figure
    %     subplot(221)
    %     imagesc(jpsth)
    %     xlim([150 800])
    %     ylim([150 800])
    %     title(['Experiment ' num2str(this_exp) ': ' num2str(trace(jpsth))])

    %     subplot(222)
    %     imagesc(raw_jpsth)
    %     xlim([150 800])
    %     ylim([150 800])
    %     title(trace(raw_jpsth))

        D = raw_jpsth - jpsth;
        D1 = raw_psth1 - jpsth1;
        D2 = raw_psth2 - jpsth2;
        D1(D1 < 0) = 0;
        D2(D2 < 0) = 0;
    %     subplot(223)
    %     imagesc(D)
    %     xlim([150 800])
    %     ylim([150 800])
    %     title(trace(D))

        % diagD = diag(D);
        % diagD(diagD < 0) = 0;

        normalization_matrix = sqrt(diag(D1)*diag(D2)');
        D_norm = D./normalization_matrix;
    %     D_norm(isnan(D_norm)) = 0;
    %     D_norm(D_norm < 0) = 0;
    %     D_norm = sqrt(D_norm);
        D_norm_all_alllocs(:,:,this_exp) = D_norm;
    %     E = D_norm*eff_mat;
    %     C = D_norm*con_mat;
    %     E(isnan(E)) = 0;
    %     C(isnan(C)) = 0;

%         common_input_score_norm_intersect(this_exp) = 0;
%         offsets = 0;
%         for k = offsets
%             D_norm(D_norm < 0) = 0;
%             this_diag = diag(sqrt(D_norm),k);
% 
%             common_input_score_norm_intersect(this_exp) = common_input_score_norm_intersect(this_exp) + nanmean(this_diag(150-floor(abs(k)/2):450-floor(abs(k)/2)));%/(301*length(offsets))
%         end
    %     subplot(224)
    %     imagesc(D_norm)
    %     xlim([150 800])
    %     ylim([150 800])
    %     title(common_input_score_norm(this_exp))
    else
%         common_input_score_norm_intersect(this_exp) = 0;
        D_norm_all_alllocs(:,:,this_exp) = zeros(1999);
    end
%     disp('result:')
%     this_exp
%     pair_id(this_exp)
%     common_input_score_norm_intersect(this_exp)
    
    
end

%%


exps_to_run = [1:5 7:8 11:12 14:17 19:25 27:30];
% exps_to_run = [16:17 19:25 27:30];
exps_to_run = 18;

for ii = 1:length(exps_to_run)

    this_exp = exps_to_run(ii)
    
    events1 = samples_psths_new{this_exp}{1};
    events2 = samples_psths_new{this_exp}{2};

    disp('got events')
%     input_locs1 = find(glm_out_newer(this_exp).ch1.glmnet_fit.beta(2:end,find(glm_out_newer(this_exp).ch1.lambda == ...
%             glm_out_newer(this_exp).ch1.lambda_1se)+2));%+5*(-pair_id(this_exp) + 1)
%     input_locs2 = find(glm_out_newer(this_exp).ch2.glmnet_fit.beta(2:end,find(glm_out_newer(this_exp).ch2.lambda == ...
%         glm_out_newer(this_exp).ch2.lambda_1se)+2+5*(-pair_id(this_exp)+1)));
% 
%     all_locs = union(input_locs1,input_locs2);
%     all_locs = intersect(input_locs1,input_locs2);
    all_locs = 1:441;
%     all_locs = all_locs(26)
    [js,is] = ind2sub([21 21],all_locs);
    
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


    iters = 500;
    if ~isempty(all_locs)



        length(all_locs)
%         all_locs_save(this_exp) = length(all_locs);
        for i = 1:length(all_locs)
            if mod(i,10) == 0
                i
            end
    %         psth1 = psth1 + sum(events1{is(i),js(i)},1)/norm_factor;
    %         psth2 = psth2 + sum(events2{is(i),js(i)},1)/norm_factor;
            exp_shuffle_stats(this_exp).locs(i).null_dist = zeros(iters,1);
            exp_shuffle_stats(this_exp).locs(i).loc_ind = [is(i) js(i)];
            
            this_loc_events1 = zeros(size(events1{is(i),js(i)},1),65)';
            this_loc_events2 = zeros(size(events1{is(i),js(i)},1),65)';

            for j= 1:size(events1{is(i),js(i)},1)
    %             raw_jpsth = raw_jpsth + events1{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
    %             raw_psth1 = raw_psth1 + events1{is(i),js(i)}(j,:)'/norm_factor*events1{is(i),js(i)}(j,:)/norm_factor;
    %             raw_psth2 = raw_psth2 + events2{is(i),js(i)}(j,:)'/norm_factor*events2{is(i),js(i)}(j,:)/norm_factor;
                [~, event_times] = findpeaks(events1{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events1{is(i),js(i)}(j,20:end)),'MinPeakDistance',10);
                this_loc_events1(ceil((event_times-1)/10),j) = 1;

                [~, event_times] = findpeaks(events2{is(i),js(i)}(j,150:800),'MinPeakHeight',0.25*std(events2{is(i),js(i)}(j,20:end)),'MinPeakDistance',10);
                this_loc_events2(ceil((event_times-1)/10),j) = 1;
            end
            
            concat1_data = this_loc_events1(:);
            concat2_data = this_loc_events2(:);
            exp_shuffle_stats(this_exp).locs(i).zero_lag_corr = concat1_data(1:end-1)'*concat2_data(1:end-1)+ ...
                    concat1_data(2:end)'*concat2_data(1:end-1) + ...
                    concat1_data(1:end-1)'*concat2_data(2:end);
            
            tmp = repmat(1:ceil(length(concat1_data)/12),12,1);
            tmp = tmp(:);
            tmp = tmp(1:length(concat1_data));
            poisson_rate1 = accumarray(tmp(:),concat1_data);
            poisson_rate2 = accumarray(tmp(:),concat2_data);  
            
            for k = 1:iters
                event_bins1 = find(poisson_rate1);
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
                exp_shuffle_stats(this_exp).locs(i).null_dist(k) = shuffled_events1(1:end-1)'*shuffled_events2(1:end-1) + ...
                    shuffled_events1(2:end)'*shuffled_events2(1:end-1) + ...
                    shuffled_events1(1:end-1)'*shuffled_events2(2:end);
            end
            [f,x] = ecdf(exp_shuffle_stats(this_exp).locs(i).null_dist);
            idx = find(x < exp_shuffle_stats(this_exp).locs(i).zero_lag_corr,1,'last');
            if ~isempty(idx)
                exp_shuffle_stats(this_exp).locs(i).p_val = 1 - f(idx);
            else
                exp_shuffle_stats(this_exp).locs(i).p_val = 1;
            end
               
        end

        exp_shuffle_stats(this_exp).input_map = ones(21,21);
        for i = 1:length(exp_shuffle_stats(this_exp).locs)
            inds = exp_shuffle_stats(this_exp).locs(i).loc_ind;
            exp_shuffle_stats(this_exp).input_map(inds(1),inds(2)) = exp_shuffle_stats(this_exp).locs(i).p_val;
        end
        
        figure;
        imagesc(exp_shuffle_stats(this_exp).input_map < .05)
        title(['Experiment: ' num2str(this_exp)])
       
    else
%         common_input_score_norm_intersect(this_exp) = 0;
        exp_shuffle_stats(this_exp).locs(i).null_dist = [];
        exp_shuffle_stats(this_exp).locs(i).loc_ind = [];
        exp_shuffle_stats(this_exp).locs(i).zero_lag_corr = 0;
    end

    
end


%%

exps_to_run = [1:5 7:8 11:12 14:17 19:25 27:30];
% exps_to_run = 30;

for ii = 1:length(exps_to_run)
    this_exp = exps_to_run(ii)
%     common_input_score_alllocs(this_exp) = 0;
    offsets = -40:40;
    for k = offsets
        this_D = D_norm_all_intersect(:,:,this_exp);
        this_D(this_D < 0) = 0;
        this_diag = diag(sqrt(this_D),k);
        norm_xcorrs_intersect(this_exp,k-offsets(1)+1) = nanmean(this_diag(150-floor(abs(k)/2):450-floor(abs(k)/2)));
         norm_xcorrs_intersect(isnan(norm_xcorrs)) = 0;
%           common_input_score_alllocs(this_exp) = nanmean(norm_xcorrs(this_exp,30:50));
%         common_input_score_norm_intersect(this_exp) = common_input_score_norm_intersect(this_exp) + nanmean(this_diag(150-floor(abs(k)/2):450-floor(abs(k)/2)));%/(301*length(offsets))
    end
end

%%
for ii = 1:length(exps_to_run)
    this_exp = exps_to_run(ii)
%         figure;
%         imagesc(exp_shuffle_stats(this_exp).input_map < .05/(21*441)); caxis([0 1])
%         colorbar
%         title(['Experiment: ' num2str(this_exp) ', Type: ' id_names{pair_id(this_exp) + 1}])
    num_common_locs(this_exp) = sum(sum( exp_shuffle_stats(this_exp).input_map < .05/(21*441)));
end
%%

for i = 1:size(norm_xcorrs,1)
   
    common_input_score_norm(i) = nanmean(norm_xcorrs(i,30:50))
end
%%
cis_tmp_null = num_shared_locs_null./num_event_locs_null;
cis_tmp_null(skip_exps,:) = [];

cis_tmp = num_shared_locs./num_event_locs;
cis_tmp(skip_exps) = [];

figure;
subplot(121)
null_dist = cis_tmp_null(~pair_id_tmp,:);
histogram(null_dist(:),0:.01:.6);
hold on
counts = histcounts(null_dist(:),0:.01:.6);
scatter(cis_tmp(~pair_id_tmp),max(counts)/2*ones(sum(~pair_id_tmp),1))
cdf95 = quantile(null_dist(:),1 - .05/length(pair_id_tmp));
plot([cdf95 cdf95],[0 max(counts)])
title('L4-L5 Pairs, Spatial Coincidence')
legend({'Shuffled Dist','Data','p = .05 corrected'})
ylabel('Counts')
xlabel('Coincidence Probability')

subplot(122)
null_dist = cis_tmp_null(pair_id_tmp,:);
histogram(null_dist(:),0:.01:.6);
hold on
counts = histcounts(null_dist(:),0:.01:.6);
scatter(cis_tmp(pair_id_tmp),max(counts)/2*ones(sum(pair_id_tmp),1))
cdf95 = quantile(null_dist(:),1 - .05/length(pair_id_tmp));
plot([cdf95 cdf95],[0 max(counts)])
title('L5-L5 Pairs, Spatial Coincidence')
legend({'Shuffled Dist','Data','p = .05 corrected'})
ylabel('Counts')
xlabel('Coincidence Probability')

%%

bin_mat = [zeros(1,2) ones(1,5) zeros(1,74)
           zeros(1,7) ones(1,5) zeros(1,69)
           zeros(1,12) ones(1,5) zeros(1,64)
           zeros(1,17) ones(1,5) zeros(1,59)
           zeros(1,22) ones(1,5) zeros(1,54)So now
           zeros(1,27) ones(1,5) zeros(1,49)
           zeros(1,32) ones(1,5) zeros(1,44)
           zeros(1,37) ones(1,5) zeros(1,39)
           zeros(1,42) ones(1,5) zeros(1,34)
           zeros(1,47) ones(1,5) zeros(1,29)
           zeros(1,52) ones(1,5) zeros(1,24)
           zeros(1,57) ones(1,5) zeros(1,19)
           zeros(1,62) ones(1,5) zeros(1,14)
           zeros(1,67) ones(1,5) zeros(1,9)
           zeros(1,72) ones(1,5) zeros(1,4)
           ];
binned_xcorr = bin_mat * norm_xcorrs_intersect';
% 
% for i = 1:30
%     figure;
%     stairs([binned_xcorr(:,i); binned_xcorr(end,i)])
%     xlim([1 16])
%     set(gca,'xtick',[1.5 8.5 15.5])
%     set(gca,'xticklabels',{'-1.75','0.0','1.75'})
% end
% figure;
% subplot(131)
% imagesc(sqrt(D_norm))
% subplot(132)
% imagesc(E)
% subplot(133)
% imagesc(C)





%%
binned_xcorr(isnan(binned_xcorr)) = 0;
cis_tmp = num_common_locs;
cis_tmp(skip_exps) = [];

figure
% binned_xcorr(:,skip_exps) = [];
% 
% figure;
% subplot(211)
% % good_exps = setdiff(1:30,skip_exps);
% % pair_id_good = pair_id(good_exps);
% mean_xcorr_1 = mean(binned_xcorr(:,~pair_id_tmp),2);
% stairs([mean_xcorr_1(1) mean_xcorr_1'],'Linewidth',2)
% hold on
% mean_xcorr_2 = mean(binned_xcorr(:,pair_id_tmp),2);
% stairs([mean_xcorr_2(1) mean_xcorr_2'],'Linewidth',2) 
% hold off
% xlim([1 16])
% set(gca,'xtick',[1.5 8.5 15.5])
% set(gca,'xticklabels',{'-1.75','0.0','1.75'})
% title('Normalized Correlogram')
% xlabel('Lag (msec)')
% ylabel('Synchrony Probability')
% legend({'L4-L5 Pairs','L5-L5 Pairs'})

% subplot(2,1,2)
scatter(ones(sum(~pair_id_tmp),1),cis_tmp(~pair_id_tmp),'k' ,'jitter','on', 'jitterAmount',0.05);
hold on
scatter(1.25,mean(cis_tmp(~pair_id_tmp)),100,[0 0 0],'filled')
hold on
plot(1.25*[1 1],mean(cis_tmp(~pair_id_tmp)) + std(cis_tmp(~pair_id_tmp))/sqrt(length(cis_tmp(~pair_id_tmp)))*[-1 1],'k','Linewidth',2)
hold on
scatter(2*ones(sum(pair_id_tmp),1),cis_tmp(pair_id_tmp),'k','jitter','on', 'jitterAmount',0.05);
hold on
scatter(2.25,mean(cis_tmp(pair_id_tmp)),100,[0 0 0],'filled')
hold on
plot(2.25*[1 1],mean(cis_tmp(pair_id_tmp)) + std(cis_tmp(pair_id_tmp))/sqrt(length(cis_tmp(pair_id_tmp)))*[-1 1],'k','Linewidth',2)
hold on
xlim([0 3])
set(gca,'Xtick',[1 2])
set(gca,'xticklabels',{['L4-L5 Pairs (N = ' num2str(sum(~pair_id_tmp)) ')'],['L5-L5 Pairs (N = ' num2str(sum(pair_id_tmp)) ')']})
ylabel('Common Input Locations')

[p,h] = ranksum(cis_tmp(~pair_id_tmp),cis_tmp(pair_id_tmp),'tail','left');
title(['p = ' num2str(p)])

%%
cis_tmp_null = num_shared_locs_null./num_event_locs_null;
cis_tmp_null(skip_exps) = [];

cis_tmp = num_shared_locs./num_event_locs;
cis_tmp(skip_exps) = [];

figure

scatter(2.5*ones(sum(pair_id_tmp),1),cis_tmp(pair_id_tmp),'k');
hold on
scatter(2.5,mean(cis_tmp(pair_id_tmp)),100,[0 0 0],'filled')
hold on
plot(2.5*[1 1],mean(cis_tmp(pair_id_tmp)) + std(cis_tmp(pair_id_tmp))/sqrt(length(cis_tmp(pair_id_tmp)))*[-1 1],'k','Linewidth',2)
hold on

scatter(3*ones(sum(pair_id_tmp),1),cis_tmp_null(pair_id_tmp),'k');
hold on
scatter(3,mean(cis_tmp_null(pair_id_tmp)),100,[0 0 0],'filled')
hold on
plot([3 3],mean(cis_tmp_null(pair_id_tmp)) + std(cis_tmp_null(pair_id_tmp))/sqrt(length(cis_tmp_null(pair_id_tmp)))*[-1 1],'k','Linewidth',2)
hold on







scatter(1.0*ones(sum(~pair_id_tmp),1),cis_tmp(~pair_id_tmp),'k');
hold on
scatter(1,mean(cis_tmp(~pair_id_tmp)),100,[0 0 0],'filled')
hold on
plot(1*[1 1],mean(cis_tmp(~pair_id_tmp)) + std(cis_tmp(~pair_id_tmp))/sqrt(length(cis_tmp(~pair_id_tmp)))*[-1 1],'k','Linewidth',2)
hold on

scatter(1.5*ones(sum(~pair_id_tmp),1),cis_tmp_null(~pair_id_tmp),'k');
hold on
scatter(1.5,mean(cis_tmp_null(~pair_id_tmp)),100,[0 0 0],'filled')
hold on
plot([1.5 1.5],mean(cis_tmp_null(~pair_id_tmp)) + std(cis_tmp_null(~pair_id_tmp))/sqrt(length(cis_tmp_null(~pair_id_tmp)))*[-1 1],'k','Linewidth',2)
hold on

xlim([0 4])
set(gca,'Xtick',[1 1.5 2.5 3])
set(gca,'xticklabels',{['L4-L5 Pairs (N = ' num2str(sum(~pair_id_tmp)) ')'],'L4-L5 Shuffle',['L5-L5 Pairs (N = ' num2str(sum(pair_id_tmp)) ')'],'L5-L5 Shuffle'})
ylabel('Common Input Probability')

[pl4,h] = ranksum(cis_tmp(~pair_id_tmp),cis_tmp_null(~pair_id_tmp),'tail','right')
[pl5,h] = ranksum(cis_tmp(pair_id_tmp),cis_tmp_null(pair_id_tmp),'tail','right')
title(['p_l4 = ' num2str(pl4) ', p_l5 = ' num2str(pl5)])

% xlim([0 3])
% set(gca,'Xtick',[1 2])
% set(gca,'xticklabels',{['L4-L5 Pairs (N = ' num2str(sum(~pair_id_tmp)) ')'],['L5-L5 Pairs (N = ' num2str(sum(pair_id_tmp)) ')']})
% ylabel('Common Input Probability')
% 
% [p,h] = ranksum(cis_tmp(~pair_id_tmp),cis_tmp(pair_id_tmp),'tail','left');
% title(['p = ' num2str(p)])

%%



exps_to_run = 8;

num_experiments = length(exps_to_run);


for i = 1:num_experiments

    this_exp = exps_to_run(i);

    data_file = ...
        [dates{this_exp} '_slice' num2str(slice_nums(this_exp)) '_cell' ...
        num2str(cell_nums(this_exp)) '' tags{this_exp} '.mat'];
    load(data_file)
    [map_ch1,map_ch2] = see_grid(data,trials{this_exp},map_index,0);
    
    figure;
    compare_trace_stack_grid_overlap({map_ch1,samples_psths_new{this_exp}{1}},Inf,1,[],0,{'L4','L5'},1)

    figure;
    %     plot_trace_stack_grid(samples_psths{this_exp}{1},Inf,1,1);
    imagesc(cellfun(@(x) max(mean(x(:,150:450),1)),samples_psths_new{this_exp}{1}))
    
    figure;
    compare_trace_stack_grid_overlap({map_ch2,samples_psths_new{this_exp}{2}},Inf,1,[],0,{'L4','L5'},1)
    
    figure;
    %     plot_trace_stack_grid(samples_psths{this_exp}{2},Inf,1,1);
    imagesc(cellfun(@(x) max(mean(x(:,150:450),1)),samples_psths_new{this_exp}{2}))
    
    figure;
    compare_trace_stack_grid_overlap({samples_psths_new{this_exp}{1},samples_psths_new{this_exp}{2}},Inf,1,[],0,{'L4','L5'},1)

end


%%
%%


dates = {'10_31','10_31','10_31','10_31','10_31','11_1','11_1','11_1','11_1','11_1','11_1','11_1','11_1'};

slice_nums = [4,5,5,6,7,1,1,3,4,5,6,7,8];

cell_nums =[1,1,2,1,1,1,2,1,1,2,1,1,1]; 

tags = {'','','','','','','','','','','','',''};

trials = {[15 16],[9 10],[1 2 3],[6 7 8],[6 7],[6 7],[5 6],[9 10 11],[7 8],[5 6],[5 6],[5 6],[5 6]};% third from last

tracedir = '/media/shababo/Layover/projects/mapping/data/som-mapping-st/som-big-job-02-dist';

paramdir = '/media/shababo/Layover/projects/mapping/data/som-mapping-st/som-big-job-02-distparams';

mkdir(tracedir)
mkdir(paramdir)

% load('/media/shababo/Layover/projects/mapping/data/som-mapping-st/map_index.mat')

build_many_jobs(dates,slice_nums,cell_nums,tags,trials,params,map_index2,tracedir,paramdir)

%%

dates = {'11_2','11_2','11_2','11_2'};

slice_nums = [1,2,3,4];

cell_nums =[2,1,1,1,]; 

tags = {'','','next',''};

trials = {[5 6 7 8],[5 6],[5 6 7 8],[9 10]};% third from last

tracedir = '/media/shababo/Layover/projects/mapping/data/som-mapping-st/som-big-job-02-11_2-dist';

paramdir = '/media/shababo/Layover/projects/mapping/data/som-mapping-st/som-big-job-02-11_2-distparams';

mkdir(tracedir)
mkdir(paramdir)

% load('/media/shababo/Layover/projects/mapping/data/som-mapping-st/map_index.mat')

build_many_jobs(dates,slice_nums,cell_nums,tags,trials,params,map_index2,tracedir,paramdir)

%% som  spike detection

%% som spike maps

dates = {'11_3','11_3','11_3','11_3','11_3','11_3','11_3','11_3','11_3',...
    '11_4','11_4','11_4','11_4','11_4','11_4','11_4','11_4','11_4',...
    '11_4','11_4','11_4','11_4','11_4','11_4','11_4','11_4','11_4','11_4'}; %

slice_nums = [1 1 1 1 1 2 2 2 2 ...
    1 1 1 2 2 2 2 2 2 ...
    3 3 3 3 3 3 3 3 3 3];

cell_nums = [1 1 1 1 1 1 3 4 4 ...
    1 1 1 1 1 2 2 2 2 ...
    1 1 4 4 5 5 6 6 7 7]; 

tags = {'','','','','','next','','','',...
    '','', '', '', '','', '', '', '',...
    '','','','','','','','','',''}; %, '', '', 'cont', 'cont', '', '', ''};% 

trials = {17:19,20:22,24:25,26:28,29,1:2,1:3,1:3,4:6,...
    6:8,17:18,19:22,11:13,14:16,10:12,13:15,16:18,19:21,...
    5:7,8:10,13:15,16:18,21:23,24:26,5:7,8:10,9:11,12:14}; %,[1 4 5 6],[1 2 3],[4 5] ,[1 2]};% third from last
trials_detection = {}; %,[1:6], [7:12], [12:18],[1:3], [4:12],[1:9],[1:6] ,[1:6]};

% 1 cell-attached, 2 current clamp, 3 voltage clamp
patch_type = [2 2 3 2 3 1 1 1 1 ...
    2 3 2 1 1 1 1 1 1 ...
    1 1 1 1 1 1 1 1 1 1];

row_inds = {4:8,4:8,4:8,4:8,4:8,4:6,4:6,4:6,4:6,...
    4:6,4:6,4:6,1:11,1:11,1:5,1:5,3:7,3:7,...
    1:11,1:11,4:6,4:6,1:11,1:11,4:6,4:6,3:7,3:7};
col_inds = {4:8,4:8,4:8,4:8,4:8,5:7,5:7,5:7,5:7,...
    5:7,5:7,5:7,1:11,1:11,1:5,1:5,4:8,4:8,...
    1:11,1:11,5:7,5:7,1:11,1:11,5:7,5:7,4:8,4:8};

% 1 is 200 mW @ 4 msec and 2 is 250 mW @ 10 msec
power = [2 2 2 1 1 2 2 2 1 ...
    2 2 2 2 1 2 1 1 2 ...
    2 1 2 1 1 2 2 1 1 2];


% id_names = {'L4-L5','L5-L5'};
% 
skip_exps = [];
% pair_id_tmp = pair_id;
% pair_id_tmp(skip_exps) = [];
% 
exps_to_run = 1:length(dates); %  13 26
exps_to_run = setdiff(exps_to_run,skip_exps)
%%
close all
%issues 13,19,20 7 14 15
% 16 too low
exps_to_run = [ 16];
detection_grids = cell(1,length(dates));

thresholds = [40 50 0 450 0 6.5 7.5 15 15 30 0 30 5 6 6 8 10 ...
    10 5 5 15 15 10 5 15 15 15 15];
for ii = 1:length(exps_to_run)
    
    this_exp = exps_to_run(ii);
    data_file = ...
        [dates{this_exp} '_slice' num2str(slice_nums(this_exp)) '_cell' ...
        num2str(cell_nums(this_exp)) '' tags{this_exp} '.mat'];
    load(data_file)
    figure
    subplot(121)
    [map_ch1,~] = see_grid(data,trials{this_exp},map_index_spikes,0,363);
    plot_trace_stack_grid(map_ch1,Inf,1,0);
    title(['Experiment ' num2str(this_exp)])
    
    this_row_i = row_inds{this_exp};
    this_col_j = col_inds{this_exp};
    detection_results = cell(length(this_row_i),length(this_col_j));
    
    for i = 1:length(this_row_i)
        for j = 1:length(this_col_j)
            switch patch_type(this_exp)
                case 1
                    detection_results{i,j} = detect_peaks(...
...%                         -1.0*bsxfun(@minus,map_ch1{i,j},median(map_ch1{i,j}(:,1:100),2)),5,20,1,1,0,1)*70;
                            -1.0*bsxfun(@minus,map_ch1{this_row_i(i),this_col_j(j)},map_ch1{this_row_i(i),this_col_j(j)}(:,1)),thresholds(this_exp),80,1,1,0,0,1)*70;
                case 2
                    detection_results{i,j} = detect_peaks(...
                        1.0*bsxfun(@minus,map_ch1{this_row_i(i),this_col_j(j)},map_ch1{this_row_i(i),this_col_j(j)}(:,1)),thresholds(this_exp),80,1,1,0,0,0)*70;
                case 3
                    detection_results{i,j} = mean(map_ch1{this_row_i(i),this_col_j(j)});
            end
        end
    end
    
    detection_grids{this_exp} = detection_results;
    subplot(122)
    plot_trace_stack_grid(detection_results,Inf,1,0);
end
    

