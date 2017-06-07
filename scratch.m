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
    '11_2','11_2','11_2','11_2','11_2','11_2','10_12','10_12'}; %,'10_7','10_7','10_7','10_8','10_8','10_8' , '10_8' , '10_12'}; %

slice_nums = [5 3 4 4 2 3 1 3 3 3 ...
    4 5 5 6 7 1 1 3 4 5 6 7 8 1 1 2 3 3 4 3 5]; %3 3 3 3 3 5 3 2];

cell_nums = [1 1 1 1 1 1 1 1 1 1 ...
    1 1 2 1 1 1 2 1 1 2 1 1 1 2 2 1 1 1 1 1 1]; % 1 1 1 1 1 1 1 1]; 

tags = {'', '', '', '', '', 'cont', '', '', '', '',...
    '','','','','','','','','','','','','','','','','next','next','', '',''}; %, '', '', 'cont', 'cont', '', '', ''};% 

trials = {[1 2 3],[4 5],[1 2],[3 4],[1 2],[4 5 6],[6],[3 4],[5 6],[7 8],...
    [15 16],[9 10],[1 2 3],[6 7 8],[6 7],[6 7],[5 6],[9 10 11],[7 8],[5 6],[5 6],[5 6],[5 6],...
    [5 6],[ 7 8],[5 6],[5 6],[ 7 8],[9 10],[1 2],[1 2]}; %,[1 4 5 6],[1 2 3],[4 5] ,[1 2]};% third from last
trials_detection = {[1:9],[1:6],[1:6],[7:12],1:6,4:12,1:3,1:6,7:12,13:18,...
    1:6,1:6,1:9,1:9,1:6,1:6,1:6,1:9,1:6,1:6,1:6,1:6,1:6,...
    1:6,7:12,1:6,1:6,7:12,1:6,1:6,1:6}; %,[1:6], [7:12], [12:18],[1:3], [4:12],[1:9],[1:6] ,[1:6]};

pair_id = logical([0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 1 0 0 0 1 1 1 1 1 1 0 1]);

map_index_id = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1];

id_names = {'L4-L5','L5-L5'};

skip_exps = [2 3 7 9 10 24 28];
pair_id_tmp = pair_id;
pair_id_tmp(skip_exps) = [];

exps_to_run = setdiff(1:length(pair_id),skip_exps); %  13 26

%%

exps_to_run = [1:31];
% exps_to_run = [20 28];
% delete(gcp('nocreate'))
% pool = parpool(3);
% poisson_prob_maps = [];


% for i = 1:length(exps_to_run)
%     this_exp = exps_to_run(i)
%     try
%      [~,~,time_posteriors_2017{this_exp},new_event_detection_2017{this_exp}] = ...
% ...%         [glm_out_vb_good_onetimebin3{this_exp}] = ...
%            glm_xcorr_all_pairs('/media/shababo/data/new-data',dates(this_exp),slice_nums(this_exp),...
%            cell_nums(this_exp),tags(this_exp),trials_detection(this_exp),glm_dummy);
% ...%             glm_from_map_est(map_estimates_full(this_exp),glm_dummy);
% ...                ROI_VB_som_data(map_estimates(this_exp),[1 pair_id(this_exp)]);
%     catch e
%         disp([num2str(this_exp) ' fail'])
%     end
% end

for i = 1:length(exps_to_run)
    this_exp = exps_to_run(i)
    try
%      [~,~,time_posteriors{this_exp},new_event_detection{this_exp}] = ...
         [poisson_prob_maps_new_detection_2017_wider_window{this_exp}] = ...
...%            glm_xcorr_all_pairs('/media/shababo/data/old-data',dates(this_exp),slice_nums(this_exp),...
...%            cell_nums(this_exp),tags(this_exp),trials_detection(this_exp),glm_dummy);
             glm_from_map_est(new_event_detection_2017(this_exp),glm_dummy);
...                ROI_VB_som_data(map_estimates(this_exp),[1 pair_id(this_exp)]);
    catch e
        disp([num2str(this_exp) ' fail'])
    end
end

% delete(pool);

%%

% exps_to_run = [1:5 7:12 14:17 19:25 27:29];
% exps_to_run = [16:23];
exps_to_run = [22];
% close all

glm_to_plot = glm_out_vb_good_onetimebin3;
resp_glm = glm_out_time_bins_6;

num_experiments = length(exps_to_run);
for ii = 1:length(exps_to_run)
    
    this_exp = exps_to_run(ii);
    
%     this_exp = 7;

    ch1_glm_map = zeros(21,21);
    ch2_glm_map = zeros(21,21);
    this_result = glm_to_plot{this_exp};
%     for j = length(this_result(1).sub_vb):-1:1
% %         ch1_glm_map = ch1_glm_map + reshape(glm_to_plot{this_exp}.ch1(j).glmnet_fit.beta(2:end,...
% %             find(glm_to_plot{this_exp}.ch1(j).lambda == ...
% %             glm_to_plot{this_exp}.ch1(j).lambda_min)),21,21)';
% %         
% %         ch2_glm_map = ch2_glm_map + reshape(glm_to_plot{this_exp}.ch2(j).glmnet_fit.beta(2:end,...
% %             find(glm_to_plot{this_exp}.ch2(j).lambda == ...
% %             glm_to_plot{this_exp}.ch2(j).lambda_min)),21,21)';
%         struct_tmp = this_result(1).sub_vb(j);
%         ch1_glm_map = ch1_glm_map + (1-ch1_glm_map).*reshape(struct_tmp.alpha,21,21)';
%         struct_tmp = this_result(2).sub_vb(j);
%         ch2_glm_map = ch2_glm_map + (1-ch2_glm_map).*reshape(struct_tmp.alpha,21,21)';
%     end
    if length(this_result(1).sub_vb) == 3
        for j = 1:length(this_result(1).sub_vb)
            struct_tmp = this_result(1).sub_vb(j);
            others = setdiff(1:length(this_result(1).sub_vb),j);
            other1 = this_result(1).sub_vb(others(1));
            other2 = this_result(1).sub_vb(others(2));
            ch1_glm_map = ch1_glm_map + reshape(struct_tmp.alpha,21,21)'.*(1-reshape(other1.alpha,21,21)').*(1-reshape(other2.alpha,21,21)');
            ch1_glm_map = ch1_glm_map + (1-reshape(struct_tmp.alpha,21,21)').*reshape(other1.alpha,21,21)'.*reshape(other2.alpha,21,21)';
            other1 = this_result(2).sub_vb(others(1));
            other2 = this_result(2).sub_vb(others(2));
            struct_tmp = this_result(2).sub_vb(j);
            ch2_glm_map = ch2_glm_map + reshape(struct_tmp.alpha,21,21)'.*(1-reshape(other1.alpha,21,21)').*(1-reshape(other2.alpha,21,21)');
            ch2_glm_map = ch2_glm_map + (1-reshape(struct_tmp.alpha,21,21)').*reshape(other1.alpha,21,21)'.*reshape(other2.alpha,21,21)';
        end

        ch1_glm_map = ch1_glm_map + reshape(struct_tmp.alpha,21,21)'.*reshape(other1.alpha,21,21)'.*reshape(other2.alpha,21,21)';
        ch2_glm_map = ch2_glm_map + reshape(struct_tmp.alpha,21,21)'.*reshape(other1.alpha,21,21)'.*reshape(other2.alpha,21,21)';
    elseif length(this_result(1).sub_vb) == 1
        ch1_glm_map = reshape(this_result(1).sub_vb.alpha,21,21)';
        ch2_glm_map = reshape(this_result(2).sub_vb.alpha,21,21)';
    elseif length(this_result(1).sub_vb) == 2
        struct_tmp1 = this_result(1).sub_vb(1);
        struct_tmp2 = this_result(1).sub_vb(2);
        ch1_glm_map = reshape(struct_tmp1.alpha,21,21)'.*(1-reshape(struct_tmp2.alpha,21,21)') + ...
            (1 - reshape(struct_tmp1.alpha,21,21)').*reshape(struct_tmp2.alpha,21,21)' + ...
            reshape(struct_tmp1.alpha,21,21)'.*reshape(struct_tmp2.alpha,21,21)';
        struct_tmp1 = this_result(2).sub_vb(1);
        struct_tmp2 = this_result(2).sub_vb(2);
        ch2_glm_map = reshape(struct_tmp1.alpha,21,21)'.*(1-reshape(struct_tmp2.alpha,21,21)') + ...
            (1 - reshape(struct_tmp1.alpha,21,21)').*reshape(struct_tmp2.alpha,21,21)' + ...
            reshape(struct_tmp1.alpha,21,21)'.*reshape(struct_tmp2.alpha,21,21)';
    end
        
%         ch1_glm_map = ch1_glm_map + reshape(glm_to_plot(this_exp).ch1.glmnet_fit.beta(2:end,...
%             find(glm_to_plot(this_exp).ch1.lambda == ...
%             glm_to_plot(this_exp).ch1.lambda_min)),21,21)';
%         
%         ch2_glm_map = ch2_glm_map + reshape(glm_to_plot(this_exp).ch2.glmnet_fit.beta(2:end,...
%             find(glm_to_plot(this_exp ).ch2.lambda == ...
%             glm_to_plot(this_exp).ch2.lambda_min)),21,21)';
        
        
%         imagesc(reshape(...
%             glm_to_plot{this_exp}.ch1(j).glmnet_fit.beta(2:end,...
%             find(glm_to_plot{this_exp}.ch1(j).lambda == ...
%             glm_to_plot{this_exp}.ch1(j).lambda_1se)),21,21)'); colormap gray
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
%         num_nonzeros = unique(sum(glm_to_plot{this_exp}.ch1(j).glmnet_fit.beta > 0));
%         target_num = num_nonzeros(end-1);
%         lambda_ind = find(sum(glm_to_plot{this_exp}.ch1(j).glmnet_fit.beta > 0) == target_num,1,'last');
%         ch1_glm_map = reshape(...
%             glm_to_plot{this_exp}.ch1(j).glmnet_fit.beta(2:end,lambda_ind),21,21)';
%         ch1_glm_map(ch1_glm_map < .25) = 0;

        
        
        
%                 subplot(223)
%         imagesc(resp_glm{this_exp}.ch1.resp_map); colormap gray
%         title(['Experiment ' num2str(this_exp) ': CH1'])
%         axis off
        
%         axis image
        
%         imagesc(reshape(...
%             glm_to_plot{this_exp}.ch2(j).glmnet_fit.beta(2:end,...
%             find(glm_to_plot{this_exp}.ch2(j).lambda == ...
%             glm_to_plot{this_exp}.ch2(j).lambda_1se)),21,21)'); colormap gray
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
%         num_nonzeros = unique(sum(glm_to_plot{this_exp}.ch2(j).glmnet_fit.beta > 0));
%         target_num = num_nonzeros(end-1);
%         lambda_ind = find(sum(glm_to_plot{this_exp}.ch2(j).glmnet_fit.beta > 0) == target_num,1,'last');
%         ch2_glm_map = reshape(...
%             glm_to_plot{this_exp}.ch2(j).glmnet_fit.beta(2:end,lambda_ind),21,21)';
%         ch2_glm_map(ch2_glm_map < .25) = 0;
% 
%         figure;
%         subplot(121)
%         imagesc(ch1_glm_map >= .1); colormap gray
%         title(['Experiment ' num2str(this_exp) ': CH1'])
%         axis off
%         subplot(122)
%         imagesc(ch2_glm_map >= .1); colormap gray
%         title(['Experiment ' num2str(this_exp) ': CH2'])
%         axis off
%         
% %         
% 
%         figure
%         subplot(121)
%         imagesc(poisson_prob_maps{this_exp}.ch1.resp_map >= find(poisscdf(1:50,poisson_prob_maps{this_exp}.ch1.rate) > (1 - .0000005),1,'first'))
%         caxis([0 1])
%         subplot(122)
%         imagesc(poisson_prob_maps{this_exp}.ch2.resp_map >= find(poisscdf(1:50,poisson_prob_maps{this_exp}.ch2.rate) > (1 - .0000005),1,'first'))
%         caxis([0 1])
%         title(['Experiment: ' num2str(this_exp)])


%         figure
%         subplot(121)
%         imagesc((poisson_prob_maps{this_exp}.ch1.resp_map >= find(poisscdf(1:50,poisson_prob_maps{this_exp}.ch1.rate) > (1 - .05),1,'first')) - (ch1_glm_map >= .1))
%         caxis([0 1])
%         subplot(122)
%         imagesc((poisson_prob_maps{this_exp}.ch2.resp_map >= find(poisscdf(1:50,poisson_prob_maps{this_exp}.ch2.rate) > (1 - .05),1,'first')) - (ch2_glm_map >= .1))
%         caxis([0 1])

        input_maps_vb(this_exp).ch1 = ch1_glm_map;
        input_maps_vb(this_exp).ch2 = ch2_glm_map;

        
        
%         axis image
%     end
%     subplot(133)
%     imagesc(xcorr_images_new{this_exp})
%     title(['Experiment ' num2str(this_exp)])
%         subplot(224)
%         imagesc(resp_glm{this_exp}.ch2.resp_map); colormap gray
%         title(['Experiment ' num2str(this_exp) ': CH1'])
%         axis off
% 
    this_exp = 8; 
    data_file = ...
         [dates{this_exp} '_slice' num2str(slice_nums(this_exp)) '_cell' ...
        num2str(cell_nums(this_exp)) '' tags{this_exp} '.mat']
    load(data_file)
    [map_ch1,map_ch2] = see_grid(data,trials{this_exp},map_indices{map_index_id(this_exp)},1,1323);
    title(['Experiment ' num2str(this_exp)])
%     
%     figure
%     imagesc(exp_shuffle_stats_map_est2_intersect(this_exp).input_map)
%     title(['Experiment: ' num2str(this_exp) ', ' id_names{pair_id(this_exp) + 1}])
     
    
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
%FDR
% close all
input_maps_2017_ww_q10 = struct();
exps_to_run = 1:31;
% exps_to_run = [1 3 20 28];
for ii = 1:length(exps_to_run)
    
    this_exp = exps_to_run(ii);

    
%     poiss_cdf_map_ch1{this_exp} =  poisscdf(poisson_prob_maps{this_exp}.ch1.resp_map,mean(poisson_prob_maps{this_exp}.ch1.baseline_map(:)));
%     poiss_cdf_map_ch2{this_exp} =  poisscdf(poisson_prob_maps{this_exp}.ch2.resp_map,mean(poisson_prob_maps{this_exp}.ch2.baseline_map(:)));

    fdr_rate = 1/10;
    test_vec = (1:441)/441 * fdr_rate;
%     tect_vec = .05*ones(size(test_vec));
    
    ch1_pvals = poisson_prob_maps_new_detection_2017_wider_window{this_exp}.ch1.p_val;
    pvals_vec = ch1_pvals(:);
    [pvals_vec,inds] = sort(pvals_vec);
    
    thresh_ind = find(pvals_vec' <= test_vec,1,'last');
    detected_inds = inds(1:thresh_ind);
    [is,js] = ind2sub([21 21],detected_inds);
    input_maps_2017_ww_q10(this_exp).ch1 = zeros(21,21);
    for i = 1:length(is)
        input_maps_2017_ww_q10(this_exp).ch1(is(i),js(i)) = 1;
    end

    ch2_pvals = poisson_prob_maps_new_detection_2017_wider_window{this_exp}.ch2.p_val;
    pvals_vec = ch2_pvals(:);
    [pvals_vec,inds] = sort(pvals_vec);
    
    thresh_ind = find(pvals_vec' <= test_vec,1,'last');
    detected_inds = inds(1:thresh_ind);
    [is,js] = ind2sub([21 21],detected_inds);
    input_maps_2017_ww_q10(this_exp).ch2 = zeros(21,21);
    for i = 1:length(is)
        input_maps_2017_ww_q10(this_exp).ch2(is(i),js(i)) = 1;
    end

%     alpha_thresh = .25;    
    
%     figure
%     colormap parula
%     
%     subplot(241)
%     imagesc(poisson_prob_maps_new_detection_2017_wider_window{this_exp}.ch1.baseline_map)
% %     caxis([0 12])
%     axis off
%     colormap parula
%     title(sprintf(['Experiment ' num2str(this_exp) '\nCell 1\nBG Map w/ Est Null Rate: ' num2str(poisson_prob_maps_new_detection_2017_wider_window{this_exp}.ch1.rate)]))
%     subplot(245)
%     imagesc(poisson_prob_maps_new_detection_2017_wider_window{this_exp}.ch2.baseline_map)
% %     caxis([0 12])
%     axis off
% %     colormap hot
%     title(sprintf(['Cell 2\nBG Map w/ Est Null Rate: ' num2str(poisson_prob_maps_new_detection_2017_wider_window{this_exp}.ch2.rate)]))
%     subplot(242)
%     imagesc(poisson_prob_maps_new_detection_2017_wider_window{this_exp}.ch1.resp_map)
% %     caxis([0 12])
%     axis off
% %     colormap hot
%     title('Response Map')
%     subplot(246)
%     imagesc(poisson_prob_maps_new_detection_2017_wider_window{this_exp}.ch2.resp_map)
% %     caxis([0 12])
%     axis off
%     title('Response Map')
%     subplot(243)
% %     colormap hot
%     
%     imagesc(input_maps_2017_wider_window_q10(this_exp).ch1)
%     axis off
%     caxis([0 1])
%     title(sprintf(['FDR Detection\nFDR Rate (q) = ' num2str(fdr_rate)]))
% %     colormap gray
%     subplot(247)
%     imagesc(input_maps_2017_wider_window_q10(this_exp).ch2)
%     axis off
% %     colormap gray
%     caxis([0 1])
%     title(sprintf(['FDR Detection\nFDR Rate (q) = ' num2str(fdr_rate)]))
%     subplot(244)
% %     imagesc(input_maps_vb(this_exp).ch1 > alpha_thresh)
%     imagesc(input_maps_2017(this_exp).ch1)
%     axis off
% %     colormap gray
%     caxis([0 1])
%     title(sprintf(['Previous Detection']))
%     subplot(248)
% %     imagesc(input_maps_vb(this_exp).ch2 > alpha_thresh)
%     imagesc(input_maps_2017(this_exp).ch2)
%     axis off
% %     colormap gray
%     caxis([0 1])
%     title(sprintf(['Previous Detection']))
    
%     set(gcf,'position',maxwinsize);
    
%     export_fig supplemntal-rate-maps.pdf -painters
%     close all
%     data_file = ...
%         [dates{this_exp} '_slice' num2str(slice_nums(this_exp)) '_cell' ...
%         num2str(cell_nums(this_exp)) '' tags{this_exp} '.mat']
%     load(data_file)
%     [map_ch1,map_ch2] = see_grid(data,trials{this_exp},map_indices{map_index_id(this_exp)},1,1323);
%     title(['Experiment ' num2str(this_exp)])
end
%%
exp = 21;
figure
l4_resp_map = poisson_prob_maps_new_detection_2017{exp}.ch2.resp_map;
l4_input_map = input_maps_2017_ww_q10(exp).ch2;
l5_resp_map = poisson_prob_maps_new_detection_2017{exp}.ch1.resp_map;
l5_input_map = input_maps_2017_ww_q10(exp).ch1;
l5_resp_map(l5_input_map == 0) = 0; %only sites that exceed 'background rate'
l4_resp_map(l4_input_map == 0) = 0;
[Acol Arow]=ind2sub([21 21],find(l4_resp_map));%);
[Bcol Brow]=ind2sub([21 21],find(l5_resp_map));%;
% subplot(121)
% hold on
numbers = repmat(1:21,21,1);
numbers_t = numbers';

scatter(numbers(:),numbers_t(:),10,[0 0 0],'.')
hold on
hs = scatter(Arow,21 - Acol + 1,200*l4_resp_map(find(l4_resp_map)),[247 158 29]/255,'.','LineWidth',2);%find(l4_resp_map)
% alpha(hs,.5)
xlim([1 21])
ylim([1 21])
axis off
% subplot(122)
scatter(numbers(:),numbers_t(:),10,[0 0 0],'.')
hold on
hs = scatter(Brow,21 - Bcol + 1,l5_resp_map(find(l5_resp_map))*20,[39 170 225]/255,'o','LineWidth',2);
% alpha(hs,.5)
xlim([1 21])
ylim([1 21])
axis off
axis image


%%
exps_to_run = 2
;setdiff(1:31,exclude);
for i = 1:length(exps_to_run)
    exp = exps_to_run(i);
%     exp = 21;
    figure
    l4_resp_map = poisson_prob_maps_new_detection_2017{exp}.ch2.resp_map;
    l4_input_map = input_maps_2017_ww_q10(exp).ch2;
    l5_resp_map = poisson_prob_maps_new_detection_2017{exp}.ch1.resp_map;
    l5_input_map = input_maps_2017_ww_q10(exp).ch1;
    l5_resp_map(l5_input_map == 0) = 0; %only sites that exceed 'background rate'
    l4_resp_map(l4_input_map == 0) = 0;
    [Acol Arow]=ind2sub([21 21],find(l4_resp_map));%);
    [Bcol Brow]=ind2sub([21 21],find(l5_resp_map));%;
    subplot(121)
    % hold on
    numbers = repmat(1:21,21,1);
    numbers_t = numbers';

    scatter(numbers(:),numbers_t(:),10,[0 0 0],'.')
    hold on
    if any(exp == l4_pairs)
        color = [247 158 29]/255;
    else
        color = [39 170 225]/255;
    end
    hs = scatter(Arow,21 - Acol + 1,200*l4_input_map(find(l4_resp_map))*6*4,color,'.','LineWidth',2);%find(l4_resp_map)
    % alpha(hs,.5)
    xlim([1 21])
    ylim([1 21])
    axis off
    axis image
    subplot(122)
    scatter(numbers(:),numbers_t(:),10,[0 0 0],'.')
    hold on
    
    hs = scatter(Brow,21 - Bcol + 1,l5_input_map(find(l5_resp_map))*6*4*200,[39 170 225]/255,'.','LineWidth',2);
    % alpha(hs,.5)
    xlim([1 21])
    ylim([1 21])
    axis off
    axis image
    set(gcf,'position',maxwinsize);
    export_fig(['input_map_scatter_2017_ww_q10_exp' num2str(exp) '.eps'])
%     saveas(gcf,['input_map_scatter_2017_ww_q10_exp' num2str(exp) '.eps'])
%     close(gcf)
end

%%
clc
exps_to_run = l5_pairs;

for i = 1:length(exps_to_run)
    this_exp = exps_to_run(i)
    [dates{this_exp} '_slice' num2str(slice_nums(this_exp)) '_cell' ...
        num2str(cell_nums(this_exp)) '' tags{this_exp} '.mat']
end
%%
exps_to_run = setdiff(1:31,exclude);
% exps_to_run = [20];
for ii = 1:length(exps_to_run)
    
    this_exp = exps_to_run(ii);

    
%     poiss_cdf_map_ch1{this_exp} =  poisscdf(poisson_prob_maps{this_exp}.ch1.resp_map,mean(poisson_prob_maps{this_exp}.ch1.baseline_map(:)));
%     poiss_cdf_map_ch2{this_exp} =  poisscdf(poisson_prob_maps{this_exp}.ch2.resp_map,mean(poisson_prob_maps{this_exp}.ch2.baseline_map(:)));

%     fdr_rate = .1;
%     test_vec = (1:441)/441 * fdr_rate;
%     
%     ch1_pvals = poisson_prob_maps_fix_check{this_exp}.ch1.p_val;
%     pvals_vec = ch1_pvals(:);
%     [pvals_vec,inds] = sort(pvals_vec);
%     
%     thresh_ind = find(pvals_vec' <= test_vec,1,'last');
%     detected_inds = inds(1:thresh_ind);
%     [is,js] = ind2sub([21 21],detected_inds);
%     input_maps_fix(this_exp).ch1 = zeros(21,21);
%     for i = 1:length(is)
%         input_maps_fix(this_exp).ch1(is(i),js(i)) = 1;
%     end
% 
%     ch2_pvals = poisson_prob_maps_fix_check{this_exp}.ch2.p_val;
%     pvals_vec = ch2_pvals(:);
%     [pvals_vec,inds] = sort(pvals_vec);
%     
%     thresh_ind = find(pvals_vec' <= test_vec,1,'last');
%     detected_inds = inds(1:thresh_ind);
%     [is,js] = ind2sub([21 21],detected_inds);
%     input_maps_fix(this_exp).ch2 = zeros(21,21);
%     for i = 1:length(is)
%         input_maps_fix(this_exp).ch2(is(i),js(i)) = 1;
%     end

%     alpha_thresh = .25;    
%     
    figure
    subplot(121)
    imagesc(input_maps_2017_ww_q10(this_exp).ch1)
    caxis([0 1])
%     caxis([0 12])
    axis off
    axis image
%     title('Input Map Ch 1')
    subplot(122)
    imagesc(input_maps_2017_ww_q10(this_exp).ch2)
    caxis([0 1])
%     caxis([0 12])
    colormap gray
    axis off
    axis image
%     title('Input Map Ch 2')
    set(gcf,'position',maxwinsize);
%     colormap grey
%     export_fig supplemntal-rate-maps.pdf -painters
    saveas(gcf,['input_map_207_ww_q10_exp' num2str(this_exp) '.jpg'])
    
%     close all
% %     
%     data_file = ...
%         [dates{this_exp} '_slice' num2str(slice_nums(this_exp)) '_cell' ...
%         num2str(cell_nums(this_exp)) '' tags{this_exp} '.mat']
%     load(data_file)
%     [map_ch1,map_ch2] = see_grid(data,trials{this_exp},map_indices{map_index_id(this_exp)},1,1323);
%     title(['Experiment ' num2str(this_exp)])
end

%% load data file
this_exp = 6;
    data_file = ...
        [dates{this_exp} '_slice' num2str(slice_nums(this_exp)) '_cell' ...
        num2str(cell_nums(this_exp)) '' tags{this_exp} '.mat']
    [map_ch1,~] = see_grid(data,trials{this_exp},map_indices{map_index_id(this_exp)},1,1323);
figure
    plot_trace_stack_grid(map_ch1,6,5,0);
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
% 
% [glm_out_doubleA,xcorr_images_doubleA,samples_psths_doubleA] = ...
%         glm_xcorr_all_pairs('/media/shababo/data/old-data',dates(this_exp),slice_nums(this_exp),...
%         cell_nums(this_exp),tags(this_exp),trials_detection(this_exp));

%%

exps_to_run = [1:5 7:8 11:12 14:17 19:25 27:30];
% exps_to_run = 24;
exps_to_run = 1:30;

for ii = 1:length(exps_to_run)

    this_exp = exps_to_run(ii)
    
    try

        events1 = samples_psths_new{this_exp}{1};
        events2 = samples_psths_new{this_exp}{2};

        disp('got events')
        input_locs1 = find(glm_out_newer(this_exp).ch1.glmnet_fit.beta(2:end,find(glm_out_newer(this_exp).ch1.lambda == ...
                glm_out_newer(this_exp).ch1.lambda_1se)+2));%+5*(-pair_id(this_exp) + 1)
        input_locs2 = find(glm_out_newer(this_exp).ch2.glmnet_fit.beta(2:end,find(glm_out_newer(this_exp).ch2.lambda == ...
            glm_out_newer(this_exp).ch2.lambda_1se)+2+5*(-pair_id(this_exp)+1)));

    %     all_locs = union(input_locs1,input_locs2);
    %     all_locs = intersect(input_locs1,input_locs2);
%         all_locs = 1:441;
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
    catch e
        disp([num2str(this_exp) ' fail'])
    end
%     disp('result:')
%     this_exp
%     pair_id(this_exp)
%     common_input_score_norm_intersect(this_exp)
    
    
end

%%

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
    input_locs1 = find(input_maps_fix2(this_exp).ch1(:) > .1);
    input_locs2 = find(input_maps_fix2(this_exp).ch2(:) > .1);
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

%%
% close all
for i = 20
    this_exp = i;
    try
        figure
        histogram(exp_shuffle_stats_map_est4_intersect(this_exp).locs(1).null_dist,1000)
        title(['Experiment: ' num2str(this_exp) ', ' id_names{pair_id(this_exp) + 1}])
    catch e
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
scatter(cis_tmp(pair_id_tmp),max(counts)/2*ones(length(l4_pairs),1))
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
           zeros(1,22) ones(1,5) zeros(1,54)
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
% ylabel('Synchrony Probability')close a
% legend({'L4-L5 Pairs','L5-L5 Pairs'})

% subplot(2,1,2)
scatter(ones(sum(~pair_id_tmp),1),cis_tmp(~pair_id_tmp),'k' ,'jitter','on', 'jitterAmount',0.05);
hold on
scatter(1.25,mean(cis_tmp(~pair_id_tmp)),100,[0 0 0],'filled')
hold on
plot(1.25*[1 1],mean(cis_tmp(~pair_id_tmp)) + std(cis_tmp(~pair_id_tmp))/sqrt(length(cis_tmp(~pair_id_tmp)))*[-1 1],'k','Linewidth',2)
hold on
scatter(2*ones(length(l4_pairs),1),cis_tmp(pair_id_tmp),'k','jitter','on', 'jitterAmount',0.05);
hold on
scatter(2.25,mean(cis_tmp(pair_id_tmp)),100,[0 0 0],'filled')
hold on
plot(2.25*[1 1],mean(cis_tmp(pair_id_tmp)) + std(cis_tmp(pair_id_tmp))/sqrt(length(cis_tmp(pair_id_tmp)))*[-1 1],'k','Linewidth',2)
hold on
xlim([0 3])
set(gca,'Xtick',[1 2])
set(gca,'xticklabels',{['L4-L5 Pairs (N = ' num2str(sum(~pair_id_tmp)) ')'],['L5-L5 Pairs (N = ' num2str(length(l4_pairs)) ')']})
ylabel('Common Input Locations')

[p,h] = ranksum(cis_tmp(~pair_id_tmp),cis_tmp(pair_id_tmp),'tail','left');
title(['p = ' num2str(p)])

%%
cis_tmp_null = num_shared_locs_null./num_event_locs_null;
cis_tmp_null(skip_exps) = [];

cis_tmp = num_shared_locs./num_event_locs;
cis_tmp(skip_exps) = [];

figure

scatter(2.5*ones(length(l4_pairs),1),cis_tmp(pair_id_tmp),'k');
hold on
scatter(2.5,mean(cis_tmp(pair_id_tmp)),100,[0 0 0],'filled')
hold on
plot(2.5*[1 1],mean(cis_tmp(pair_id_tmp)) + std(cis_tmp(pair_id_tmp))/sqrt(length(cis_tmp(pair_id_tmp)))*[-1 1],'k','Linewidth',2)
hold on

scatter(3*ones(length(l4_pairs),1),cis_tmp_null(pair_id_tmp),'k');
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
set(gca,'xticklabels',{['L4-L5 Pairs (N = ' num2str(sum(~pair_id_tmp)) ')'],'L4-L5 Shuffle',['L5-L5 Pairs (N = ' num2str(length(l4_pairs)) ')'],'L5-L5 Shuffle'})
ylabel('Common Input Probability')

[pl4,h] = ranksum(cis_tmp(~pair_id_tmp),cis_tmp_null(~pair_id_tmp),'tail','right')
[pl5,h] = ranksum(cis_tmp(pair_id_tmp),cis_tmp_null(pair_id_tmp),'tail','right')
title(['p_l4 = ' num2str(pl4) ', p_l5 = ' num2str(pl5)])

% xlim([0 3])
% set(gca,'Xtick',[1 2])
% set(gca,'xticklabels',{['L4-L5 Pairs (N = ' num2str(sum(~pair_id_tmp)) ')'],['L5-L5 Pairs (N = ' num2str(length(l4_pairs)) ')']})
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

%NOTE:

% 
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
skip_exps = [1 10 15];
% pair_id_tmp = pair_id;
% pair_id_tmp(skip_exps) = [];
% 
exps_to_run = 1:length(dates); %  13 26
exps_to_run = setdiff(exps_to_run,skip_exps)
%%
% close all
%issues 13,19,20 7 14 15
% 16 too low
% exps_to_run = [1:length(dates)];

% detection_grids = cell(1,length(dates));

thresholds = [40 50 0 450 0 10 7.5 15 15 30 0 30 5 6 6 8 10 ...
    10 15 15 15 15 20 20 15 15 15 15];
exps_to_run = [19];
figure
for ii = 1:length(exps_to_run)
    
    this_exp = exps_to_run(ii);
    data_file = ...
        [dates{this_exp} '_slice' num2str(slice_nums(this_exp)) '_cell' ...
        num2str(cell_nums(this_exp)) '' tags{this_exp} '.mat'];
    load(data_file)
%    figure
     subplot(1,1,ii)
    [map_ch1,~] = see_grid(data,trials{this_exp},map_index_spikes,0,363);
    plot_trace_stack_grid(map_ch1,3,10,0);
 %   title(['Experiment ' num2str(this_exp)])
%     
%      this_row_i = row_inds{this_exp};
%     this_col_j = col_inds{this_exp};
%     detection_results = cell(length(this_row_i),length(this_col_j));
%     
%     for i = 1:length(this_row_i)
%         for j = 1:length(this_col_j)
%             switch patch_type(this_exp)
%                 case 1
%                     detection_results{i,j} = detect_peaks(...
% ...%                         -1.0*bsxfun(@minus,map_ch1{i,j},median(map_ch1{i,j}(:,1:100),2)),5,20,1,1,0,1)*70;
%                             -1.0*bsxfun(@minus,map_ch1{this_row_i(i),this_col_j(j)},median(map_ch1{this_row_i(i),this_col_j(j)},2)),thresholds(this_exp),80,1,1,0,0,0)*70;
%                 case 2
%                     detection_results{i,j} = detect_peaks(...
%                         1.0*bsxfun(@minus,map_ch1{this_row_i(i),this_col_j(j)},map_ch1{this_row_i(i),this_col_j(j)}(:,1)),thresholds(this_exp),80,1,1,0,0,0)*70;
%                 case 3
%                     detection_results{i,j} = mean(map_ch1{this_row_i(i),this_col_j(j)});
%             end
%         end
%     end
%     
%     detection_grids{this_exp} = detection_results;
%     subplot(122)
%     plot_trace_stack_grid(detection_results,Inf,1,0);
end

%%
exps_to_run = setdiff(1:length(dates),union([6],find(patch_type == 3)));

centers = zeros(length(dates),2);
centers(17,:) = [4 3];
centers(13:14,:) = [5 6; 5 6];
spike_counts = zeros(11,11,length(dates));
first_spike_time = nan(11,11,length(dates));
first_spike_time_trials = cell(11,11,length(dates));
jitters = nan(11,11,length(dates));
for ii = 1:length(exps_to_run)
%     ii = 1
    this_exp = exps_to_run(ii)
    
    if isequal(centers(this_exp,:),[0 0])
        center = ceil(size(detection_grids{this_exp})/2);
    else
        center = centers(this_exp,:);
    end
    spike_counts_trials(this_exp) = size(detection_grids{this_exp}{1,1},1);
    for i = 1:size(detection_grids{this_exp},1)
        for j = 1:size(detection_grids{this_exp},2)
            
            main_i = i - center(1) + ceil(size(spike_counts,1)/2);
            main_j = j - center(2) + ceil(size(spike_counts,2)/2);
            if main_i > 11 || main_j > 11
                break
            end
            
            spike_counts(main_i,main_j,this_exp) = ...
                sum(sum(detection_grids{this_exp}{i,j}(:,100:end)/70,2));
            spike_counts_trials(this_exp) = size(detection_grids{this_exp}{i,j},1);
%             first_spike_time(main_i,main_j,this_exp) = ...
%                 arrayfun(@(x) find(x,1,'first'),detection_grids{this_exp}{i,j}(:,100:end));
            if spike_counts(main_i,main_j,this_exp) > 0
                count = 1;
                for k = 1:size(detection_grids{this_exp}{i,j},1)
                    time = find(detection_grids{this_exp}{i,j}(k,100:end),1,'first');
                    if ~isempty(time)
                        if count == 1
                            first_spike_time(main_i,main_j,this_exp) = 0;
                            first_spike_time_trials{main_i,main_j,this_exp} = time;
                        else
                            first_spike_time_trials{main_i,main_j,this_exp} = [first_spike_time_trials{main_i,main_j,this_exp} time];
                        end
                        first_spike_time(main_i,main_j,this_exp) = ...
                            first_spike_time(main_i,main_j,this_exp) + time;
                        count = count + 1;
                    end
                end
                first_spike_time(main_i,main_j,this_exp) = ...
                    first_spike_time(main_i,main_j,this_exp)/count;
                if length(first_spike_time_trials{main_i,main_j,this_exp}) > 4
%                     jitters(main_i,main_j,this_exp) = 2*sqrt(2*log(2))*std(first_spike_time_trials{main_i,main_j,this_exp});
                      jitters(main_i,main_j,this_exp) = max(first_spike_time_trials{main_i,main_j,this_exp}) - min(first_spike_time_trials{main_i,main_j,this_exp}) ;
                      
%                     if jitters(main_i,main_j,this_exp) < 10
%                         return
%                     end
                else
                    first_spike_time(main_i,main_j,this_exp) = NaN;
                end
            end
            
        end
    end
    

end

%%

for ii = 1:length(high_power_inds) 
    this_exp = high_power_inds(ii);
    prob_spike(:,:,this_exp) = cellfun(@(x) length(x),squeeze(first_spike_time_trials(:,:,this_exp)))/spike_counts_trials(this_exp);
end

%%

for ii = 1:length(exps_to_run) 
    this_exp = exps_to_run(ii);
    spike_count_means(:,:,this_exp) = spike_counts(:,:,this_exp)/spike_counts_trials(this_exp);
end

%%

figure
x = (1:11)';
mean_vec = squeeze(mean(spike_count_means(:,6,high_pow_exps),3));
error = squeeze(std(spike_count_means(:,6,high_pow_exps),[],3))/sqrt(length(high_pow_exps));  
fill([x;flipud(x)],[(mean_vec-error);flipud(mean_vec+error)],[1 1 1]*.75,'linestyle','none');
hold on
plot(x,mean_vec,'Color',[0 0 0]/255)

figure
x = (1:11)';
mean_vec = squeeze(mean(spike_count_means(6,:,high_power_inds),3))';
error = squeeze(std(spike_count_means(6,:,high_power_inds),[],3))'/sqrt(length(high_power_inds));  
fill([x;flipud(x)],[(mean_vec-error);flipud(mean_vec+error)],[1 1 1]*.75,'linestyle','none');
hold on
plot(x,mean_vec,'Color',[0 0 0]/255)

%%

figure
x = (1:11)';
mean_vec = squeeze(mean(prob_spike(:,6,high_power_inds),3));
error = squeeze(std(prob_spike(:,6,high_power_inds),[],3))/sqrt(length(high_power_inds));  
fill([x;flipud(x)],[(mean_vec-error);flipud(mean_vec+error)],[1 1 1]*.75,'linestyle','none');
hold on
plot(x,mean_vec,'Color',[0 0 0]/255)

figure
x = (1:11)';
mean_vec = squeeze(mean(prob_spike(6,:,high_power_inds),3))';
error = squeeze(std(prob_spike(6,:,high_power_inds),[],3))'/sqrt(length(high_power_inds));  
fill([x;flipud(x)],[(mean_vec-error);flipud(mean_vec+error)],[1 1 1]*.75,'linestyle','none');
hold on
plot(x,mean_vec,'Color',[0 0 0]/255)

%%
for ii = 1:length(high_power_inds)
    this_exp = high_power_inds(ii);
    num_spike_locations(this_exp) = sum(sum(squeeze(prob_spike(:,:,this_exp) >= .75)));
end

figure
histogram(num_spike_locations(high_power_inds))

%%

figure
% subplot(121)
scatter(first_spike_time(6,6,high_power_inds),spike_count_means(6,6,high_power_inds))
subplot(122)
scatter(jitters(6,6,high_power_inds),spike_count_means(6,6,high_power_inds))


%%
figure;

jitters_red = squeeze(jitters([5 7],[5 7],:));
jitters_red = jitters_red(~isnan(jitters_red));
ecdf(jitters_red(:))



%%

for ii = 1:length(exps_to_run)
%     ii = 1
    this_exp = exps_to_run(ii)
    
    num_spike_locations(this_exp) = sum(sum(spike_counts(:,:,this_exp) >= .75));
    
end
    
%%
exclude = [4 2 9 10 24 27 7 15 17 1 31]; % added 6 for now
l4_pairs = setdiff(find(pair_id == 0),exclude);%[1 4 5 9 17 18 19 21 22 23 30];
l5_pairs = setdiff(find(pair_id == 1),exclude);%[11 12 14 15 20 26 28];
exps = union(l4_pairs,l5_pairs);
% num_locs_exps = sum(num_locs(exps));
% 
% p_val_thresh = .05/num_locs_exps;

fdr_rate = .05;
all_p_vals_l4 = [];
all_p_vals_l5 = [];
% exps = 13
for i = exps
    
    
    
    if num_locs_union_full_ww_q10(i) > 0 && num_locs_full_ww_q10(i) > 0
        num_doubleAs = num_locs_full_ww_q10(i);
        [pvals_vec,inds] = sort([exp_shuffle_stats_event_detection_2017_ww_q10(i).locs(1:num_locs_full_ww_q10(i)).p_val]);
        test_vec = (1:num_doubleAs)/num_doubleAs * fdr_rate;
%         test_vec = .05*ones(size(test_vec));
        thresh_ind = find(pvals_vec <= test_vec,1,'last');
%         detected_inds = inds(1:thresh_ind);
%         for j = 1:length(detected_inds)
%             exp_shuffle_stats_event_detection_2017_ww_q10(i).locs(detected_inds(j)).detection = 1;
%         end
%         for j = setdiff(1:length(exp_shuffle_stats_event_detection_2017_ww_q10(i).locs),detected_inds)
%             exp_shuffle_stats_event_detection_2017_ww_q10(i).locs(j).detection = 0;
%         end
%         if thresh_ind > 0
%             ci_score(i) = thresh_ind/num_locs_union_full_ww_q10(i);
%             
%         else
%             ci_score(i) = 0;
%         end
        ci_score2(i) = num_locs_full_ww_q10(i)/num_locs_union_full_ww_q10(i);
%         if any(l4_pairs == i)
%             all_p_vals_l4 = [all_p_vals_l4 pvals_vec];
%         else
%             all_p_vals_l5 = [all_p_vals_l5 pvals_vec];
%         end
%         ci_score(i) = sum(sum(exp_shuffle_stats_map_est6_intersect(i).input_map < p_val_thresh))/num_locs_union(i);
    else
%         ci_score(i) = 0;
        ci_score2(i) = 0;
    end
    
end

cis_tmp = ci_score;
% cis_tmp(29:30) = 0;
% pair_id_tmp = zeros(1,30);
% pair_id_tmp(exps) = 1;
% pair_id_tmp = logical(pair_id_tmp);
figure

scatter(2.75*ones(length(cis_tmp(l4_pairs)),1),cis_tmp(l4_pairs),'k','jitter','on', 'jitterAmount',0.05)%,'k' ,'jitter','on', 'jitterAmount',0.25);
hold on
% scatter(1.25*ones(length(ci_score2(l4_pairs)),1),ci_score2(l4_pairs),'k' ,'jitter','on', 'jitterAmount',0.1);
% hold on
% plot([1.75;1.25],[cis_tmp(l4_pairs);ci_score2(l4_pairs)],'k')
hold on
scatter(2.5,mean(cis_tmp(l4_pairs)),100,[0 0 0],'filled')
hold on
plot(2.5*[1 1],mean(cis_tmp(l4_pairs)) + std(cis_tmp(l4_pairs))/sqrt(length(cis_tmp(l4_pairs)))*[-1 1],'k','Linewidth',2)
hold on
% scatter(1.0,mean(ci_score2(l4_pairs)),100,[0 0 0],'filled')
% hold on
% plot(1.0*[1 1],mean(ci_score2(l4_pairs)) + std(ci_score2(l4_pairs))/sqrt(length(ci_score2(l4_pairs)))*[-1 1],'k','Linewidth',2)
% hold on
% scatter(3*ones(length(l4_pairs),1),cis_tmp_null(pair_id_tmp),'k');
% hold on
% scatter(3,mean(cis_tmp_null(pair_id_tmp)),100,[0 0 0],'filled')
% hold on
% plot([3 3],mean(cis_tmp_null(pair_id_tmp)) + std(cis_tmp_null(pair_id_tmp))/sqrt(length(cis_tmp_null(pair_id_tmp)))*[-1 1],'k','Linewidth',2)
% hold on



scatter(3.5*ones(length(cis_tmp(l5_pairs)),1),cis_tmp(l5_pairs),'k','jitter','on', 'jitterAmount',0.05)% ,'jitter','on', 'jitterAmount',0.25);
hold on
% scatter(2.0*ones(length(ci_score2(l5_pairs)),1),ci_score2(l5_pairs),'k' ,'jitter','on', 'jitterAmount',0.05);
% hold on
% plot([3.25;2.75],[cis_tmp(l5_pairs);ci_score2(l5_pairs)],'k')
hold on
scatter(3.25,mean(cis_tmp(l5_pairs)),100,[0 0 0],'filled')
hold on
plot(3.25*[1 1],mean(cis_tmp(l5_pairs)) + std(cis_tmp(l5_pairs))/sqrt(length(cis_tmp(l5_pairs)))*[-1 1],'k','Linewidth',2)
% hold on
% scatter(1.75,mean(ci_score2(l5_pairs)),100,[0 0 0],'filled')
% hold on
% plot(1.75*[1 1],mean(ci_score2(l5_pairs)) + std(ci_score2(l5_pairs))/sqrt(length(ci_score2(l5_pairs)))*[-1 1],'k','Linewidth',2)

%  set(gca,'yscale','log');
% scatter(1.5*ones(sum(~pair_id_tmp),1),cis_tmp_null(~pair_id_tmp),'k');
% hold on
% scatter(1.5,mean(cis_tmp_null(~pair_id_tmp)),100,[0 0 0],'filled')
% hold on
% plot([1.5 1.5],mean(cis_tmp_null(~pair_id_tmp)) + std(cis_tmp_null(~pair_id_tmp))/sqrt(length(cis_tmp_null(~pair_id_tmp)))*[-1 1],'k','Linewidth',2)
% hold on

xlim([2.25 4])
set(gca,'Xtick',[])
% set(gca,'xticklabels',{['L4-L5 Pairs (N = ' num2str(length(l4_pairs)) ')'],['L5-L5 Pairs (N = ' num2str(length(l5_pairs)) ')']})
ylabel('Coincidence Probability')

[p,h] = ranksum(cis_tmp(l4_pairs),cis_tmp(l5_pairs),'tail','left');
title(['p = ' num2str(p)])

% figure;
% edges = logspace(-5,0,20);
% histogram(all_p_vals_l5,edges,'Normalization','probability');
% hold on
% histogram(all_p_vals_l4,edges,'Normalization','probability');
% set(gca,'xscale','log')

%% do FDR across all p-vals

l4_pairs = setdiff(find(pair_id == 0),exclude);%[1 4 5 9 17 18 19 21 22 23 30];
l5_pairs = setdiff(find(pair_id == 1),exclude);%[11 12 14 15 20 26 28];
exps = union(l4_pairs,l5_pairs);
% num_locs_exps = sum(num_locs_ww_q10(exps));

% p_val_thresh = .05/num_locs_exps;

fdr_rate = .05;
all_p_vals_l4 = [];
all_p_vals_l5 = [];
% exps = 13
all_pvals = [];
for i = exps
    
    
    
    if num_locs_union_full_ww_q10(i) > 0 && num_locs_full_ww_q10(i) > 0
        num_doubleAs = num_locs_full_ww_q10(i);
        pvals_vec = [exp_shuffle_stats_event_detection_2017_ww_q10(i).locs(1:num_locs_full_ww_q10(i)).p_val];
        pval_inds{i} = (1:length(pvals_vec)) + length(all_pvals);
        
        all_pvals = [all_pvals pvals_vec];
        
        if any(l4_pairs == i)
            all_p_vals_l4 = [all_p_vals_l4 pvals_vec];
        else
            all_p_vals_l5 = [all_p_vals_l5 pvals_vec];
        end
        
    end
end
[all_pvals,inds] = sort(all_pvals);

test_vec = (1:length(all_pvals))/length(all_pvals) * fdr_rate;


thresh_ind = find(all_pvals <= test_vec,1,'last');
detected_inds = inds(1:thresh_ind);


for i = exps
    
    if num_locs_union_full_ww_q10(i) > 0 && num_locs_full_ww_q10(i) > 0
        detect_count(i) = 0;
        for j = 1:length(pval_inds{i})
            this_ind = pval_inds{i}(j) - min(pval_inds{i}) + 1;
            if any(pval_inds{i}(j) == detected_inds)
                exp_shuffle_stats_event_detection_2017_ww_q10(i).locs(this_ind).detection = 1;
                detect_count(i) = detect_count(i) + 1;
            else
                exp_shuffle_stats_event_detection_2017_ww_q10(i).locs(this_ind).detection = 0;
                
            end
        end
        


        if thresh_ind > 0
            ci_score(i) = detect_count(i)/num_locs_union_full_ww_q10(i);
            
        else
            ci_score(i) = 0;
        end
       
    else
        ci_score(i) = 0;

    end
    
    
    
end

cis_tmp = ci_score;
% cis_tmp(29:30) = 0;
% pair_id_tmp = zeros(1,30);
% pair_id_tmp(exps) = 1;
% pair_id_tmp = logical(pair_id_tmp);
figure

scatter(2.75*ones(length(cis_tmp(l4_pairs)),1),cis_tmp(l4_pairs),100,'k','jitter','on', 'jitterAmount',0.2)%,'k' ,'jitter','on', 'jitterAmount',0.25);
hold on
% scatter(1.25*ones(length(ci_score2(l4_pairs)),1),ci_score2(l4_pairs),100,'k' ,'jitter','on', 'jitterAmount',2.0);
% hold on
% plot([1.75;1.25],[cis_tmp(l4_pairs);ci_score2(l4_pairs)],'k')
hold on
scatter(2.5,mean(cis_tmp(l4_pairs)),100,[0 0 0],'filled')
hold on
plot(2.5*[1 1],mean(cis_tmp(l4_pairs)) + std(cis_tmp(l4_pairs))/sqrt(length(cis_tmp(l4_pairs)))*[-1 1],'k','Linewidth',2)
hold on
% scatter(1.0,mean(ci_score2(l4_pairs)),100,[0 0 0],'filled')
% hold on
% plot(1.0*[1 1],mean(ci_score2(l4_pairs)) + std(ci_score2(l4_pairs))/sqrt(length(ci_score2(l4_pairs)))*[-1 1],'k','Linewidth',2)
% hold on
% scatter(3*ones(length(l4_pairs),1),cis_tmp_null(pair_id_tmp),'k');
% hold on
% scatter(3,mean(cis_tmp_null(pair_id_tmp)),100,[0 0 0],'filled')
% hold on
% plot([3 3],mean(cis_tmp_null(pair_id_tmp)) + std(cis_tmp_null(pair_id_tmp))/sqrt(length(cis_tmp_null(pair_id_tmp)))*[-1 1],'k','Linewidth',2)
% hold on



scatter(3.5*ones(length(cis_tmp(l5_pairs)),1),cis_tmp(l5_pairs),100,'k','jitter','on', 'jitterAmount',0.05)% ,'jitter','on', 'jitterAmount',0.25);
hold on
% scatter(2.0*ones(length(ci_score2(l5_pairs)),1),ci_score2(l5_pairs),'k' ,'jitter','on', 'jitterAmount',0.05);
% hold on
% plot([3.25;2.75],[cis_tmp(l5_pairs);ci_score2(l5_pairs)],'k')
hold on
scatter(3.25,mean(cis_tmp(l5_pairs)),100,[0 0 0],'filled')
hold on
plot(3.25*[1 1],mean(cis_tmp(l5_pairs)) + std(cis_tmp(l5_pairs))/sqrt(length(cis_tmp(l5_pairs)))*[-1 1],'k','Linewidth',2)
hold on
% scatter(1.75,mean(ci_score2(l5_pairs)),100,[0 0 0],'filled')
% hold on
% plot(1.75*[1 1],mean(ci_score2(l5_pairs)) + std(ci_score2(l5_pairs))/sqrt(length(ci_score2(l5_pairs)))*[-1 1],'k','Linewidth',2)

%  set(gca,'yscale','log');
% scatter(1.5*ones(sum(~pair_id_tmp),1),cis_tmp_null(~pair_id_tmp),'k');
% hold on
% scatter(1.5,mean(cis_tmp_null(~pair_id_tmp)),100,[0 0 0],'filled')
% hold on
% plot([1.5 1.5],mean(cis_tmp_null(~pair_id_tmp)) + std(cis_tmp_null(~pair_id_tmp))/sqrt(length(cis_tmp_null(~pair_id_tmp)))*[-1 1],'k','Linewidth',2)
% hold on

xlim([2.25 4])
set(gca,'Xtick',[])
% set(gca,'xticklabels',{['L4-L5 Pairs (N = ' num2str(length(l4_pairs)) ')'],['L5-L5 Pairs (N = ' num2str(length(l5_pairs)) ')']})
ylabel('Coincidence Probability')

[p,h] = ranksum(cis_tmp(l4_pairs),cis_tmp(l5_pairs),'tail','left');
title(['p = ' num2str(p)])

figure;
edges = logspace(-4,0,20);
histogram(all_p_vals_l5,edges,'Normalization','probability');
hold on
histogram(all_p_vals_l4,edges,'Normalization','probability');
set(gca,'xscale','log')

%%

figure;
histogram([num_locs1_full_ww_q10(l4_pairs) num_locs1_full_ww_q10(l5_pairs) num_locs2_full_ww_q10(l5_pairs)],0:5:120,'Normalization','probability'); hold on
histogram(num_locs2_full_ww_q10(l4_pairs),0:5:80,'Normalization','probability'); hold off
% axis off

figure

scatter(ones(length(l4_pairs),1),num_locs2_full_ww_q10(l4_pairs),'k')
hold on
scatter(1+ones(length([num_locs1_full_ww_q10(l4_pairs) num_locs1_full_ww_q10(l5_pairs) num_locs2_full_ww_q10(l5_pairs)]),1),...
    [num_locs1_full_ww_q10(l4_pairs) num_locs1_full_ww_q10(l5_pairs) num_locs2_full_ww_q10(l5_pairs)],'k')
xlim([0.5 2.5])

sum([num_locs1_full(l4_pairs) num_locs1_full(l5_pairs) num_locs2_full(l5_pairs)])/length([num_locs1_full(l4_pairs) num_locs1_full(l5_pairs) num_locs2_full(l5_pairs)])

sum(num_locs2_full(l4_pairs))/length(num_locs2_full(l4_pairs))

mean([num_locs1_full(l4_pairs) num_locs1_full(l5_pairs) num_locs2_full(l5_pairs)])
mean(num_locs2_full(l4_pairs))
std([num_locs1_full(l4_pairs) num_locs1_full(l5_pairs) num_locs2_full(l5_pairs)])/sqrt(length([num_locs1_full(l4_pairs) num_locs1_full(l5_pairs) num_locs2_full(l5_pairs)]))
std(num_locs2_full(l4_pairs))/sqrt(length(num_locs2_full(l4_pairs)))
%%
% close all
for ii = 1:length(z_spike_traces)
    
    map_ch = z_spike_traces{ii};
    
    for i = 1:length(map_ch)
       
        z_detection_results_tmp{i} = detect_peaks(...
                -1.0*bsxfun(@minus,map_ch{i},map_ch{i}(:,1)),20,80,1,1,0,0,1)*70;

    end
    
    z_detection{ii} = z_detection_results_tmp';
    clear z_detection_results_tmp
    figure
    subplot(121)
    plot_trace_stack_grid(map_ch,Inf,1,0);
    subplot(122)
    plot_trace_stack_grid(z_detection{ii},Inf,1,0);
end



%%
z_first_spike_time_trials = {};
centers = [9 11 10 11];
% spike_counts = zeros(15,1,4);
z_first_spike_time = nan(15,4);
z_first_spike_time_trials = cell(15,4);
% jitters = nan(15,1,4);
z_spike_counts = zeros(13,4);
z_spike_trials = zeros(13,4);
z_trial_counts = zeros(15,4);
for ii = 1:length(z_detection)

    
    for i = 1:length(z_detection{ii})
            
        offset_position = abs(i - centers(ii)) + 1;

        z_spike_counts(offset_position,ii) = z_spike_counts(offset_position,ii) + ...
            sum(sum(z_detection{ii}{i}(:,100:end)/70,2));
        z_spike_trials(offset_position,ii) = z_spike_trials(offset_position,ii) + size(z_detection{ii}{i},1);
%             first_spike_time(main_i,main_j,this_exp) = ...
%                 arrayfun(@(x) find(x,1,'first'),detection_grids{this_exp}{i,j}(:,100:end));
        z_trial_counts(offset_position,ii) = z_trial_counts(offset_position,ii) + size(z_detection{ii}{i},1);
        if z_spike_counts(offset_position,ii) > 0
            count = 1;
            for k = 1:size(z_detection{ii}{i},1)
                time = find(z_detection{ii}{i}(k,100:end),1,'first');
                if ~isempty(time)
%                     if count == 1
                    if isempty(z_first_spike_time_trials{offset_position,ii})
% %                         z_first_spike_time(offset_position,ii) = 0;
                        z_first_spike_time_trials{offset_position,ii} = time;
                    else
                        z_first_spike_time_trials{offset_position,ii} = [z_first_spike_time_trials{offset_position,ii} time];
                    end
                    z_first_spike_time(offset_position,ii) = ...
                        z_first_spike_time(offset_position,ii) + time;
                    
                end
                    count = count + 1;
            end
            z_first_spike_time(offset_position,ii) = ...
                z_first_spike_time(offset_position,ii)/count;
            if length(z_first_spike_time_trials{offset_position,ii}) > 4
%                     jitters(main_i,main_j,this_exp) = 2*sqrt(2*log(2))*std(first_spike_time_trials{main_i,main_j,this_exp});
                  z_jitters(offset_position,ii) = max(z_first_spike_time_trials{offset_position,ii}) - min(z_first_spike_time_trials{offset_position,ii}) ;
%                     if jitters(main_i,main_j,this_exp) < 10
%                         return
%                     end
            end
        end

    end
    

end
z_spike_counts_mean = z_spike_counts./z_spike_trials;
z_spike_counts_mean(isnan(z_spike_counts_mean)) = 0;

z_first_spike_time = cellfun(@(x) mean(x),z_first_spike_time_trials);

z_jitters = cellfun(@(x) 2*sqrt(2*log(2))*(sqrt(var(x)+25^2)),z_first_spike_time_trials);

%%
z_spike_counts_thresh = z_spike_counts_mean;
z_spike_counts_thresh(z_spike_counts_thresh > 1) = 1;

figure; plot(0:10:120,mean(z_spike_counts_mean,2))

figure; plot(0:10:120,nanmean(z_first_spike_time(:,[1 2 3 4]),2))

%%

test_jitters = z_jitters(:);
% test_jitters = test_jitters(:);
% spike_counts_mean = z_spike_counts_mean(2:10,:);
% spike_counts_mean = spike_counts_mean(:);
% test_jitters(spike_counts_mean < .5) = NaN;
test_jitters = test_jitters(~isnan(test_jitters));
jitters_red = squeeze(jitters(:));
jitters_red = jitters_red(~isnan(jitters_red));
all_jitters = [test_jitters];
figure; 
subplot(121)
histogram(test_jitters(:)/20,50)
subplot(122)
ecdf(test_jitters(:)/20)
% spike_times_tmp = z_first_spike_time(2:10,2:4);
% spike_counts_mean = z_spike_counts_mean(2:10,2:4);
% spike_times_tmp = spike_times_tmp(spike_counts_mean < .5)
% spike_times_tmp = spike_times_tmp(~isnan(spike_times_tmp));
% spike_times_tmp = spike_times_tmp(:);
% hold on; ecdf(spike_times_tmp(:)/20)
% hold on; plot([5 5],[0 1],'k--')

%%

figure;
subplot(121)
jitter_mean = nanmean(jitters,3);
jitter_mean(isnan(jitter_mean)) = max(jitter_mean(:)) + min(jitter_mean(:));
imagesc(jitter_mean/20)
axis off
colorbar
subplot(122)
errorbar(0:10:40,nanmean(z_jitters(1:5,:)/20,2),std(z_jitters(1:5,:)/20,[],2)/sqrt(size(z_jitters,2)))
xlim([-5 45])

%%

figure;
subplot(121)
jitter_mean = nanmean(first_spike_time,3);
jitter_mean(isnan(jitter_mean)) = max(jitter_mean(:)) + min(jitter_mean(:));
imagesc(jitter_mean/20)
axis off
colorbar
subplot(122)
errorbar(0:10:40,nanmean(z_first_spike_time(1:5,2:4)/20,2),std(z_first_spike_time(1:5,2:4)/20,[],2)/sqrt(size(z_first_spike_time,2)-1))
xlim([-5 45])

%%


for this_exp = 1:4
    
    z_prob_spike(:,this_exp) = cellfun(@(x) length(x),squeeze(z_first_spike_time_trials(1:13,this_exp)))./z_spike_trials(:,this_exp);
end


%%
figure
subplot(211)
x = (0:10:120)';
mean_vec = mean(z_spike_counts_mean,2);
error = bsxfun(@rdivide,squeeze(std(z_spike_counts_mean,[],2)),sqrt(sum(z_spike_trials > 0,2)));  
% fill([x(1:10);flipud(x(1:10))],[(mean_vec(1:10)-error(1:10));flipud(mean_vec(1:10)+error(1:10))],[1 1 1]*.75,'linestyle','none');
% hold on
% plot(x(1:10),mean_vec(1:10),'Color',[0 0 0]/255)
errorbar(x(1:10),mean_vec(1:10),error(1:10));
xlim([0 90])

subplot(212)
x = (0:10:120)';
mean_vec = nanmean(z_prob_spike,2);
error = bsxfun(@rdivide,squeeze(nanstd(z_prob_spike,[],2)),sqrt(sum(z_spike_trials > 0,2)));  
% fill([x(1:10);flipud(x(1:10))],[(mean_vec(1:10)-error(1:10));flipud(mean_vec(1:10)+error(1:10))],[1 1 1]*.75,'linestyle','none');
% hold on
% plot(x(1:10),mean_vec(1:10),'Color',[0 0 0]/255)
errorbar(x(1:10),mean_vec(1:10),error(1:10));
xlim([0 90])

%%

i = 10; j = 9; 
this_exp = 20;

figure; 
subplot(211)
plot_trace_stack(time_posteriors{this_exp}{1}{1}{i,j},4000,zeros(size(time_posteriors{this_exp}{1}{1}{i,j},1),3),'-',new_event_detection{this_exp}{1}{1}{i,j})
subplot(212)
plot_trace_stack(time_posteriors{this_exp}{1}{2}{i,j},4000,zeros(size(time_posteriors{this_exp}{1}{2}{i,j},1),3),'-',new_event_detection{this_exp}{1}{2}{i,j})

%%
exps_to_run = exps
psth_cells = zeros(100,30,2);
grand_psth = zeros(100,1);
count = 0;
for i = 1:length(exps_to_run)
    this_exp = exps_to_run(i)
    for j = 1:2
        for k = 1:numel(new_event_detection_2017{this_exp}{1}{j})
            if input_maps_2017_ww_q10(this_exp).(['ch' num2str(j)])(k) == 1
                for m = 1:length(new_event_detection_2017{this_exp}{1}{j}{k})
                    psth_cells(new_event_detection_2017{this_exp}{1}{j}{k}{m}.times,this_exp,j) = ...
                        psth_cells(new_event_detection_2017{this_exp}{1}{j}{k}{m}.times,this_exp,j) + 1/length(new_event_detection_2017{this_exp}{1}{j}{k});
                    grand_psth(new_event_detection_2017{this_exp}{1}{j}{k}{m}.times) = grand_psth(new_event_detection_2017{this_exp}{1}{j}{k}{m}.times) + 1;
                    count = count + 1;
                end
            end    
        end
    end
end
grand_psth = grand_psth'/count;

% exps_to_run = [1 3 20 28]


full_psth = zeros(100,1);
all_psths = zeros(length(exps_to_run)*2,100);
count = 1;
% figure
for i = 1:length(exps_to_run)
    this_exp = exps_to_run(i);
    full_psth = full_psth + sum(psth_cells(:,this_exp,:),3);
    all_psths(count,:) = psth_cells(:,this_exp,1);
    count = count + 1;
    all_psths(count,:) = psth_cells(:,this_exp,2);
    count = count + 1;
%     plot(psth_cells(:,this_exp,1))
%     hold on
%     plot(psth_cells(:,this_exp,2))
%     hold on
end
% figure;
% plot_trace_stack(all_psths,10,zeros(size(all_psths,1),3),'-')
% full_psth = full_psth/(length(exps_to_run)*2);
mean_vec = mean(all_psths);
error = grand_psth.*(1-grand_psth)/sqrt(length(exps_to_run)*2);
x = 1:length(error);
figure
% plot(2:99,full_psth(2:end-1))
% hold on
fill([x';flipud(x')],[(grand_psth'-error');flipud(grand_psth'+error')],[1 1 1]*.75,'linestyle','none');
hold on
plot(x,grand_psth,'Color',[0 0 0]/255)
xlim([1 100])
% errorbar(mean(all_psths),error)
%%

this_exp = 6;
i = 6; j = 10; 
figure; 
subplot(121); 
plot_trace_stack(time_posteriors_2017{this_exp}{1}{1}{i,j},4000,zeros(size(time_posteriors_2017{this_exp}{1}{1}{i,j},1),3),'-',new_event_detection_2017{this_exp}{1}{1}{i,j}); 
subplot(122); 
plot_trace_stack(time_posteriors_2017{this_exp}{1}{2}{i,j},4000,zeros(size(time_posteriors_2017{this_exp}{1}{2}{i,j},1),3),'-',new_event_detection_2017{this_exp}{1}{2}{i,j})

%%

this_exp = 6;
i = 6; j = 10; 
figure; 
subplot(121); 
plot_trace_stack(map_ch1{i,j},300,zeros(size(time_posteriors_2017{this_exp}{1}{1}{i,j},1),3),'-',new_event_detection_2017{this_exp}{1}{1}{i,j}); 
subplot(122); 
plot_trace_stack(map_ch2{i,j},300,zeros(size(time_posteriors_2017{this_exp}{1}{2}{i,j},1),3),'-',new_event_detection_2017{this_exp}{1}{2}{i,j})


%%

i = 28;
sweeps_window = [3000 5000];
time_window = [0 Inf];
amps_window = [0 Inf];

matchstr =  [dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
    num2str(cell_nums(i)) tags{i} '_ch1']
rebuildmap_file = ...
    ['old-data/' dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
    num2str(cell_nums(i)) tags{i} '_ch1_trace_grid-rebuildmap.mat'];
resultsdir = '/media/shababo/data/new-data';
results_grid_ch1 = rebuild_detection_grid(...
resultsdir,matchstr,...
rebuildmap_file);


results_grid_trunc_ch1 = ...
cellfun(@(x) arrayfun(@(y) ...
truncate_samples(y,sweeps_window,time_window,amps_window),x),...
results_grid_ch1,'UniformOutput',0);

matchstr =  [dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
    num2str(cell_nums(i)) tags{i} '_ch2']
rebuildmap_file = ...
    ['old-data/' dates{i} '_slice' num2str(slice_nums(i)) '_cell' ...
    num2str(cell_nums(i)) tags{i} '_ch2_trace_grid-rebuildmap.mat'];
resultsdir = '/media/shababo/data/new-data';
results_grid_ch2 = rebuild_detection_grid(...
resultsdir,matchstr,...
rebuildmap_file);


results_grid_trunc_ch2 = ...
cellfun(@(x) arrayfun(@(y) ...
truncate_samples(y,sweeps_window,time_window,amps_window),x),...
results_grid_ch2,'UniformOutput',0);

detection_viewer(map_ch1,results_grid_trunc_ch1);
detection_viewer(map_ch2,results_grid_trunc_ch2);

%%
index_key = spots_key/20 + 9;
all_targets(end,:) = [];
all_indices = zeros(size(all_targets,2),2,size(all_targets,1));

for i = 1:size(all_targets,1)
    
    all_indices(:,:,i) = index_key(all_targets(i,:),1:2);

end

%%


all_detect_ests = {detect_ests_t19,detect_ests_t18,detect_ests_t20};
all_mean_events = {mean_events_image_t19,mean_events_image_t18,mean_events_image_t20};
mean_events_thresh = zeros([size(detect_ests_t19) length(all_detect_ests)]);

mean_amp_thresh = ones([size(mean_events_thresh)]);
threshes = [.45 1.2 .4];

for k = 1:length(all_detect_ests)
    tmp = zeros(size(all_detect_ests{k}));
    tmp(all_mean_events{k} > threshes(k)) = all_mean_events{k}(all_mean_events{k} > threshes(k));
    mean_events_thresh(:,:,k) = tmp;
    for i = 1:size(mean_events_thresh,1)
        for j = 1:size(mean_events_thresh,2)
            if ~isempty(all_detect_ests{k}{i,j})
                mean_amp_thresh(i,j,k) = mean([all_detect_ests{k}{i,j}{1}.amp]);
                if isnan(mean_amp_thresh(i,j,k))
                    mean_amp_thresh(i,j,k) = 1;
                end
            end
        end
    end
end



[X,Y,Z] = meshgrid(-160:20:160,-160:20:160,(50:-50:-50)*2);
cmap = hot(1000);

figure; scatter3(X(:),Y(:),Z(:),5,[0 0 0],'filled'); axis off; axis image
% [X,Y] = meshgrid(-160:20:160,-160:20:160);
hold on; scatter3(X(:),Y(:),Z(:),sign(mean_events_thresh(:))*100+.001,cmap(round(mean_amp_thresh(:)/(max(mean_amp_thresh(:))+10)*1000),:),'filled')





%%

vb_ed_3D_doubleAx = ROI_VB_3D_grid(detect_ests_3D_fix8,all_indices_3D,[zeros(1,1000) 50*ones(1,1000) 100*ones(1,1000)],[]);


map1_ed_doubleAx = reshape(vb_ed_3D_doubleAx.alpha(1:289),17,17);
map2_ed_doubleAx= reshape(vb_ed_3D_doubleAx.alpha(1+289:289+289),17,17);
map3_ed_doubleAx = reshape(vb_ed_3D_doubleAx.alpha(1+289*2:289+289*2),17,17);
figure; 
subplot(321); imagesc(map1_ed_doubleA); title('z = 0 um'); caxis([0 1]); axis off
subplot(322);imagesc(map1_ed_doubleAx); title('z = 0 um');   caxis([0 1]); axis off
subplot(323); imagesc(map2_ed_doubleA); title('z = 50 um');  caxis([0 1]); axis off
subplot(324);imagesc(map2_ed_doubleAx); title('z = 50 um');  caxis([0 1]); axis off
subplot(325); imagesc(map3_ed_doubleA); title('z = 100 um');  caxis([0 1]); axis off
subplot(326);imagesc(map3_ed_doubleAx); title('z = 100 um'); caxis([0 1]); axis off

%%
% close all
% all_trials = [8];
% trial_names = {'z0','z40','z80','z120'};
% base_name = 'data/04142017_s2c1_%dmw_%s_traces';
% for j = 1:length(all_trials)
%     trials = all_trials(j)
%     
%     this_seq = cell(length(trials),1);
%     this_stim_key = cell(length(trials),1);
%     power_curve_num = cell(length(trials),1);
%     stims_per_trial = zeros(length(trials),1);
%     for i = 1:length(trials)
%         cur_trial = trials(i);
%         this_seq{i} = data.trial_metadata(cur_trial).sequence;
%         stims_per_trial(i) = length(this_seq{i});
%         this_stim_key{i} = data.trial_metadata(cur_trial).stim_key;
%         power_curve_num{i} = unique([this_seq{i}.target_power]);
%     end

% <<<<<<< HEAD
close all
trials = [5];
this_seq = cell(length(trials),1);
this_stim_key = cell(length(trials),1);
power_curve_num = cell(length(trials),1);

for i = 1:length(trials)
    cur_trial = trials(i);
    this_seq{i} = data.trial_metadata(cur_trial).sequence;
    stims_per_trial(i) = length(this_seq{i});
    this_stim_key{i} = data.trial_metadata(cur_trial).stim_key;
    power_curve_num{i} = unique([this_seq{i}.target_power]);
end
power_curve_num = unique([power_curve_num{:}]);
maps = cell(length(power_curve_num),1);
this_seq = [this_seq{:}];
max_trial = length(this_seq);
% max_trial = 1200;
for i = 1:length(power_curve_num)
    [traces_ch1,traces_ch2] = ...
    get_stim_stack(data,trials,...
        stims_per_trial,[this_seq.start]);
    this_seq = this_seq(1:max_trial);
    traces_pow{1} = traces_ch1([this_seq.target_power] == power_curve_num(i),:);
    traces_pow{2} = traces_ch2([this_seq.target_power] == power_curve_num(i),:);
    this_seq_power = this_seq([this_seq.target_power] == power_curve_num(i));
    [maps{i}, map_index] = see_grid_multi(traces_pow,[],this_seq_power,this_stim_key{1},5,1);
    title(['Power = ' num2str(power_curve_num(i)) ' mW'])
% =======
%     power_curve_num = unique([power_curve_num{:}]);
%     power_curve_num = [25 75];
%     maps = cell(length(power_curve_num),1);
%     this_seq = [this_seq{:}];
%     clear traces
%     for i = 1:length(power_curve_num)
%         [traces_ch1,traces_ch2] = ...
%         get_stim_stack(data,trials,...
%             sum(stims_per_trial));
% 
%         traces_pow{1} = traces_ch1([this_seq.target_power] == power_curve_num(i),:);
%         traces_pow{2} = traces_ch2([this_seq.target_power] == power_curve_num(i),:);
%         this_seq_power = this_seq([this_seq.target_power] == power_curve_num(i));
%     %     [maps{i}, map_index] = see_grid_multi(traces_pow,this_seq_power,this_stim_key{1},15,1);
%         name = sprintf(base_name, power_curve_num(i),trial_names{j});
% %         build_trace_param_files(traces_pow{1},get_params(),'/media/shababo/data/04142017_s2c1_traces','/media/shababo/data/04142017_s2c1_params',name);
%     %     title(['Power = ' num2str(power_curve_num(i)) ' mW'])
%     end
% >>>>>>> origin/master
end
%%
output = ROI_VB_3D_charge(map_index([this_seq_power.precomputed_target_index],:,:),zeros(1000,1),-traces_pow{1});
figure; imagesc(reshape(output.alpha,33,33))
%%
figure;
subplot(121)
plot_trace_stack(maps{1}{1}{26,4},40,'k')
subplot(122)
plot_trace_stack(maps{2}{1}{26,4},40,'k')

%%
x = 12; y = 21;
map_index_trial = map_index([this_seq_power.precomputed_target_index],:,:);
trials = find((map_index_trial(:,1,1) == x & map_index_trial(:,2,1) == y) | ...
    (map_index_trial(:,1,2) == x & map_index_trial(:,2,2) == y) | ...
    (map_index_trial(:,1,3) == x & map_index_trial(:,2,3) == y))
figure; plot(traces_pow{1}([trials],:)')
%%
figure; plot(traces_cent([trials],:)')
sum(traces_cent(trials,140:800),2)
Y_n(trials)

%%

subset_inds = randi(size(traces,1),20,1);
figure; plot_trace_stack(traces(subset_inds,:),100,'-',events_map(subset_inds))

%% plot nuclear detect with map

% close all
trials = 3:10;
% trials = 4;

this_seq = cell(length(trials),1);
this_stim_key = cell(length(trials),1);
power_curve_num = cell(length(trials),1);
stim_starts = cell(length(trials),1);

full_stim_key = [];
clear full_seq
% full_seq = struct();
for i = 1:length(trials)
    cur_trial = trials(i);
    this_seq{i} = data.trial_metadata(cur_trial).sequence;
    stims_per_trial(i) = length(this_seq{i});
    this_stim_key{i} = data.trial_metadata(cur_trial).stim_key;
    power_curve_num{i} = unique([this_seq{i}.target_power]);
    stim_starts{i} = [data.trial_metadata(cur_trial).sequence.start];
    for j = 1:length(this_seq{i})
        if i == 1 && j == 1
            full_seq(1) = this_seq{i}(j);
        else
            full_seq(end+1) = this_seq{i}(j);
        end
        full_seq(end).precomputed_target_index = ...
            full_seq(end).precomputed_target_index + size(full_stim_key,1);
    end
    full_stim_key = [full_stim_key; this_stim_key{i}];
end
power_curve_num = unique([power_curve_num{:}]);
maps = cell(length(power_curve_num),1);
this_seq = [this_seq{:}];
max_trial = length(this_seq);
% max_trial = 1200;
[traces_ch1,traces_ch2] = ...
    get_stim_stack(data,trials,...
        stims_per_trial,stim_starts);
stim_inds = [full_seq.precomputed_target_index];
on_cell_trials = isnan(full_stim_key(stim_inds,1,2));
on_cell_trials = ones(size(on_cell_trials));
% power_curve_num = 100;
traces = [];
stim_pow = [];
target_locs = [];
stim_inds = [];
deorder = [];
num_trials = 0;
for i = 1:length(power_curve_num)
    
    
%     this_seq = this_seq(1:max_trial);
    traces_pow{1} = traces_ch1(on_cell_trials' & [full_seq.target_power] == power_curve_num(i),:);
%     traces = [traces; traces_pow{1}];
    deorder = [deorder find(on_cell_trials' & [full_seq.target_power] == power_curve_num(i))]; 
    traces_pow{2} = traces_ch2(on_cell_trials' & [full_seq.target_power] == power_curve_num(i),:);
    this_seq_power = full_seq(on_cell_trials' & [full_seq.target_power] == power_curve_num(i));
%     mpp_pow = mpp(on_cell_trials' & [full_seq.target_power] == power_curve_num(i));
    mpp_pow = mpp(num_trials+(1:length(this_seq_power)));
    num_trials = num_trials + length(this_seq_power);
    [maps{i}, mpp_maps{i}] = see_grid_multi(traces_pow,mpp_pow,this_seq_power,full_stim_key,5,1);
    title(['Power = ' num2str(power_curve_num(i)) ' mW'])
    get_mean_events_image(mpp_maps{i}, 2000, 1, 1);
    title(['Event Counts, Power = ' num2str(power_curve_num(i)) ' mW'])
end















