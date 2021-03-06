function plot_trace_stack(traces, offset_step, linespec, varargin)


offset = 0;
stim_start = 0;
time_after_stim = 1; %1000

trial_length = size(traces,2)

for trial = 1:size(traces,1)
    
    offset = offset - 2 * min(traces(trial,:));
    
end

stim_top = 2*max(traces(1,:));
stim_bottom = -offset;

% num_stims = length(find(diff(stims(1,:))))/2;
% change_points = find(diff(stims(1,:)));

% if ~isempty(change_points)
% 
%     trial_length = change_points(end) - change_points(1) + stim_start + time_after_stim;
% 
%     for i = 1:num_stims
% 
%         stim_length = change_points(2*i) - change_points(2*i - 1);
%         this_start = change_points(2*i-1)- change_points(1) + stim_start;
%         rectangle('Position', [this_start stim_bottom stim_length stim_top-stim_bottom],'FaceColor','b','EdgeColor','b')
%         hold on
%     end
% else
%     trial_length = default_length;
% 
% end

offset = 0;

if ~isempty(varargin) && ~isempty(varargin{1})
    events = varargin{1};
% else
%     events = [];
end
if length(varargin) > 1 && ~isempty(varargin{2})
    vert_offset = varargin{2};
else
    vert_offset = 0;
end

if length(varargin) > 2 && ~isempty(varargin{3})
    linewidth = varargin{3};
else
    linewidth = 2;
end
    
if length(varargin) > 3 && ~isempty(varargin{4})
    plot_mean = varargin{4};
else
    plot_mean = 0;
end

for trial = 1:size(traces,1)
    
%     if isempty(change_points)
%         this_trial_start = 1;
%     else
%         this_trial_start = find(stims(trial,:),1,'first') - stim_start;
%     end
    time = (0:trial_length-1);
    trace_to_plot = traces(trial,:);
    %median(trace_to_plot)
    plot(time,trace_to_plot - offset - trace_to_plot(1) + vert_offset,linespec,'LineWidth',linewidth,'color',[0 0 0])
    hold on
    
    
    
%     if ~isempty(events)
%         these_events.times(these_events.times < 40) = [];
%         if iscell(events)
%             these_events = events{trial};
%         else
%             these_events = events(trial);
%         end
%         these_events.times = ceil(these_events.times-1);
%         events_mat = [these_events.amp' these_events.tau1' these_events.tau2' these_events.times'];
%         denoised_curve = build_curve(events_mat,0,trial_length,1,length(trace_to_plot));
%         size(denoised_curve)
%         size(trace_to_plot)
%         plot((0:trial_length-1),-denoised_curve - offset + vert_offset,linespec,'LineWidth',linewidth)
%         scatter((these_events.times - stim_start),(max(trace_to_plot) - offset - trace_to_plot(1) + vert_offset + offset_step/5)*ones(size(these_events.times)),[],[0 0 0],'filled')
%         hold on
%     end

    if exist('events','var')
        these_events = events(trial);
        event_times = [];
        event_pos = [];
        if ~ isempty(these_events)
%                     length(offsets)
%             for ii = 1:length(these_events)
%                 if iscell(these_events)
%                     this_struct = these_events{ii};
%                 else
%                     this_struct = these_events(ii);                        
%                 end
                event_times = these_events.times;
                event_pos = (-offset - trace_to_plot(1) + vert_offset)*ones(size(these_events.times));
%             end
%         event_times = these_events.times;
        event_pos(event_times > length(time)-20) = [];
        event_times(event_times > length(time)-20) = [];
        event_pos(event_times < 10) = [];
        event_times(event_times < 10) = [];
%         event_pos
%         event_times
        scatter(time(round(event_times)),event_pos+40,20*ones(size(event_pos)),[.0 .5 1],'filled')
        end
        hold on;
    end
    
    offset = offset + offset_step;
    
    
end

if plot_mean
    meantrace = mean(traces);
    plot(time,meantrace - meantrace(1),linespec,'LineWidth',4,'color',[0 0 0]);
end
% 
% if ~isempty(varargin)
%     bar_limits = varargin{1};
% 
%     if ~isempty(bar_limits)
% 
%         bar_corner_time = trial_length/10;
%         bar_corner_y = -offset + vert_offset;
% 
%         plot([bar_corner_time; bar_corner_time], bar_corner_y + [0; bar_limits(2)], '-k',  bar_corner_time + [0; bar_limits(1)], [bar_corner_y; bar_corner_y], '-k', 'LineWidth', 2)
%         text(bar_corner_time - bar_limits(1)/2,bar_corner_y + bar_limits(2)/2, [num2str(bar_limits(2)) ' pA'], 'HorizontalAlignment','right')
%         text(bar_corner_time + bar_limits(1)/2,bar_corner_y - bar_limits(2)/2, [num2str(bar_limits(1)*1000) ' ms'], 'HorizontalAlignment','center')
%     end
% end

axis tight
axis off

hold off

