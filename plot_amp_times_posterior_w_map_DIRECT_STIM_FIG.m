function time_amp_posteriors = plot_amp_times_posterior_w_map_DIRECT_STIM_FIG(results_file, trace_offset, bin_edges, varargin)

load(results_file)

try
    load(params.traces_filename)
catch
    [pathname, filename] = fileparts(params.traces_filename);
    pathname = '/media/shababo/Layover/projects/mapping/code/psc-detection/data/for-paper/';
    load([pathname filename])
end

% load('/home/shababo/Projects/Mapping/code/psc-detection/data/simulated-data-longer-traces-epsc.mat')
% load('/home/shababo/Desktop/simulated-data-longer-traces-epsc.mat')

if ~isfield(params,'start_ind')
    params.start_ind = 1;
end

if ~isfield(params,'duration')
    params.duration = size(traces,2);
end
params.duration = 2700;
traces = traces(:,params.start_ind:(params.start_ind + params.duration - 1));

if isfield(params,'traces_ind')
    traces = traces(params.traces_ind,:);
    
end

if ~isempty(varargin) && ~isempty(varargin{1})
    traces_ind = varargin{1};
    traces = traces(traces_ind,:);
    if exist('true_signal','var')
        true_signal = true_signal(traces_ind,:);
    end
    results = results(traces_ind);
end

if length(varargin) >= 2 && ~isempty(varargin{2})
    max_sample = varargin{2};
end


[num_traces, T] = size(traces);
num_amp_bins = length(bin_edges);

if isfield(params,'event_samples')
    event_samples = params.event_samples;
else
    event_samples = 4000;
end

time_amp_posteriors = zeros(size(traces,1),length(bin_edges),T);
direct_stim_exp_val = zeros(1,T);

stim_start = 500;
stim_duration = .05*20000;
stim_in = [zeros(1,stim_start) ones(1,stim_duration) zeros(1,T-stim_start-stim_duration)];
t = 0:T-1;

burn_in = 500;
for i = 1:size(traces,1)
    
    if ~exist('max_sample','var')
        map_i = results(i).map_ind;
    else
        [map_val, map_i] = max(results(i).trials.obj(1:max_sample));
    end
    this_time_amp_posterior = zeros(length(bin_edges),T);
    
    for j = 1:length(results(i).trials.times)
        
        for k = 1:length(results(i).trials.times{j})
            this_time_ind = ceil(results(i).trials.times{j}(k));
            this_amp_ind = find(results(i).trials.amp{j}(k) < bin_edges,1,'first');
            
            if this_time_ind <= T
                this_time_amp_posterior(this_amp_ind,this_time_ind) =...
                    this_time_amp_posterior(this_amp_ind,this_time_ind) + 1/length(results(i).trials.times);
            end
            
        end
        
        if j > burn_in
            stim_tau_rise = results(i).trials.stim_tau_rise{j}; % values for chr2 from lin et al 2009 (biophysics)
            stim_tau_fall = results(i).trials.stim_tau_fall{j};
            stim_amp = results(i).trials.stim_amp{j};

            stim_decay = exp(-t/stim_tau_fall);
            stim_rise = -exp(-t/stim_tau_rise);
            stim_kernel = (stim_decay + stim_rise)/sum(stim_decay + stim_rise);
            stim_response = conv(stim_in,stim_kernel);
            direct_stim_exp_val = direct_stim_exp_val + stim_amp*stim_response(1:T)/max(stim_response(1:T));
        end
        
        
    end
    
    direct_stim_exp_val = direct_stim_exp_val/(length(results(i).trials.times)-burn_in);
    
    time_amp_posteriors(i,:,:) = this_time_amp_posterior;
    this_curve = zeros(1,T);
    for j = 1:length(results(i).trials.times{map_i})
        
%         if results(i).trials.tau{map_i}{j}(2) < 300
        
            ef = genEfilt_ar(results(i).trials.tau{map_i}{j},event_samples);
            [~, this_curve, ~] = addSpike_ar(results(i).trials.times{map_i}(j),...
                                                this_curve, 0, ef, ...
                                                results(i).trials.amp{map_i}(j),...
                                                results(i).trials.tau{map_i}{j},...
                                                traces(i,:), ...
                                                results(i).trials.times{map_i}(j), ...
                                                2, 1, 1);
%         end
                                        
    end
    map_curves(i,:) = this_curve;% + results(i).trials.base{map_i};
end

ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[.3 .1 .6 .8]);

descr = {'Parameters:'
    ['a_{min} = ' num2str(params.a_min)];
    ['tau^1_{min} = ' num2str(params.tau1_min)];
    ['tau^1_{max} = ' num2str(params.tau1_max)];
    ['tau^2_{min} = ' num2str(params.tau2_min)];
    ['tau^2_{max} = ' num2str(params.tau2_max)];
    ['p_{spike} = ' num2str(params.p_spike)]
    ['num sweeps = ' num2str(params.num_sweeps)]};


axes(ax1) % sets ax1 to current axes
text(.025,0.6,descr)

axes(ax2)
% plot_trace_stack(traces,trace_offset,bsxfun(@plus,zeros(length(traces),3),[1 .4 .4]),'-',[.005 25],0)
% hold on

offset = 140;
plot_trace_stack(traces,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 0]),'-',[],offset)
hold on
plot_scatter_stack(time_amp_posteriors,trace_offset,bin_edges,-35,2000)
hold on
if exist('true_signal','var')
    disp('...')
    times_vec = zeros(1,size(time_amp_posteriors,3));
%     times_vec(ceil(true_event_times{1})) = max(max(max(time_amp_posteriors)))+.1;

%     hold on
    plot_trace_stack(true_signal(:,1:T),trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 0]),'-',[.01 25],0,2)
    hold on
    
    scatter(true_event_times{1}(find(true_event_times{1} < params.duration))/20000,-true_amplitudes{1}(find(true_event_times{1} < params.duration))-35,'xr','LineWidth',3)
    hold on
end

stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
stim_tau_fall = .013*20000;
stim_amp = 50;
stim_start = 500;
stim_duration = .05*20000;
stim_in = [zeros(1,stim_start) ones(1,stim_duration) zeros(1,T-stim_start-stim_duration)];
t = 0:T-1;
stim_decay = exp(-t/stim_tau_fall);
stim_rise = -exp(-t/stim_tau_rise);
stim_kernel = (stim_decay + stim_rise)/sum(stim_decay + stim_rise);
stim_response = conv(stim_in,stim_kernel);
stim_response = stim_amp*stim_response(1:T)/max(stim_response(1:T));
plot((0:length(stim_response)-1)/20000,-stim_response+offset/2,'r-','linewidth',2)
plot((0:length(stim_response)-1)/20000,-stim_response+offset/2,'b--','linewidth',2)


% plot_trace_stack(params.event_sign*map_curves,trace_offset,bsxfun(@plus,zeros(length(traces),3),[1 0 0]),'--',[.01 5],0,1)
hold on

leg_dummy_1 = plot(0,0,'k','LineWidth',2,'Visible','off');
leg_dummy_2 = scatter(0,0,'ob','filled','Visible','off');
leg_dummy_3 = plot(0,0,'--k','LineWidth',2,'Visible','off');
leg_dummy_4 = scatter(0,0,'xr','linewidth',3,'Visible','off');
legend([leg_dummy_1, leg_dummy_2, leg_dummy_3, leg_dummy_4],{'observation','Bayesian time-amplitude posterior','true current','true time-amplitude coordinate'})


hold off

title(strrep(results_file,'_','-'))

if length(varargin) > 2 && ~isempty(varargin{3})
    [dir,name,~] = fileparts(results_file);
    savefig([dir '/' name '.fig'])
end

% map = params.event_sign*map_curves;