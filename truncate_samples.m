function posterior_trunc = truncate_samples(posterior,sweep_bounds,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    time_bounds = varargin{1};
else
    time_bounds = [0 Inf];
end

if length(varargin) > 1 && ~isempty(varargin{2})
    amp_bounds = varargin{2};
else
    amp_bounds = [-Inf Inf];
end


posterior_trunc = struct();
% num_new_events = sum(posterior.num_events(sweep_bounds(1):sweep_bounds(2)));
event_sample_start = max(sum(posterior.num_events(1:sweep_bounds(1)-1)) + 1,1);
event_sample_end = sum(posterior.num_events(1:sweep_bounds(2)));
if event_sample_end > 0
    posterior_trunc.amp = posterior.amp(event_sample_start:event_sample_end);
else
    posterior_trunc.amp = [];
end
posterior_trunc.base = posterior.base(sweep_bounds(1):sweep_bounds(2));
if event_sample_end > 0
    posterior_trunc.tau1 = posterior.tau1(event_sample_start:event_sample_end);
    posterior_trunc.tau2 = posterior.tau2(event_sample_start:event_sample_end);
else
    posterior_trunc.tau1 = [];
    posterior_trunc.tau2 = [];
end
posterior_trunc.num_events = posterior.num_events(sweep_bounds(1):sweep_bounds(2));
posterior_trunc.phi = posterior.phi(sweep_bounds(1):sweep_bounds(2),:);
posterior_trunc.noise = posterior.noise(sweep_bounds(1):sweep_bounds(2));
posterior_trunc.obj = posterior.obj(sweep_bounds(1):sweep_bounds(2));
if event_sample_end > 0
    posterior_trunc.times = posterior.times(event_sample_start:event_sample_end);
else
    posterior_trunc.times = [];
end

sweep_count = 1;
event_count = 1;
bad_samps = [];
num_events_tmp = posterior_trunc.num_events;


for i = 1:length(posterior_trunc.times)
    
%     i
    

    while sweep_count <= length(posterior_trunc.num_events) && ...
            event_count > posterior_trunc.num_events(sweep_count)
        sweep_count = sweep_count + 1;
        event_count = 1;
    end
    
%     event_count
%     sweep_count
%      sum(posterior_trunc.num_events(1:sweep_count)) 
%     posterior_trunc.num_events(sweep_count)
%     
%     posterior_trunc.times(i)
%     posterior_trunc.amp(i)
    
    if posterior_trunc.times(i) < time_bounds(1) || posterior_trunc.times(i) > time_bounds(2) || ...
            posterior_trunc.amp(i) < amp_bounds(1) || posterior_trunc.amp(i) > amp_bounds(2)
        bad_samps = [bad_samps i];
        num_events_tmp(sweep_count) = num_events_tmp(sweep_count) - 1;
    end
    
    event_count = event_count + 1;
    
    
    
end

% length(bad_samps)
posterior_trunc.num_events = num_events_tmp;
posterior_trunc.times(bad_samps) = [];
posterior_trunc.amp(bad_samps) = [];
posterior_trunc.tau1(bad_samps) = [];
posterior_trunc.tau2(bad_samps) = [];

        