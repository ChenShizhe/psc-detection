function [target_feature_mat, trace_times, all_times] = plot_detection_results(results,var_names)

traces = [];
event_times = [];
target_feature_mat = [];

min_tau1 = 0;
min_tau2 = 0;
max_tau1 = Inf;
max_tau2 = Inf;


min_amp = -Inf;
min_time = 0;
max_time = Inf;%.015*20000;
hot_time_min = .0095;
hot_time_max = .012;
% hot_times = ([.0043 .0177 .0324 .0467 .0610 .0754 .0896 .1040 .1187 .1333 .1470 .1612 .1761 .1891 .2046])*20000;
% hot_time_min = .0047;
% hot_time_max = .007;
hot_amp_min = -Inf;
% hot_amp_min = 20;
hot_amp_max = Inf;
hot_min_tau2 = 0;%50;

% results = results(1:end-1);

for j = 1:length(results)
    
%     if ~any(j == [18 21 31 32])
        
%         if isfield(results(j).trials,'curves')
%             traces(j,:) = results(j).trials.curves{results(j).map_ind}(min_time:end);
%         else
%             traces(j,:) = build_curve(results,j, .05*20000);
%         end
        
        taus = zeros(0,2);
        time_constants = zeros(0,2);
        good_events = [];
        
        for k = 1:length(results(j).trials.tau{results(j).map_ind})
            if ~isempty(results(j).trials.tau{results(j).map_ind}{k})
%                 if all(results(j).trials.tau{results(j).map_ind}{k} > [min_tau1 min_tau2]) && all(results(j).trials.tau{results(j).map_ind}{k} < [max_tau1 max_tau2]) && results(j).trials.times{results(j).map_ind}(k) > min_time && results(j).trials.times{results(j).map_ind}(k) < max_time && results(j).trials.amp{results(j).map_ind}(k) > min_amp
%                 if all(results(j).trials.tau{results(j).map_ind}{k} > [min_tau1 min_tau2]) && all(results(j).trials.tau{results(j).map_ind}{k} < [max_tau1 max_tau2]) && any(abs(hot_times - results(j).trials.times{results(j).map_ind}(k)) < 25) && results(j).trials.amp{results(j).map_ind}(k) > min_amp
                    good_events = [good_events k];
                    taus = [taus; results(j).trials.tau{results(j).map_ind}{k}];
                    
                    tau = results(j).trials.tau{results(j).map_ind}{k};
                    this_curve = zeros(1,800);
                    ef = genEfilt_ar(tau,400);
                    [~, this_curve, ~] = addSpike_ar(1000,...
                        this_curve, 0, ef, ...
                        1,...
                        tau,...
                        zeros(size(this_curve)), ...
                        1, ...
                        2, 1, 1);
                    

                    rise_time = find(this_curve > .67*max(this_curve),1,'first');
                    [~,this_max_time] = max(this_curve); 
                    decay_time = find(this_curve(rise_time:end) < .33*max(this_curve),1,'first') + rise_time - this_max_time;
                    
                    time_constants = [time_constants; rise_time decay_time];% this_max_time];
                    
                    
%                 end
            end
        end
        
        
        times = results(j).trials.times{results(j).map_ind}(good_events);
        event_times = [event_times times];
        

        new_features = [results(j).trials.amp{results(j).map_ind}(good_events)' taus results(j).trials.times{results(j).map_ind}(good_events)'/20000];
        target_feature_mat = [target_feature_mat; new_features];

end

size(target_feature_mat)
groups = [];

for i = 1:size(target_feature_mat,1)
    
    this_amp = target_feature_mat(i,1);
    this_time = target_feature_mat(i,end);
    if this_time > hot_time_min && this_time < hot_time_max && ...
            this_amp > hot_amp_min && this_amp < hot_amp_max && ...
            target_feature_mat(i,3) > hot_min_tau2
        disp('hit')
        groups = [groups 2];
    else
        groups = [groups 1];
    end
end

figure
gplotmatrix(target_feature_mat,target_feature_mat,groups)
% 
figure
plotmatrix(target_feature_mat)

%%

all_times = [];
trace_times = cell(length(results),1);
for i = 1:length(results)
    trace_times{i} = [];
    this_trial_results = results(i).trials.times;
    for j = 1:length(this_trial_results)
        trace_times{i} = [trace_times{i} this_trial_results{j}];
        all_times = [all_times this_trial_results{j}];
    end  
end
% 
% figure;
% for i = 1:length(results)
%     subplot(length(results),1,i)
%     bar(hist(trace_times{i},100:100:2100) / results(1).params.nsweeps)
%     ylim([0 1])
%     axis off
% end
% 
% figure;
% hist(all_times,20)

