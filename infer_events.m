function infer_events(params)

if params.rand == 1
    rng(params.seed)
end

if params.cluster
    addpath('/vega/stats/users/bms2156/psc-detection/functions')
end


if ~isfield(params,'start_ind')
    params.start_ind = 1;
end

try
    load(params.traces_filename,'traces')
catch e
    disp('bad file')
    params.traces_filename
    params.traces_filename = 'data/for-paper/direct-stim-w-events-real.mat';
    load('data/for-paper/direct-stim-w-events-real.mat')
end



% the file were your traces are, traces should be saved in this mat file in
% an N x T matrix called 'traces' where N = number of traces and T = number
% of samples
% load(params.traces_filename,'traces');

if iscell(traces)
    params.is_grid = 1;
else
    params.is_grid = 0;
end

if params.is_grid
    
    if isfield(params,'grid_reduce') && params.grid_reduce
        edge_num = params.grid_reduce_count;
        traces_reduce = cell(edge_num);
        i_picks = randsample(1:size(traces,1),edge_num);
        j_picks = randsample(1:size(traces,2),edge_num);
        for i = 1:edge_num
            for j = 1:edge_num
                traces_reduce{i,j} = traces{i_picks(i),j_picks(j)};
            end
        end
        traces = traces_reduce;
    end 
    [traces, rebuild_map] = stack_traces(traces);
    params.rebuild_map = rebuild_map;
end

% assignin('base','params',params)
% return

if ~isfield(params,'duration')
    params.duration = size(traces,2);
end

% grab section in time
traces = traces(:,params.start_ind:(params.start_ind + params.duration - 1));

% grab subset of traces
if isfield(params,'traces_ind')
    traces = traces(params.traces_ind,:);
end

if params.tau1_min >= params.tau1_max || params.tau2_min >= params.tau2_max
    results = 'infeasible parameter set';
    params.savename = [params.savename(1:end-4) '-z.mat'];
    save(params.savename,'results')
    return
end

results = struct();
disp(['About to run inference on: ' num2str(size(traces,1)) ' traces...']);

if params.par
    
    if ~params.cluster
        delete(gcp('nocreate'))
        this_pool = parpool();
    else
        jobdir = ['/vega/stats/users/bms2156/matlab-pool/' params.jobid];
        mkdir(jobdir)
        pc = parcluster('local'); 
	pc.JobStorageLocation = jobdir;
        this_pool = parpool(pc,16);
        
    end
    addAttachedFiles(this_pool,{'functions/sampleParams_ARnoise_splittau.m',...
				'functions/add_base_ar.m','functions/addSpike_ar.m',...
				'functions/genEfilt_ar.m','functions/predAR.m',...
				'functions/remove_base_ar.m','functions/removeSpike_ar.m'});
    
    load_struct = load(params.init_method.template_file);
    template = load_struct.template;
    if params.event_sign == -1.0
        template = -template;
    end
    
    parfor trace_ind = 1:size(traces,1)
    %     
        disp(['Starting trace #' num2str(trace_ind)])
        trace = params.event_sign*traces(trace_ind,:);
        trace = trace - min(trace);

%         event_times_init_old = template_matching(-1*params.event_sign*traces(trace_ind,:), params.dt,...
%             params.init_method.tau, params.init_method.amp_thresh, params.init_method.conv_thresh);
        
        
        
        nfft = length(trace) + length(template);
        
        [filtered_trace, event_times_init,event_sizes_init] = wiener_filter(trace,template,params.init_method.ar_noise_params,...
           nfft, params.dt, params.init_method.theshold, params.init_method.min_interval);
%         [event_sizes_init, event_times_init] = findpeaks(trace,'minpeakheight',30,'minpeakprominence',10);
        event_times_init
        event_sizes_init
        %event_times_init = [];
%         filtered_trace = [];
        %event_sizes_init = [];

%         assignin('base','event_times_init_old',event_times_init_old)
%         assignin('base','event_times_init',event_times_init)
        results(trace_ind).event_times_init = event_times_init;
        results(trace_ind).filtered_trace = filtered_trace;
        results(trace_ind).event_sizes_init = event_sizes_init;
        
        tau = [mean([params.tau1_min params.tau1_max]) ...
            mean([params.tau2_min params.tau2_max])]/params.dt;
        
        if isfield(params,'init_only') && ~params.init_only
            if params.direct_stim
    %             event_times_init = ceil(length(trace)*rand(1,length(trace)*params.p_spike));
                [results(trace_ind).trials, results(trace_ind).mcmc]  = ...
                    sampleParams_ar_2taus_directstim(trace,tau,event_times_init,params);
            else

    %             event_times_init = template_matching(-1*params.event_sign*traces(trace_ind,:), params.dt,...
    %                 params.init_method.tau, params.init_method.amp_thresh, params.init_method.conv_thresh);
                [results(trace_ind).trials, results(trace_ind).mcmc]  = ...
                    sampleParams_ARnoise_splittau(trace,tau,event_times_init,params);
            end
        end
        disp(['Finished trace #' num2str(trace_ind)])
    end
    
    
        delete(this_pool)
        if params.cluster
            rmdir(jobdir)
   	end

else
    
            
    load_struct = load(params.init_method.template_file);
    template = load_struct.template;

    for trace_ind = 1:size(traces,1)
        disp(['Starting trace #' num2str(trace_ind)])
        trace = params.event_sign*traces(trace_ind,:);
        trace = trace - min(trace);

%         event_times_init_old = template_matching(-1*params.event_sign*traces(trace_ind,:), params.dt,...
%             params.init_method.tau, params.init_method.amp_thresh, params.init_method.conv_thresh);


        
        nfft = length(trace) + length(template);
        [~, event_times_init,event_sizes_init] = wiener_filter(trace,template,params.init_method.ar_noise_params,...
            nfft, params.dt, params.init_method.theshold, params.init_method.min_interval);
        event_times_init
        event_sizes_init
        results(trace_ind).event_times_init = event_times_init;
%         results(trace_ind).filtered_trace = filtered_trace;
        results(trace_ind).event_sizes_init = event_sizes_init;
%         assignin('base','event_times_init_old',event_times_init_old)
%         assignin('base','event_times_init',event_times_init)
        
        tau = [mean([params.tau1_min params.tau1_max]) mean([params.tau2_min params.tau2_max])]/params.dt;
        
        if isfield(params,'init_only') && ~params.init_only

            if params.direct_stim
    %             event_times_init = ceil(length(trace)*rand(1,length(trace)*params.p_spike));
                [results(trace_ind).trials, results(trace_ind).mcmc]  = sampleParams_ar_2taus_directstim(trace,tau,event_times_init,params);
            else
                [results(trace_ind).trials, results(trace_ind).mcmc]  = sampleParams_ARnoise_splittau(trace,tau,event_times_init,params);
            end
        end
    end
end

if isfield(params,'init_only') && ~params.init_only
    disp('finding min err...')
    % map sample
    for trace_ind = 1:size(traces,1)

        [results(trace_ind).map, results(trace_ind).map_ind] = max(results(trace_ind).trials.obj);

    end
end

% if params.is_grid
%     results_grid = unstack_results(results, rebuild_map);
% end

% results = results_grid;

disp('saving...')
save(params.full_save_string,'results','params','traces')

disp('done')


