function plot_times_posterior(results_file, posterior_file, trace_offset, varargin)

load(results_file,'params')
load(posterior_file,'posterior')
params.traces_filename
try
    load(params.traces_filename)
catch
    load('data/for-paper/all-evoked-ipscs.mat')
end

load(results_file,'traces')
% load('/home/shababo/Projects/Mapping/code/psc-detection/data/simulated-data-longer-traces-epsc.mat')
% load('/home/shababo/Desktop/simulated-data-longer-traces-epsc.mat')

if ~isfield(params,'start_ind')
    params.start_ind = 1;
end

if ~isfield(params,'duration')
    params.duration = size(traces,2);
end

traces = traces(:,params.start_ind:(params.start_ind + params.duration - 1));

if isfield(params,'traces_ind')
    traces = traces(params.traces_ind,:);
end

if length(varargin) > 2 && ~isempty(varargin{3})
    params.traces_ind = varargin{3};
else
    params.traces_ind = 1:size(traces,1);
end
traces = traces(params.traces_ind,:);

if ~isempty(varargin) && ~isempty(varargin{1})
    traces_ind = varargin{1};
    traces = traces(traces_ind,:);
    true_signal = true_signal(traces_ind,:);
end

if length(varargin) >= 2
    max_sample = varargin{2};
end

T = size(traces,2);

if isfield(params,'event_samples')
    event_samples = params.event_samples;
else
    event_samples = 4000;
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
    plot_trace_stack(traces,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 0]),'-',[.005 10],0)
    hold on
    if exist('true_signal','var')
        times_vec = zeros(size(posterior));
        for i = 1:size(posterior,1)        
            times_vec(i,ceil(true_event_times{i})) = max(max(posterior))+.1;
        end
        plot_scatter_stack(times_vec,trace_offset,[0 0],20,100,[0 0 0])
        hold on
        %plot_trace_stack(true_signal,trace_offset,bsxfun(@plus,zeros(length(traces),3),[0 0 1]),'-',[],80)
        %hold on
    end
    plot_scatter_stack(posterior,trace_offset,[0 0],5,100,[0 0 1])
    hold off

    title(strrep(results_file,'_','-'))

    % if length(varargin) > 1 && varargin{2}
    %     [dir,name,~] = fileparts(results_file);
    %     savefig([dir '/' name '.fig'])
    % end

    % map = params.event_sign*map_curves;




