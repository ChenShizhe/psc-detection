function varargout = detection_viewer(varargin)
% DETECTION_VIEWER MATLAB code for detection_viewer.fig
%      DETECTION_VIEWER, by itself, creates a new DETECTION_VIEWER or raises the existing
%      singleton*.
%
%      H = DETECTION_VIEWER returns the handle to a new DETECTION_VIEWER or the handle to
%      the existing singleton*.
%
%      DETECTION_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETECTION_VIEWER.M with the given input arguments.
%
%      DETECTION_VIEWER('Property','Value',...) creates a new DETECTION_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before detection_viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to detection_viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help detection_viewer

% Last Modified by GUIDE v2.5 18-Oct-2016 10:15:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @detection_viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @detection_viewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before detection_viewer is made visible.
function detection_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to detection_viewer (see VARARGIN)

% Choose default command line output for detection_viewer
handles.output = hObject;

handles.data.trace_grid = varargin{1};
handles.data.posterior_grid = varargin{2};
if length(varargin) == 3
    handles.data.events_grid = varargin{3};
else
    handles.data.events_grid = cellfun(@kmeans_estimate,handles.data.posterior_grid,handles.data.trace_grid,'UniformOutput',0);
    assignin('base','new_events_grid',handles.data.events_grid)
end

handles.current_row = 1;
handles.current_column = 1;
[handles.row_max, handles.col_max] = size(handles.data.trace_grid);

plot_trace_stack_grid(handles.data.trace_grid,Inf,1,0,[],handles.traces_grid_axes)

set(handles.grid_size,'String',['Grid Size: ' num2str(handles.row_max) ...
    ' x ' num2str(handles.col_max)])

% Update handles structure
guidata(hObject, handles);

update_results_axes(handles)


% UIWAIT makes detection_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = detection_viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function update_results_axes(handles)

axes(handles.detection_axes)

offset = 100; % expose
burn_in = 1;
colors_groups = hsv(100);
colors_groups = colors_groups(randperm(100),:);
colors_lines =  lines(100);

traces = handles.data.trace_grid{handles.current_row,handles.current_column};
posteriors = handles.data.posterior_grid{handles.current_row,handles.current_column};
events_by_trace = handles.data.events_grid{handles.current_row,handles.current_column};
num_traces = length(posteriors);




for i = 1:num_traces
    
    posterior = truncate_samples(posteriors(i),[burn_in length(posteriors(i).num_events)]);
    event_feature_means = events_by_trace{i};
    k = size(event_feature_means,1);
    
        
        if ~isfield(handles,'event_sign')
            recon_corr = corr([build_curve(event_feature_means,0,size(traces,2)/20000,1/20000,2000)' ...
                traces(i,:)' - traces(i,end)]);
            handles.event_sign = sign(recon_corr(2));
            guidata(handles.up,handles)
        end
        
        labels = ones(size(posterior.times));
%         length(posterior.times)
%         labels = 1:100:length(posterior.times);
%         size(labels)
%         labels = repmat(labels,105,1);
%         size(labels)
%         labels = labels(:);
%         size(labels)
%         labels = labels(1:length(posterior.times));
%         assignin('base','labels',labels)
        gscatter(posterior.times,-posterior.amp,labels,colors_lines(i,:),[],1,0)
        hold on
        gscatter(posterior.times,posterior.tau1,labels,colors_lines(i,:),[],1,0)
        hold on
        gscatter(posterior.times,posterior.tau2,labels,colors_lines(i,:),[],1,0)
        hold on
    
     
    if  k > 0 && length(posterior.amp) >= k
        scatter(event_feature_means(:,4), -event_feature_means(:,1),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        scatter(event_feature_means(:,4), event_feature_means(:,2),100,colors_lines(i,:),'x','LineWidth',2)
        hold on
        scatter(event_feature_means(:,4), event_feature_means(:,3),100,colors_lines(i,:),'x','LineWidth',2)
        hold on

        trace = -1.0*traces(i,:);
        trace = trace - min(trace);
        
        
        plot(1:size(traces,2),handles.event_sign*trace - 200 - offset*(i-1),'color',colors_lines(i,:))
        plot(1:size(traces,2), handles.event_sign*build_curve(event_feature_means,mean(posterior.base),size(traces,2)/20000,1/20000,2000) - 200 - offset*(i-1),'color',colors_lines(i,:),'Linewidth',2)
        xlim([1 size(traces,2)])

    else

        plot(1:size(traces,2),traces(i,:) - traces(i,end) - 200 - offset*(i-1),'color',colors_lines(i,:))
        xlim([1 size(traces,2)])
    end
end

hold off


% --- Executes on button press in down.
function down_Callback(hObject, eventdata, handles)
% hObject    handle to down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current_row_tmp = handles.current_row + 1;
if current_row_tmp > handles.row_max
    return
else
    handles.current_row = current_row_tmp;
    set(handles.goto_row,'String',num2str(handles.current_row));
    guidata(hObject, handles);
    update_results_axes(handles);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function goto_row_Callback(hObject, eventdata, handles)
% hObject    handle to goto_row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of goto_row as text
%        str2double(get(hObject,'String')) returns contents of goto_row as a double

current_row_tmp = str2double(get(hObject,'String'));
if current_row_tmp > handles.row_max || current_row_tmp < 1
    set(hObject,'String',num2str(handles.current_row));
    return
else
    handles.current_row = current_row_tmp;
    guidata(hObject, handles);
    update_results_axes(handles);
end



% --- Executes during object creation, after setting all properties.
function goto_row_CreateFcn(hObject, eventdata, handles)
% hObject    handle to goto_row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function goto_column_Callback(hObject, eventdata, handles)
% hObject    handle to goto_column (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of goto_column as text
%        str2double(get(hObject,'String')) returns contents of goto_column as a double

current_col_tmp = str2double(get(hObject,'String'));
if current_col_tmp > handles.col_max || current_col_tmp < 1
    set(hObject,'String',num2str(handles.data.current_column));
    return
else
    handles.current_column = current_col_tmp;
    guidata(hObject, handles);
    update_results_axes(handles);
end


% --- Executes during object creation, after setting all properties.
function goto_column_CreateFcn(hObject, eventdata, handles)
% hObject    handle to goto_column (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in up.
function up_Callback(hObject, eventdata, handles)
% hObject    handle to up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current_row_tmp = handles.current_row - 1;
if current_row_tmp < 1
    return
else
    handles.current_row = current_row_tmp;
    set(handles.goto_row,'String',num2str(handles.current_row));
    guidata(hObject, handles);
    update_results_axes(handles);
end


% --- Executes on button press in left.
function left_Callback(hObject, eventdata, handles)
% hObject    handle to left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current_col_tmp = handles.current_column - 1;
if current_col_tmp < 1
    return
else
    handles.current_column = current_col_tmp;
    set(handles.goto_column,'String',num2str(handles.current_column));
    guidata(hObject, handles);
    update_results_axes(handles);
end


% --- Executes on button press in right.
function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current_col_tmp = handles.current_column + 1;
if current_col_tmp > handles.col_max
    return
else
    handles.current_column = current_col_tmp;
    set(handles.goto_column,'String',num2str(handles.current_column));
    guidata(hObject, handles);
    update_results_axes(handles);
end
