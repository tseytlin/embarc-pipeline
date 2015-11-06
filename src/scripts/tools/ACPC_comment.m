function varargout = ACPC_comment(varargin)
% ACPC_COMMENT MATLAB code for ACPC_comment.fig
%      ACPC_COMMENT, by itself, creates a new ACPC_COMMENT or raises the existing
%      singleton*.
%
%      H = ACPC_COMMENT returns the handle to a new ACPC_COMMENT or the handle to
%      the existing singleton*.
%
%      ACPC_COMMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ACPC_COMMENT.M with the given input arguments.
%
%      ACPC_COMMENT('Property','Value',...) creates a new ACPC_COMMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ACPC_comment_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ACPC_comment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ACPC_comment

% Last Modified by GUIDE v2.5 15-Dec-2014 19:46:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ACPC_comment_OpeningFcn, ...
                   'gui_OutputFcn',  @ACPC_comment_OutputFcn, ...
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

% --- Executes just before ACPC_comment is made visible.
function ACPC_comment_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ACPC_comment (see VARARGIN)

% Choose default command line output for ACPC_comment
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ACPC_comment wait for user response (see UIRESUME)
% uiwait(handles.ACPC_comment);
global directory;
if nargin
    directory = varargin{end};
else
    directory = pwd();
end

% check if there is an existing file selected
global file;
set(handles.directory_text,'String',directory);
if ~isempty(dir([directory '/PASS.txt']))
    file = [directory '/PASS.txt'];
    set(handles.radio_pass,'Value',1);
elseif ~isempty(dir([directory '/WARN.txt']))
    file = [directory '/WARN.txt'];
    set(handles.radio_warn,'Value',1); 
elseif ~isempty(dir([directory '/FAIL*.txt']))
    name = dir([directory '/FAIL*.txt']);
    name = name.name;
    file = [directory '/' name];
    set(handles.radio_fail,'Value',1);    
end
set(handles.directory_text,'String',getFile(handles));

% load content of text file int comment field
if exist(file,'file')
    text = fileread(file);
    set(handles.comment_text,'String',text); 
end

% -- get currently selected file
function file = getFile(handles)
    file = get(handles.directory_text,'String');
    if exist(file,'dir')
        file = [file '/' get(get(handles.comment_group,'SelectedObject'),'String') '.txt'];
   end    
    
    
% --- Outputs from this function are returned to the command line.
function varargout = ACPC_comment_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ok_button.
function ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% remove previous file
global file;
if exist(file,'file')
  delete(file);  
end

% write new comment file
text = get(handles.comment_text,'String');
fid = fopen(getFile(handles),'w');
[n m] = size(text);
for i=1:n
    fprintf(fid,'%s\n',text(i,:));
end
fclose(fid);

delete(gcf());

% --- Executes on button press in cancel_button.
function cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcf());


function comment_text_Callback(hObject, eventdata, handles)
% hObject    handle to comment_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of comment_text as text
%        str2double(get(hObject,'String')) returns contents of comment_text as a double


% --- Executes during object creation, after setting all properties.
function comment_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to comment_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in radio_fail.
function radio_fail_Callback(hObject, eventdata, handles)
% hObject    handle to radio_fail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_fail
global directory;
set(handles.directory_text,'String',[directory '/FAIL.txt']);

% --- Executes on button press in radio_warn.
function radio_warn_Callback(hObject, eventdata, handles)
% hObject    handle to radio_warn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_warn
global directory;
set(handles.directory_text,'String',[directory '/WARN.txt']);

% --- Executes on button press in radio_pass.
function radio_pass_Callback(hObject, eventdata, handles)
% hObject    handle to radio_pass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_pass
global directory;
set(handles.directory_text,'String',[directory '/PASS.txt']);
