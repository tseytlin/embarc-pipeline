function varargout = ACPC(varargin)
% ACPC MATLAB code for ACPC.fig
%      ACPC, by itself, creates a new ACPC or raises the existing
%      singleton*.
%
%      H = ACPC returns the handle to a new ACPC or the handle to
%      the existing singleton*.
%
%      ACPC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ACPC.M with the given input arguments.
%
%      ACPC('Property','Value',...) creates a new ACPC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ACPC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ACPC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ACPC

% Last Modified by GUIDE v2.5 18-Sep-2014 10:51:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ACPC_OpeningFcn, ...
                   'gui_OutputFcn',  @ACPC_OutputFcn, ...
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

% get subject name if possible
function subject = getSubject()
	subject = regexp(regexp(pwd(),'/embarc_.*fmriraw_\d+','match'),'[A-Z]+\d+_\dR\d','match');
	sz = size(subject);
	if sz(1) > 0
		subject = char(subject{1});
    else
        subject = regexp(regexp(pwd(),'/diamond_.*','match'),'\d+','match');
        sz = size(subject);
        if sz(1) > 0
            subject = char(subject{1}{2});
        else 
            subject = '';
        end
    end

% find image that goes with selected object
function image = getSelectedFile(lbl)
    % is it embarc study?
    image = '';
    subject = getSubject();
    if ~isempty(subject)
        % first check if there are fpm images for b0map
        image = dir(['dicom_' lbl '/fpm*.nii']);
        % check if we have vdm image
        if isempty(image)
            image = dir(['dicom_' lbl '/vdm*.nii']);
        end
        % if not, then lets view the nii image
        if isempty(image)
            image = dir(['dicom_' lbl '/' subject '_' lbl '*.nii']);
        end
        if ~isempty(image)
            image = [ 'dicom_' lbl '/' image(1).name ',1'];
        end
        if isempty(image)
            image = dir([lbl '/*.nii']);
            if ~isempty(image)
              image = [ lbl '/' image(1).name ',1'];
            end
        end
        if isempty(image)
            image = dir([lbl '/*.img']);
            if ~isempty(image)
              image = [ lbl '/' image(1).name ',1'];
            end
        end
    else
        % if not embarc, just pick a very first image that you find
        image = dir([lbl '/*.nii']);
        % if not, then lets view the img image
        if isempty(image)
            image = dir([lbl '/*.img']);
        end
        if ~isempty(image)
            image = [ lbl '/' image(1).name ',1'];
        end
    end
    
    

% load and update the list of images
function content = loadImages()
	list_content = {};
	% select directories recursively if not embarc
	subject = getSubject();
	if ~isempty(subject)
		sdir = dir('*');
		sdir = {sdir.name};
	else
		[st sdir] = dos('find -type d');
		[sdir st] = regexp(strtrim(sdir), '\n+', 'split'); 
	end
	
	for i=1: length(sdir)
		dname = char(sdir(i));
		% see if there is appropriate image file there
		if (~isempty(dir([dname '/*.nii'])) || ~isempty(dir([dname '/*.img'])))
			 % get the proper name for subdirectory
			 % if not embarc, just get directory name
			 name = regexp(dname,'dicom_(.+)','tokens');
			 if ~isempty(name)
			 	name = name{1};
			 else
			 	name = dname;
			 	% try to get last token of directory name
			 	%z = regexp(name,'/','split');
			 	%if length(z) > 0
			 	%	name = z(end);
			 	%end
			 end
			 
			 % skip special directories
			 filter = regexp(char(name),'b0map_bold_(6\.5.+|8\.5_[ri].*)','match');
			 if ~isempty(filter)
			 	continue;
			 end
			 			 
			 % add to list
			 if isempty(dir([dname '/*.acpc.csv']))
			 	list_content{end+1} = char(name);	
			 else
			 	list_content{end+1} = ['<html><font color=green>' char(name) '</font></html>'];
			 end
			 %list_content{end+1} = ['<html><font color=' color '>' char(name) '</font></html>'];	
			 %list_content{end+1} = char(name);	
		end
	end
	content = list_content;

function selection = getSelection(item)
	selection = item;
	m = regexp(item,'<html><font color=[a-z]+>(.*)</font></html>','tokens');
	if ~isempty(m)
		selection = char(m{1});
	end

% --- Executes just before ACPC is made visible.
function ACPC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ACPC (see VARARGIN)

% figure out subject name for EMBARC study
subject = getSubject();
if ~isempty(subject)
	set(handles.subject,'String',subject);
end

% now load the file list
set(handles.Images,'String',loadImages());

% initialize SPM
spm('defaults', 'FMRI');
spm('createintwin');
spm_jobman('initcfg');

% Choose default command line output for ACPC
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ACPC wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = ACPC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Images.
function Images_Callback(hObject, eventdata, handles)
% hObject    handle to Images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Images contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Images
  
% --- Executes during object creation, after setting all properties.
function Images_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in tregistration.
function tregistration_Callback(hObject, eventdata, handles)
% hObject    handle to tregistration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to registration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
files  = {};
subject = getSubject();
contents = cellstr(get(handles.Images,'String'));
lbl = getSelection(contents{get(handles.Images,'Value')});

% timeseries for fieldmap makes no sense
if length(findstr(lbl,'b0map')) > 0
	return
end


% add anatomical for embarc stuff
if ~isempty(subject) && exist(['dicom_anatomical/' subject '_anatomical.nii'])
	files{end+1} = ['dicom_anatomical/' subject '_anatomical.nii,1'];
end


% first lets figure out the N
if ~isempty(subject) && ~isempty(regexp(subject,'[A-Z]+\d+_\dR\d'))
	image = dir(['dicom_' lbl '/' subject '_' lbl '*.nii']);
	image = {image.name};
	% if single image, it maybe 4D volue
	if length(image) == 1	
		image = ['dicom_' lbl '/' char(image(1))];
		hdr = spm_vol(deblank(image));
		n = length(hdr);
		files{end+1} = sprintf('%s,%i',image,1);
		files{end+1} = sprintf('%s,%i',image,round(n/4));
		files{end+1} = sprintf('%s,%i',image,round(n/2));
		files{end+1} = sprintf('%s,%i',image,round(n*3/4));
		files{end+1} = sprintf('%s,%i',image,n);
	else
		% add first volume in 3D series
		n = length(image);
		files{end+1} = ['dicom_' lbl '/' char(image(1)) ',1'];
		files{end+1} = ['dicom_' lbl '/' char(image(round(n/4))) ',1'];
		files{end+1} = ['dicom_' lbl '/' char(image(round(n/2))) ',1'];
		files{end+1} = ['dicom_' lbl '/' char(image(round(n*3/4))) ',1'];
		files{end+1} = ['dicom_' lbl '/' char(image(n)) ',1'];
	end
else
	image = dir([lbl '/*.nii']);
	if isempty(image)
		image = dir([lbl '/*.img']);
	end
	image = {image.name};
	% if single image, it maybe 4D volue
	if length(image) == 1	
		image = [lbl '/' char(image(1))];
		hdr = spm_vol(deblank(image));
		n = length(hdr);
		files{end+1} = sprintf('%s,%i',image,1);
		files{end+1} = sprintf('%s,%i',image,round(n/4));
		files{end+1} = sprintf('%s,%i',image,round(n/2));
		files{end+1} = sprintf('%s,%i',image,round(n*3/4));
		files{end+1} = sprintf('%s,%i',image,n);
	else
		% add first volume in 3D series
		n = length(image);
		files{end+1} = [lbl '/' char(image(1)) ',1'];
		files{end+1} = [lbl '/' char(image(round(n/4))) ',1'];
		files{end+1} = [lbl '/' char(image(round(n/2))) ',1'];
		files{end+1} = [lbl '/' char(image(round(n*3/4))) ',1'];
		files{end+1} = [lbl '/' char(image(n)) ',1'];
	end
end

% now we need to get the t
if ~isempty(files)
	disp(files(1));

	clear matlabbatch;
	matlabbatch{1}.spm.util.checkreg.data = files;
	spm_jobman('serial', matlabbatch);
end


% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.Images,'String'));
lbl = getSelection(contents{get(handles.Images,'Value')});
image = getSelectedFile(lbl);
% now load the file list
set(handles.Images,'String',loadImages());

% do job
if ~isempty(image)
	clear matlabbatch;
	display(['Viewing: ' image]);
	matlabbatch{1}.spm.util.disp.data = {image};
	spm_jobman('serial', matlabbatch);
end


% --- Executes on button press in registration.
function registration_Callback(hObject, eventdata, handles)
% hObject    handle to registration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
files  = {};
subject = getSubject();
contents = cellstr(get(handles.Images,'String'));
% built image list
for i=1: length(contents)
	lbl = getSelection(contents{i});
	% is it embarc study?
	if ~isempty(subject)
		image = getSelectedFile(lbl);
        if ~isempty(image)
            files{end+1} = image; 
        end
        %image = dir(['dicom_' lbl '/' subject '_' lbl '*.nii']);
		%if length(findstr(lbl,'b0map')) > 0
		%	image = dir(['dicom_' lbl '/fpm*.nii']);
		%end
		
		%if ~isempty(image)
		%	files{end+1} = [ 'dicom_' lbl '/' char(image(1).name) ',1'];
		%end
	else
		% if not embarc, just pick a very first image that you find
		image = dir([lbl '/*.nii']);
		if isempty(image)
			image = dir([lbl '/*.img']);
		end
		if ~isempty(image)
			files{end+1} = [ lbl '/' char(image(1).name) ',1'];
		end
	end
	
end

if ~isempty(files)
	clear matlabbatch;
	matlabbatch{1}.spm.util.checkreg.data = files;
	spm_jobman('serial', matlabbatch);
end

% --- Executes on button press in display.
function comment_Callback(hObject, eventdata, handles)
    contents = cellstr(get(handles.Images,'String'));
    lbl = getSelection(contents{get(handles.Images,'Value')});
    image = getSelectedFile(lbl);
    [path name ext] = fileparts(image);
    ACPC_comment(path);
    
% --- Executes on button press in AC_Button.
function AC_Button_Callback(hObject, eventdata, handles)
% hObject    handle to AC_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st;
eval('ex=1;st.posinf;','ex=0;');
if ex
	set(handles.AC,'String',get(st.mp, 'String'));
end


function AC_Callback(hObject, eventdata, handles)
% hObject    handle to AC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AC as text
%        str2double(get(hObject,'String')) returns contents of AC as a double


% --- Executes during object creation, after setting all properties.
function AC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PC_Button.
function PC_Button_Callback(hObject, eventdata, handles)
% hObject    handle to PC_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st;
eval('ex=1;st.posinf;','ex=0;');
if ex
	set(handles.PC,'String',get(st.mp, 'String'));
end


function PC_Callback(hObject, eventdata, handles)
% hObject    handle to PC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PC as text
%        str2double(get(hObject,'String')) returns contents of PC as a double


% --- Executes during object creation, after setting all properties.
function PC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TP_Button.
function TP_Button_Callback(hObject, eventdata, handles)
% hObject    handle to TP_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st;
eval('ex=1;st.posinf;','ex=0;');
if ex
	set(handles.TP,'String',get(st.mp, 'String'));
end


function TP_Callback(hObject, eventdata, handles)
% hObject    handle to TP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TP as text
%        str2double(get(hObject,'String')) returns contents of TP as a double


% --- Executes during object creation, after setting all properties.
function TP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function right_Callback(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of right as text
%        str2double(get(hObject,'String')) returns contents of right as a double


% --- Executes during object creation, after setting all properties.
function right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function forward_Callback(hObject, eventdata, handles)
% hObject    handle to forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of forward as text
%        str2double(get(hObject,'String')) returns contents of forward as a double


% --- Executes during object creation, after setting all properties.
function forward_CreateFcn(hObject, eventdata, handles)
% hObject    handle to forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function up_Callback(hObject, eventdata, handles)
% hObject    handle to up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of up as text
%        str2double(get(hObject,'String')) returns contents of up as a double


% --- Executes during object creation, after setting all properties.
function up_CreateFcn(hObject, eventdata, handles)
% hObject    handle to up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pitch_Callback(hObject, eventdata, handles)
% hObject    handle to pitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pitch as text
%        str2double(get(hObject,'String')) returns contents of pitch as a double


% --- Executes during object creation, after setting all properties.
function pitch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function roll_Callback(hObject, eventdata, handles)
% hObject    handle to roll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roll as text
%        str2double(get(hObject,'String')) returns contents of roll as a double


% --- Executes during object creation, after setting all properties.
function roll_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yaw_Callback(hObject, eventdata, handles)
% hObject    handle to yaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yaw as text
%        str2double(get(hObject,'String')) returns contents of yaw as a double


% --- Executes during object creation, after setting all properties.
function yaw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function resizeX_Callback(hObject, eventdata, handles)
% hObject    handle to resizeX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resizeX as text
%        str2double(get(hObject,'String')) returns contents of resizeX as a double


% --- Executes during object creation, after setting all properties.
function resizeX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resizeX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function resizeY_Callback(hObject, eventdata, handles)
% hObject    handle to resizeY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resizeY as text
%        str2double(get(hObject,'String')) returns contents of resizeY as a double


% --- Executes during object creation, after setting all properties.
function resizeY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resizeY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function resizeZ_Callback(hObject, eventdata, handles)
% hObject    handle to resizeZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resizeZ as text
%        str2double(get(hObject,'String')) returns contents of resizeZ as a double


% --- Executes during object creation, after setting all properties.
function resizeZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resizeZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in get.
function get_Callback(hObject, eventdata, handles)
% hObject    handle to get (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st;
eval('ex=1;st.B;','ex=0;');
if ex
	set(handles.right,'String',st.B(1));
	set(handles.forward,'String',st.B(2));
	set(handles.up,'String',st.B(3));
	set(handles.pitch,'String',st.B(4));
	set(handles.roll,'String',st.B(5));
	set(handles.yaw,'String',st.B(6));
	set(handles.resizeX,'String',st.B(7));
	set(handles.resizeY,'String',st.B(8));
	set(handles.resizeZ,'String',st.B(9));
end

% --- Executes on button press in set.
function set_Callback(hObject, eventdata, handles)
% hObject    handle to set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st;
eval('ex=1;st.B;','ex=0;');
if ex 
	set(findobj('Callback','spm_image(''repos'',1)'),'String',get(handles.right,'String'));
	set(findobj('Callback','spm_image(''repos'',2)'),'String',get(handles.forward,'String'));
	set(findobj('Callback','spm_image(''repos'',3)'),'String',get(handles.up,'String'));
	set(findobj('Callback','spm_image(''repos'',4)'),'String',get(handles.pitch,'String'));
	set(findobj('Callback','spm_image(''repos'',5)'),'String',get(handles.roll,'String'));
	set(findobj('Callback','spm_image(''repos'',6)'),'String',get(handles.yaw,'String'));
	set(findobj('Callback','spm_image(''repos'',7)'),'String',get(handles.resizeX,'String'));
	set(findobj('Callback','spm_image(''repos'',8)'),'String',get(handles.resizeY,'String'));
	set(findobj('Callback','spm_image(''repos'',9)'),'String',get(handles.resizeZ,'String'));

	% update image
	for i=1:6
		obj = findobj('Callback',['spm_image(''repos'',' sprintf('%i',i) ')']);
		try, st.B(i) = eval(get(obj,'String')); end
    	set(obj,'String',st.B(i));
   		st.vols{1}.premul = spm_matrix(st.B);
    	spm_image('Zoom');
    	spm_image('Update');
	end
end


% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.right,'String','0');
	set(handles.forward,'String','0');
	set(handles.up,'String','0');
	set(handles.pitch,'String','0');
	set(handles.roll,'String','0');
	set(handles.yaw,'String','0');
	set(handles.resizeX,'String','1');
	set(handles.resizeY,'String','1');
	set(handles.resizeZ,'String','1');


% --- Executes on button press in Comment.
function Comment_Callback(hObject, eventdata, handles)
% hObject    handle to Comment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% comment reorientation valaues
ac = sscanf(get(handles.AC,'String'), '%g %g %g');
pc = sscanf(get(handles.PC,'String'), '%g %g %g');
tp = sscanf(get(handles.TP,'String'), '%g %g %g');

% do linear transformation for AC
if length(ac) == 3 
	% set fields
	set(handles.right,'String',sprintf('%g',-ac(1)));
	set(handles.forward,'String',sprintf('%g',-ac(2)));
	set(handles.up,'String',sprintf('%g',-ac(3)));
end

% do yaw and pitch
if length(ac) == 3  && length(pc) == 3 
	% compensate for origin
	pc(1) = pc(1) - ac(1);
	pc(2) = pc(2) - ac(2);
	pc(3) = pc(3) - ac(3);
	
	% do rotation
	%set(handles.roll,'String',sprintf('%g',-atan(pc(2)/pc(1))));
	set(handles.yaw,'String',sprintf('%g',-atan(pc(3)/pc(2))));
	set(handles.pitch,'String',sprintf('%g',-atan(pc(3)/pc(1))));
end

% do yaw and pitch
if length(ac) == 3  && length(tp) == 3 
	% compensate for origin
	tp(1) = tp(1) - ac(1);
	tp(2) = tp(2) - ac(2);
	tp(3) = tp(3) - ac(3);
	
	% do rotation
	set(handles.roll,'String',sprintf('%g',-atan(tp(2)/tp(1))));
end


% --- Executes when user attempts to close figure1.
function ACPC_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
quit();
