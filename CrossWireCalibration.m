function varargout = CrossWireCalibration(varargin)
% CROSSWIRECALIBRATION MATLAB code for CrossWireCalibration.fig
%      CROSSWIRECALIBRATION, by itself, creates a new CROSSWIRECALIBRATION or raises the existing
%      singleton*.
%
%      H = CROSSWIRECALIBRATION returns the handle to a new CROSSWIRECALIBRATION or the handle to
%      the existing singleton*.
%
%      CROSSWIRECALIBRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CROSSWIRECALIBRATION.M with the given input arguments.
%
%      CROSSWIRECALIBRATION('Property','Value',...) creates a new CROSSWIRECALIBRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CrossWireCalibration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CrossWireCalibration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CrossWireCalibration

% Last Modified by GUIDE v2.5 02-May-2018 16:36:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CrossWireCalibration_OpeningFcn, ...
                   'gui_OutputFcn',  @CrossWireCalibration_OutputFcn, ...
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


% --- Executes just before CrossWireCalibration is made visible.
function CrossWireCalibration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CrossWireCalibration (see VARARGIN)

global fhData H0 mark1 mark2 mark3 iMark1 iMark2 iMark3 LinkFM
% Choose default command line output for CrossWireCalibration
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
H0 = eye(4);
%H0(2,4) = 27.98;
%H0(3,4) = 88.588;
%{
mark1.p = zeros(4,100);
mark1.frames = zeros(4,4,100);
iMark1 = 0;
mark2.p = zeors(4,100);
mark2.frames = zeros(4,4,100);
iMark2 = 0;
mark3.p = zeors(4,100);
mark3.frames = zeros(4,4,100);
iMark3 = 0;
LinkFM = [];
%}
% UIWAIT makes CrossWireCalibration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CrossWireCalibration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnPrev.
function btnPrev_Callback(hObject, eventdata, handles)
% hObject    handle to btnPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global iFrame nFrames fhData Imax
iFrame = iFrame-1;
if iFrame == 0, iFrame = nFrames; end;

strStatus = sprintf('%i of %i frames', iFrame, nFrames);
set(handles.txtFrames, 'string', strStatus);
axes(handles.axesFrame);
colormap(gray);
imagesc(fhData.img.linx, fhData.img.liny, fhData.img.frames(:,:,iFrame),[0,Imax/2]); axis equal;

% --- Executes on button press in btnNext.
function btnNext_Callback(hObject, eventdata, handles)
% hObject    handle to btnNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global iFrame nFrames fhData Imax LinkFM iMark1 mark1 iMark2 mark2 iMark3 mark3 
iFrame = iFrame+1;
if iFrame > nFrames, iFrame = 1; end;

strStatus = sprintf('%i of %i frames', iFrame, nFrames);
set(handles.txtFrames, 'string', strStatus);
axes(handles.axesFrame);
colormap(gray);
imagesc(fhData.img.linx, fhData.img.liny, fhData.img.frames(:,:,iFrame),[0,Imax/2]); axis equal;

hold on
if LinkFM(1,iFrame)>0
    p = mark1.pos(:,(LinkFM(1,iFrame)));
    plot(p(1), p(2), 'or');
end
if LinkFM(2,iFrame)>0
    p = mark2.pos(:,(LinkFM(2,iFrame)));
    plot(p(1), p(2), 'og');
end
if LinkFM(3,iFrame)>0
    p = mark3.pos(:,(LinkFM(3,iFrame)));
    plot(p(1), p(2), 'or');
end

% --- Executes on button press in btnLoadFile.
function btnLoadFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global iFrame nFrames Hws I grid
filePath = 'D:\data\AmfiTrackCalibration';
[ Hws, I, nFrames ] = loadData(filePath);
[ grid.x, grid.y, grid.z] = linearArrayGridMsensor(
iFrame = 1;



% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnMark1.
function btnMark1_Callback(hObject, eventdata, handles)
% hObject    handle to btnMark1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fhData iFrame mark1 iMark1 LinkFM
r = getrect(handles.axesFrame);
dx = fhData.img.linx(2) - fhData.img.linx(1);
dy = fhData.img.liny(2) - fhData.img.liny(1);
r(1) = floor( (r(1)-fhData.img.linx(1))/dx);
r(2) = floor( (r(2)-fhData.img.liny(1))/dy);
r(3) = floor( r(3)/dx);
r(4) = floor( r(4)/dy);
frame = fhData.img.frames(:,:,iFrame);
rio = frame( r(2):(r(2)+r(4)), r(1):(r(1)+r(3)) );
[ymax,j] = max(rio);
[xmax,i] = max(ymax);
xmax = (i+r(1))*dx + fhData.img.linx(1);
ymax = (j(i)+r(2))*dy + fhData.img.liny(1);
axes(handles.axesFrame);
hold on
plot(xmax,ymax,'or');
hold off

iMark1 = iMark1+1;
mark1.pos(:,iMark1) = [xmax, ymax, 0, 1];
mark1.frames(:,:,iMark1) = fhData.pos.frames(:,:,iFrame);
LinkFM(1,iFrame) = iMark1;



% --- Executes on button press in btnUnmark.
function btnUnmark_Callback(hObject, eventdata, handles)
% hObject    handle to btnUnmark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fhData iFrame mark1 iMark1 mark2 iMark2 mark3 iMark3 LinkFM

if LinkFM(1,iFrame)>0
    mark1( LinkFM(1,iFrame) ) = [];
    LinkFM(1,iFrame) = 0;
    iMark1  = iMark1-1;
end
if LinkFM(2,iFrame)>0
    mark2( LinkFM(2,iFrame) ) = [];
    LinkFM(2,iFrame) = 0;
    iMark2  = iMark2-1;
end
if LinkFM(3,iFrame)>0
    mark3( LinkFM(3,iFrame) ) = [];
    LinkFM(3,iFrame) = 0;
    iMark3  = iMark3-1;
end

% --- Executes on button press in btnMark2.
function btnMark2_Callback(hObject, eventdata, handles)
% hObject    handle to btnMark2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fhData iFrame mark2 iMark2 LinkFM
r = getrect(handles.axesFrame);
dx = fhData.img.linx(2) - fhData.img.linx(1);
dy = fhData.img.liny(2) - fhData.img.liny(1);
r(1) = floor( (r(1)-fhData.img.linx(1))/dx);
r(2) = floor( (r(2)-fhData.img.liny(1))/dy);
r(3) = floor( r(3)/dx);
r(4) = floor( r(4)/dy);
frame = fhData.img.frames(:,:,iFrame);
rio = frame( r(2):(r(2)+r(4)), r(1):(r(1)+r(3)) );
[ymax,j] = max(rio);
[xmax,i] = max(ymax);
xmax = (i+r(1))*dx + fhData.img.linx(1);
ymax = (j(i)+r(2))*dy + fhData.img.liny(1);
axes(handles.axesFrame);
hold on
plot(xmax,ymax,'og');
hold off

iMark2 = iMark2+1;
mark2.pos(:,iMark2) = [xmax, ymax, 0, 1];
mark2.frames(:,:,iMark2) = fhData.pos.frames(:,:,iFrame);
LinkFM(2,iFrame) = iMark2;

% --- Executes on button press in btnMark3.
function btnMark3_Callback(hObject, eventdata, handles)
% hObject    handle to btnMark3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fhData iFrame  mark3 iMark3 LinkFM
r = getrect(handles.axesFrame);
dx = fhData.img.linx(2) - fhData.img.linx(1);
dy = fhData.img.liny(2) - fhData.img.liny(1);
r(1) = floor( (r(1)-fhData.img.linx(1))/dx);
r(2) = floor( (r(2)-fhData.img.liny(1))/dy);
r(3) = floor( r(3)/dx);
r(4) = floor( r(4)/dy);
frame = fhData.img.frames(:,:,iFrame);
rio = frame( r(2):(r(2)+r(4)), r(1):(r(1)+r(3)) );
[ymax,j] = max(rio);
[xmax,i] = max(ymax);
xmax = (i+r(1))*dx + fhData.img.linx(1);
ymax = (j(i)+r(2))*dy + fhData.img.liny(1);
axes(handles.axesFrame);
hold on
plot(xmax,ymax,'ob');
hold off

mark3.pos(:,iMark3) = [xmax, ymax, 0, 1];
mark3.frames(:,:,iMark3) = fhData.pos.frames(:,:,iFrame);
iMark3 = iMark3 + 1;
LinkFM(3,iFrame) = iMark3;


% --- Executes on button press in btnCalib.
function btnCalib_Callback(hObject, eventdata, handles)
% hObject    handle to btnCalib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnOpen3D.
function btnOpen3D_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpen3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
