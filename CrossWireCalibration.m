function varargout = CrossWireCalibration(varargin)
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
handles.output = hObject; % Choose default command line output for CrossWireCalibration
guidata(hObject, handles); % Update handles structure
% addpath(genpath('E:\libs\matlab')); % add libraries
addpath(genpath('D:\nav\libs\matlab')); % add libraries


% --- Outputs from this function are returned to the command line.
function varargout = CrossWireCalibration_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;% Get default command line output from handles structure


% --- Executes on button press in btnPrev.
function btnPrev_Callback(hObject, eventdata, handles)
moveFrame('prev');
guiUpdate(handles);

function btnNext_Callback(hObject, eventdata, handles)
moveFrame('next');
guiUpdate(handles);

function btnLoadFile_Callback(hObject, eventdata, handles)                  
%fpath = 'E:\data\AmfiTrackCalibration150718';
%fpath = 'D:\data\AmfiTrackCalibration150718';
fpath = 'I:\WORKDIR\Anton\data\AmfiTrackCalibration18062018';
loadFile(fpath);
guiUpdate(handles);


% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
saveMarkers()

% --- Executes on button press in btnMark1.
function btnMark1_Callback(hObject, eventdata, handles)
putMarker(1, handles);
guiUpdate(handles);



% --- Executes on button press in btnUnmark.
function btnUnmark_Callback(hObject, eventdata, handles)
removeAllMarkers();
guiUpdate(handles);

% --- Executes on button press in btnMark2.
function btnMark2_Callback(hObject, eventdata, handles)
putMarker(2, handles);
guiUpdate(handles);

% --- Executes on button press in btnMark3.
function btnMark3_Callback(hObject, eventdata, handles)
putMarker(3, handles);
guiUpdate(handles);

% --- Executes on button press in btnCalib.
function btnCalib_Callback(hObject, eventdata, handles)
calibrate();


% --- Executes on button press in btnShowAllPlanes3D.
function btnShowAllPlanes3D_Callback(hObject, eventdata, handles)
showAllPlanes();

% --- Executes on button press in btnShowAllPoints3D.
function btnShowAllPoints3D_Callback(hObject, eventdata, handles)
showAllPoints();

%------------Application functions-----------------------------------------
function guiUpdate(handles)
    global iFrame nFrames I mark1 mark2 mark3 linSpace
    if nFrames == 0,    return; end;
    axes(handles.axesFrame);
    Imax = max(max(I(:,:,iFrame)));
    imagesc(linSpace.y, linSpace.x, I(:,:,iFrame), [0 Imax/20]);
    hold on
    plot(mark1(2,iFrame), mark1(1,iFrame), '*g');
    plot(mark2(2,iFrame), mark2(1,iFrame),  '*y');
    plot(mark3(2,iFrame), mark3(1,iFrame),  '*r');
    hold off
    xlabel('y (m)'); ylabel('x (m)'); 
    str = sprintf('frame %d of %d', iFrame, nFrames);
    set(handles.txtFrames,'string', str);

function moveFrame(strDirection)
    global iFrame nFrames
    if nFrames == 0,    return; end;
    if strcmp(strDirection, 'next')
        iFrame = iFrame+1;
        if iFrame>nFrames, iFrame = 1; end;
    end
    if strcmp(strDirection, 'prev')
        iFrame = iFrame-1;
        if iFrame == 0, iFrame = nFrames; end;
    end

function loadFile(filePath)
    global iFrame nFrames Hws I linSpace mark1 mark2 mark3 fpath Hst
    Hst = eye(4); %Hst(1,4) = 88.608e-3; Hst(3,4) = 18.386e-3;
    fpath = filePath;
    [ Hws, I, nFrames ] = loadData(fpath);
    if nFrames == 0, return; end;
    [ nSamples, nElements ] = size(I(:,:,1));
    [ gridx, gridy, ~] = linearArrayGridMsensor(0.1953e-3, nElements, 7.813e+6, 1540, nSamples); 
    linSpace.x = gridx(1:nSamples, 1)';
    linSpace.y = gridy(1,1:nElements);

    fullFileName = [ fpath, '\', 'MARK','.mat' ];
    if exist(fullFileName, 'file') == 2
        load(fullFileName);
        mark1 = MARK{1};
        mark2 = MARK{2};
        mark3 = MARK{3};
        clear('MARK');
    else
        mark1 = ones(4,nFrames)*NaN;
        mark2 = ones(4,nFrames)*NaN;
        mark3 = ones(4,nFrames)*NaN;
    end
    iFrame = 1;
    
function putMarker(mark, handles)
    global mark1 mark2 mark3 iFrame nFrames
    if mark<1 && mark>3,        return; end;
    if nFrames == 0,            return; end;
    [x y] = getpts(handles.axesFrame);
    p = [ y(1) x(1) 0 1]'; 
    switch mark
        case 1
            mark1(:,iFrame) = p;
        case 2
            mark2(:,iFrame) = p;
        case 3
            mark3(:,iFrame) = p;
    end
    
function removeAllMarkers()
    global mark1 mark2 mark3 iFrame nFrames
    if nFrames == 0,    return; end;
    mark1(:,iFrame) = ones(4,1)*NaN;
    mark2(:,iFrame) = ones(4,1)*NaN;
    mark3(:,iFrame) = ones(4,1)*NaN;
    
function saveMarkers
    global fpath mark1 mark2 mark3 nFrames
    if nFrames == 0, return; end;
    MARK = cell(3);
    MARK{1} = mark1;
    MARK{2} = mark2;
    MARK{3} = mark3;
    fullFileName = [ fpath, '\', 'MARK','.mat' ];
    save(fullFileName,'-v7.3', 'MARK');
    
function showAllPoints()
    global Hws mark1 mark2 mark3 nFrames Hst
    global fPoints
    if nFrames == 0, return; end;
    p1 = ones(4,nFrames);
    p2 = ones(4,nFrames);
    p3 = ones(4,nFrames);
    for i = 1:nFrames
        p1(:,i) = Hws(:,:,i)*Hst*mark1(:,i);
        p2(:,i) = Hws(:,:,i)*Hst*mark2(:,i);
        p3(:,i) = Hws(:,:,i)*Hst*mark3(:,i);
    end
    if isempty(fPoints), fPoints = figure(); end;
    if ~isvalid(fPoints), fPoints = figure(); end;
    figure(fPoints);
    subplot(2,2,1); plot3(p1(1,:), p1(2,:), p1(3,:), '*g', p2(1,:), p2(2,:), p2(3,:), '*y', p3(1,:), p3(2,:), p3(3,:), '*r'); xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)'); grid on; 
    subplot(2,2,2); plot(p1(1,:), p1(2,:), '*g', p2(1,:), p2(2,:), '*y', p3(1,:), p3(2,:), '*r'); xlabel('x(m)'); ylabel('y(m)'); grid on;
    subplot(2,2,3); plot(p1(1,:), p1(3,:), '*g', p2(1,:), p2(3,:), '*y', p3(1,:), p3(3,:), '*r'); xlabel('x(m)'); ylabel('z(m)'); grid on;
    subplot(2,2,4); plot(p1(2,:), p1(3,:), '*g', p2(2,:), p2(3,:), '*y', p3(2,:), p3(3,:), '*r'); xlabel('y(m)'); ylabel('z(m)'); grid on;
    
function showAllPlanes()
        global I Hws mark1 mark2 mark3 nFrames linSpace Hst
        global fPlanes
        if nFrames == 0, return; end;
        p1 = ones(4,nFrames);
        p2 = ones(4,nFrames);
        p3 = ones(4,nFrames);
        h = max(linSpace.x) - min(linSpace.x);
        w = max(linSpace.y) - min(linSpace.y);
        T = eye(4); T(1,4) = -w/2; R = getRotationMetrix(0,0,90);
        Hcorr = R*T;
        if isempty(fPlanes), fPlanes = figure(); end;
        if ~isvalid(fPlanes), fPlanes = figure(); end;
        figure(fPlanes);
        hold on
        for i = 1:nFrames
            p1 = Hws(:,:,i)*Hst*mark1(:,i);
            p2 = Hws(:,:,i)*Hst*mark2(:,i);
            p3 = Hws(:,:,i)*Hst*mark3(:,i);
            Iscl = I(:,:,i)/( max(max(I(:,:,i))) - min(min(I(:,:,i))) )*255;
            showImageInSpace(Iscl, h, w, Hws(:,:,i)*Hst*Hcorr);
            plot3(p1(1), p1(2),p1(3),'*g', p2(1),p2(2),p2(3),'*y', p3(1), p3(2),p3(3),'*r');
        end
        hold off
        grid on; xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)'); axis equal;

function res = f1(x)
        global Hws nFrames mark1
        %R = eye(4); R(1:3,1:3) = qGetR(x(1:4));
        %T = eye(4); T(1:3,4) = x(5:7);
        %Hst = R*T;
        Hst = eye(4); Hst(1:3,1:4) = reshape(x,3,4);
        nnorms = 0;
        for i = 1:nFrames
            if sum(isnan(mark1(:,i))) == 0
                nnorms = nnorms+1;
                v = Hws(:,:,i)*Hst*mark1(:,i);
                norms(nnorms) = norm(v(1:3));
            end
        end
        res = rms(norms);
        
function calibrate()
        global Hws Hst mark1 nFrames
        nx = 12;
        lambda = ones(nx,1)*0.02;
        x = eye(4); x = x(1:3,1:4); x = x(:);
        f = f1(x);
        xnew = x + 0.01;
        ni = 100;
        for i = 1:ni
            i
            dx = xnew - x; x = xnew;
            df = f - f1(x); f = f1(x);
            gradf = df./dx;
            xnew = x - lambda.*gradf;
            
            showAllPoints();
            Hst = eye(4); Hst(1:3,1:4) = reshape(x,3,4);
            f
            pause(0.1);
        end
        %R = eye(4); R(1:3,1:3) = qGetR(x(1:4));
        %T = eye(4); T(1:3,4) = x(5:7);
        %Hst = R*T;
        
