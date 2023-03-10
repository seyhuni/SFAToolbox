function varargout = SFAToolbox(varargin)
% SFATOOLBOX MATLAB code for SFAToolbox.fig
%      SFATOOLBOX, by itself, creates a new SFATOOLBOX or raises the existing
%      singleton*.
%
%      H = SFATOOLBOX returns the handle to a new SFATOOLBOX or the handle to
%      the existing singleton*.
%
%      SFATOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SFATOOLBOX.M with the given input arguments.
%
%      SFATOOLBOX('Property','Value',...) creates a new SFATOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SFAToolbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SFAToolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SFAToolbox

% Last Modified by GUIDE v2.5 10-Mar-2023 14:45:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SFAToolbox_OpeningFcn, ...
                   'gui_OutputFcn',  @SFAToolbox_OutputFcn, ...
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


% --- Executes just before SFAToolbox is made visible.
function SFAToolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SFAToolbox (see VARARGIN)

% Choose default command line output for SFAToolbox
handles.output = hObject;
cla(handles.axes1,'reset');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SFAToolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SFAToolbox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EMGData; global SR; 

[SwallowFileName, path] = uigetfile('*.*','.mat');
cd(path)
File=fullfile(path,SwallowFileName);
SR=2000;
EMGData=load(File);
global EMGData; global SR; 

EMGData=EMGData.datas;
EMGData_=EMGData(500:end);
EMGData_=EMGData_-mean(EMGData_);
EMGData_=EMGData_/max(EMGData_);
EMGData_=EMGData_(EMGData_>-1);
SR=2000;
block_duration_datacorrected=EMGData_;
Band    = (2 / SR) * [25, 400];
[B, A]  = butter(6, Band, 'Bandpass');   
SRignal = filtfilt(B, A, double(block_duration_datacorrected));
    % figure
 Band    = (2 / SR) * [45, 55];
[B, A]  = butter(4, Band, 'stop');
block_duration_datacorrected_last = filtfilt(B, A, double(block_duration_datacorrected));
    N=length(block_duration_datacorrected_last);         %number of points
    t=(0:N-1)/SR;   %time vector
    sgf = sgolayfilt(block_duration_datacorrected_last,3,201);
    % figure
    % plot(sgf,'r')
    % title(['Savitzky-Golay Smoothed - block: ' num2str(blocks{block_num})])
    sgf=sgf(101:end-101);
    ynew=sgf/max(sgf); 
    % initialize filtered signal
    eogF = ynew;
%     figure
%     plot(eogF)
%     xlabel('Time')
%     ylabel('Amplitude')
%     xlim([1020.231193413892 260567.0467979074])
%     ylim([-0.0001904646550288437 0.0002626080637102741])
    % TKEO basic  % Teager–Kaiser energy operator to obtain EMG Bursts
    for i=2:length(eogF)-1
        eogF(i) = ynew(i)^2 - ynew(i-1)*ynew(i+1);
    end   
    % eogF=eogF(1:end-1);
    [c,l] = wavedec(eogF,7,'sym4'); % 8 level decomposition  

    % for t=1:8
    %     D(:,t)=wrcoef('d',c,l,'db6',t);
    % end

    % low frequency components to get swallow patterns, so that approximation
    % of wavelets are obtained
     clear A;
    for t=1:7
        A(:,t)=wrcoef('a',c,l,'sym4',t);
    end

    A8=A(:,7);  % A8 is the filtered and swallow pattern obtained signal

    % figure
    % plot(A8)
    A8_=A8(500:end-500);
%     figure
%     plot(A8_)
%     xlabel('Time')
%     ylabel('Amplitude')
            rmsSwallows = sqrt(movmean(A8_.^2, 2500));   % Burst Detection using RMS Value Over ‘WinLen’ Samples

    axes(handles.axes1); cla;
    plot(rmsSwallows/max(rmsSwallows))
xlim([0 length(rmsSwallows)])
xlabel('Time')
ylabel('Amplitude')

f = msgbox("The operation has been successfully completed","Success");

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EMGData; global SR; 
cla(handles.axes1,'reset');

EMGData_=EMGData(500:end);
EMGData_=EMGData_-mean(EMGData_);
EMGData_=EMGData_/max(EMGData_);
EMGData_=EMGData_(EMGData_>-1);
% EMGData_(163000:165000)=0;

% ylim([-1 1.1])
% [c,l]=wavedec(datas3,6,'sym4');
% a8 = wrcoef('a',c,l,'sym4',6);
% block_duration_datacorrected=datas3-a8;
SR=2000;
% faxis=linspace(-SR/2,SR/2,length(datas3));
% plot(faxis,fftshift(abs(fft(datas3))));
block_duration_datacorrected=EMGData_;
Band    = (2 / SR) * [25, 400];
[B, A]  = butter(6, Band, 'Bandpass');   
SRignal = filtfilt(B, A, double(block_duration_datacorrected));
    % figure
 Band    = (2 / SR) * [45, 55];
[B, A]  = butter(4, Band, 'stop');   
block_duration_datacorrected_last = filtfilt(B, A, double(block_duration_datacorrected));
    N=length(block_duration_datacorrected_last);         %number of points
    t=(0:N-1)/SR;   %time vector
    sgf = sgolayfilt(block_duration_datacorrected_last,3,201);
    % figure
    % plot(sgf,'r')
    % title(['Savitzky-Golay Smoothed - block: ' num2str(blocks{block_num})])
    sgf=sgf(101:end-101);
    ynew=sgf/max(sgf); 
    % initialize filtered signal
    eogF = ynew;
%     figure
%     plot(eogF)
%     xlabel('Time')
%     ylabel('Amplitude')
%     xlim([1020.231193413892 260567.0467979074])
%     ylim([-0.0001904646550288437 0.0002626080637102741])
    % TKEO basic  % Teager–Kaiser energy operator to obtain EMG Bursts
    for i=2:length(eogF)-1
        eogF(i) = ynew(i)^2 - ynew(i-1)*ynew(i+1);
    end   
    % eogF=eogF(1:end-1);
    [c,l] = wavedec(eogF,7,'sym4'); % 8 level decomposition  

    % for t=1:8
    %     D(:,t)=wrcoef('d',c,l,'db6',t);
    % end

    % low frequency components to get swallow patterns, so that approximation
    % of wavelets are obtained
     clear A;
    for t=1:7
        A(:,t)=wrcoef('a',c,l,'sym4',t);
    end

    A8=A(:,7);  % A8 is the filtered and swallow pattern obtained signal

    % figure
    % plot(A8)
    A8_=A8(500:end-500);
%     figure
%     plot(A8_)
%     xlabel('Time')
%     ylabel('Amplitude')
            rmsSwallows = sqrt(movmean(A8_.^2, 2500));   % Burst Detection using RMS Value Over ‘WinLen’ Samples
    axes(handles.axes1); cla;
        rmsSwallows=rmsSwallows/max(rmsSwallows);

    plot(rmsSwallows)
    button = 1;

    while sum(button) <=1   % read input with right click only for threshold
       [xx,yy,button] = ginput(1); % you can choose threshold by right clicking on plot 1 time
    end  % you can change
     thresholdValue=(yy);
     x=2000
        for k=x:(length(rmsSwallows)-x)
            interval_=rmsSwallows(k-x+1:k+x);
            if(interval_(ceil(length(interval_)/2))==max(interval_) & interval_(ceil(length(interval_)/2))>thresholdValue)
                swallowpeaks(k) = rmsSwallows(k);
            end
        end
    swallowpeaks_pos=find(swallowpeaks>thresholdValue);
    hv=ones(1,length(rmsSwallows));
    swallowoffset=sqrt(-1)*hv;
    swallowonset=swallowoffset;
    
    for k=1:length(swallowpeaks_pos)
        mt=rmsSwallows(swallowpeaks_pos(k):swallowpeaks_pos(k)+4500);
        mn=min(mt);
        swallowoffset(find(mt==mn)+swallowpeaks_pos(k)-1)=mn;
    end
    swallowoffset_last=find(swallowoffset~=sqrt(-1));

    for m=2:length(swallowpeaks_pos)
        mt2=rmsSwallows(swallowpeaks_pos(m)-5000:swallowpeaks_pos(m));
        mn2=min(mt2);
        swallowonset(swallowpeaks_pos(m)-5000-1+find(mt2==mn2))= mn2;
    end
    swallowonset_last=find(swallowonset~=sqrt(-1));
    
    hold on
    plot(swallowpeaks_pos,rmsSwallows(swallowpeaks_pos),'r+','MarkerFaceColor','r','LineWidth',1)
    hold on
    plot(swallowonset_last,rmsSwallows(swallowonset_last),'g*','MarkerFaceColor','g','LineWidth',1)
    hold on
    plot(swallowoffset_last,rmsSwallows(swallowoffset_last),'k+','MarkerFaceColor','k','LineWidth',1)
    xlabel('Time')
    ylabel('Amplitude')
    xlim([0 length(rmsSwallows)])
    ylim([0 1.4])
    swallow_amp=zeros(1,length(swallowpeaks_pos));
% format longG
for i=1:length(swallow_amp)
    swallow_amp(i)=rmsSwallows(swallowpeaks_pos(i))-rmsSwallows(swallowoffset_last(i));
end
format short g

BSAmp=0;
ASAmp=0;
BlackScreenSwallows=find(swallowpeaks_pos<60*SR);
AdvSwallows=find(swallowpeaks_pos>=60*SR);
for i=1:length(BlackScreenSwallows)
    BSAmp=BSAmp+swallow_amp(i);
end
for j=(length(BlackScreenSwallows)+1):length(BlackScreenSwallows)+length(AdvSwallows)
    ASAmp=ASAmp+swallow_amp(j);
end
MeanBlackScreen=BSAmp/length(BlackScreenSwallows);
MeanAdv=ASAmp/length(AdvSwallows);
    set(handles.edit1, 'String', length(AdvSwallows));
    set(handles.edit2, 'String', length(BlackScreenSwallows));
    set(handles.edit3, 'String', MeanAdv);
    set(handles.edit4, 'String', MeanBlackScreen);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
