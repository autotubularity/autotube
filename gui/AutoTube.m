function varargout = AutoTube(varargin)
% AUTOTUBE MATLAB code for AutoTube.fig
%      AUTOTUBE, by itself, creates a new AUTOTUBE or raises the existing
%      singleton*.
%
%      H = AUTOTUBE returns the handle to a new AUTOTUBE or the handle to
%      the existing singleton*.
%
%      AUTOTUBE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTOTUBE.M with the given input arguments.
%
%      AUTOTUBE('Property','Value',...) creates a new AUTOTUBE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AutoTube_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AutoTube_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AutoTube

% Last Modified by GUIDE v2.5 03-Feb-2018 17:46:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AutoTube_OpeningFcn, ...
                   'gui_OutputFcn',  @AutoTube_OutputFcn, ...
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


function setup

    addpath(genpath('../code/'));
    addpath(genpath('../libs/'));

% --- Executes just before AutoTube is made visible.
function AutoTube_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AutoTube (see VARARGIN)

setup;

% Choose default command line output for AutoTube
handles.output = hObject;

set(handles.imageInput, 'visible','off');
set(handles.text_inputImage, 'visible','off');

set(handles.imageDetection, 'visible','off');
set(handles.text_imageDetection, 'visible','off');

set(handles.imageOutput, 'visible','off');
set(handles.text_imageOutput, 'visible','off');

% Filling variables
handles.fileFolder = pwd;

handles.data   = [];
handles.segmentation = [];
handles.denoising = [];
handles.adjustment = [];
handles.mask = [];
handles.strParams = [];
handles.currPos = get(gcf, 'position');

handles = defaultParams(handles);
handles = loadModels(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AutoTube wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles = defaultParams(handles)

    data.fileName.input     = '_input';
    data.fileName.adjusted  = '_adjusted';
    data.fileName.illumination = '_illumination';
    data.fileName.denoised  = '_denoised';
    data.fileName.tubularity = '_tubularity';
    data.fileName.edges     = '_edges';
    data.fileName.preproc   = '_preproc';
	data.fileName.mask_bw   = '_mask_bw';
    data.fileName.clean     = '_clean';
    data.fileName.hull     = '_hull';
    data.fileName.skel      = '_skel';
    data.fileName.branch    = '_branches';
    data.fileName.rings     = '_rings';
    
    data.fileName.extension = '.tif';

%     
%     data.microscope.objectiveMagnification = 10;
%     set(handles.textEdit_ObjectiveMagnification, 'string', sprintf('%d', data.microscope.objectiveMagnification));
%     
%     data.microscope.lenseMagnification = 1;
%     set(handles.textEdit_LenseMagnification, 'string', sprintf('%d', data.microscope.lenseMagnification));
%     
%     data.microscope.CMount = 1;
%     set(handles.textEdit_CMount, 'string', sprintf('%d', data.microscope.CMount));

    data.microscope.cameraPixSize = 1.0;
    set(handles.textEdit_CameraPixelSize, 'string', sprintf('%.2f', data.microscope.cameraPixSize));
    
%     data.microscope.binning = 1;
%     set(handles.textEdit_Binning, 'string', sprintf('%d', data.microscope.binning));
    

    %# Setting default values for drop-down menues
    set(handles.popupmenuColorChannel, 'Value',2 ); %# Green Color Channel
    set(handles.popupmenuAdjustment, 'Value',1 ); %# Auto-Contrast
    set(handles.popupmenuDenoise, 'Value',1 ); %# BM3D
    set(handles.popupmenuThreshold, 'Value',1 ); %# MultiOtsu
    
    set(handles.checkbox_Adjustment, 'value', 1);
    
    adjustmentStr = get(handles.popupmenuAdjustment, 'String');
    adjustmentVal = get(handles.popupmenuAdjustment,'Value');
    handles.adjustment.type = lower(adjustmentStr{adjustmentVal});
    
    set(handles.checkbox_Illumination, 'value', 1);
    data.illumination.kernelSize = 51;
    set(handles.textEdit_IlluminationDiskSze, 'string', sprintf('%d', data.illumination.kernelSize));
    
    set(handles.checkbox_Denoising, 'value', 1);
    denoiseStr = get(handles.popupmenuDenoise, 'String');
    denoiseVal = get(handles.popupmenuDenoise,'Value');
    handles.denoising.type = lower(denoiseStr{denoiseVal});
    
    segmentStr = get(handles.popupmenuThreshold, 'String');
    segmentVal = get(handles.popupmenuThreshold,'Value');
    handles.segmentation.type = lower(segmentStr{segmentVal});
    
    set(handles.checkbox_SmallRegions, 'value', 1);
    data.minAreaPercent = 0.01;
    %data.minAreaPercent = 0.1;
    set(handles.textEdit_RegionAnalysis,'string',data.minAreaPercent);
    
    data.brushRad = 15;
    set(handles.slider_brushSize, 'min', 1);
    set(handles.slider_brushSize, 'max', 35);
    set(handles.slider_brushSize, 'sliderstep', [1/34 1/34]);
    set(handles.slider_brushSize, 'value', data.brushRad);
    set(handles.text_brushSize, 'string', sprintf('Brush size: %d', data.brushRad) );
    
    set(handles.checkbox_Spurs, 'value', 1);
    data.spurLength = 15;
    set(handles.slider_spurLength, 'min', 0);
    %set(handles.slider_spurLength, 'max', 500);
    set(handles.slider_spurLength, 'max', 100);
    set(handles.slider_spurLength, 'sliderstep', [1/100 1/100]);
    set(handles.slider_spurLength, 'value', data.spurLength);
    set(handles.text_spurLength, 'string', sprintf('Spur Length: %d', data.spurLength) );
    
    set(handles.checkbox_SpatialBranches, 'value', 1);
    data.spatialBranches = 10;
    set(handles.slider_spatialBranches, 'min', 0);
    set(handles.slider_spatialBranches, 'max', 100);
    set(handles.slider_spatialBranches, 'sliderstep', [5/100 5/100]);
    %set(handles.slider_spatialBranches, 'sliderstep', [1/100 1/100]);
    set(handles.slider_spatialBranches, 'value', data.spatialBranches);
    set(handles.text_spatialBranches, 'string', sprintf('Spatial Distance: %d', data.spatialBranches) );
    
    data.paintTimeOut = 0.15;
    
    handles.data = data;
    
    handles.mask.status = 'iddle';
    handles.mask.maskStackLimit = 30;
    handles.mask.maskColor = 0.8;
    handles.mask.maskColor_Overlay = 0.4;
    handles.mask.alphaColor = 'gray';
    handles.mask.tStart = tic; %#changed 12/02 by javier

    return

    
function handles = loadModels(handles)

    % loading edge method.
    opts = edgesTrain();                % default options (good settings)
    opts.modelDir = 'models/';          % model will be in models/forest
    opts.modelFnm = 'modelBsds';        % model name
    opts.nPos = 5e5; opts.nNeg=5e5;     % decrease to speedup training
    opts.useParfor = 0;                 % parallelize if sufficient memory
    
    handles.edge = [];
    currentdir   = pwd;
	tic, handles.edge.model=edgesTrain(opts); toc; % will load model if already trained
	cd(currentdir);
   
    handles.edge.model.opts.multiscale = 0;          % for top accuracy set multiscale=1
    handles.edge.model.opts.sharpen = 2;             % for top speed set sharpen=0
    handles.edge.model.opts.nTreesEval = 4;          % for top speed set nTreesEval=1
    handles.edge.model.opts.nThreads = 4;            % max number threads for evaluation
    handles.edge.model.opts.nms = 0;                 % set to true to enable nms


% --- Outputs from this function are returned to the command line.
function varargout = AutoTube_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
    varargout{1} = handles.output;

function trimStr = removeBlankSpaceStr(myStr)

    myStr = strtrim(myStr);
    trimStr = myStr(~isspace(myStr));
    return;
    
function handles = getHandlesDirectories(handles)

    handles.data.microscope.cameraPixSize = str2double(get(handles.textEdit_CameraPixelSize, 'String'));

    adjusmentStr  = get(handles.popupmenuAdjustment, 'String');
    adjustmentVal  = get(handles.popupmenuAdjustment, 'Value');
    adjustmentType = lower(adjusmentStr{adjustmentVal});
    adjustmentType = removeBlankSpaceStr(adjustmentType);

    doTubularity = get(handles.checkbox_Tubularity,'value');
    valTubularity = str2double(get(handles.textEdit_Tubularity, 'string'));
    tubuStr = 'noTubu';
    if doTubularity
        tubuStr = sprintf('wiTubu%02d', valTubularity);
    end
    
    segmentStr  = get(handles.popupmenuThreshold, 'String');
    segmentVal  = get(handles.popupmenuThreshold, 'Value');
    segmentType = lower(segmentStr{segmentVal});
    
    dateStr = datestr(datetime('today'));
    %dateStr = datestr(now,'dd.mm.yyyy-HH.MM.SS');
    timeNow = '';
    
    handles.strParams.tubuStr = tubuStr;
    handles.strParams.segmentType = segmentType;

    outDir      = fullfile(handles.fileFolder, 'output', sprintf('%s_%s_%s_%s%s', adjustmentType, tubuStr, segmentType, dateStr, timeNow));
    handles.outDir  = outDir;
    
    outInputDir = fullfile(outDir, '00_input');
    handles.outInputDir  = outInputDir;
    
    outAdjDir      = fullfile(outDir, '01_adjusted');
    handles.outAdjDir  = outAdjDir;
    
    outIlluDir     = fullfile(outDir, '02_illumination');
    handles.outIlluDir  = outIlluDir;
    
    outDenDir      = fullfile(outDir, '03_denoised');
    handles.outDenDir  = outDenDir;
    
    outTubularityDir   = fullfile(outDir, '04_tubularity');
    handles.outTubularityDir  = outTubularityDir;
    
    outEdgesDir    = fullfile(outDir, '04a_edges');
    handles.outEdgesDir  = outEdgesDir;
    
    outPreprocDir    = fullfile(outDir, '04b_preproc');
    handles.outPreprocDir  = outPreprocDir;
    
    outMaskBWDir   = fullfile(outDir, '05_thresholded');
    handles.outMaskBWDir  = outMaskBWDir;

    outCleanDir    = fullfile(outDir, '06a_cleaned');
    handles.outCleanDir  = outCleanDir;
    
    outHullDir    = fullfile(outDir, '06b_hull');
    handles.outHullDir  = outHullDir;
    
    outSkelDir     = fullfile(outDir, '07a_skel');
    handles.outSkelDir  = outSkelDir;
    
    outBranchesDir = fullfile(outDir, '07b_branches');
    handles.outBranchesDir  = outBranchesDir;
    
    outRingsDir = fullfile(outDir, '07c_rings');
    handles.outRingsDir  = outRingsDir;
    
    outOverlayDir = fullfile(outDir, '08_overlays');
    handles.outOverlayDir  = outOverlayDir;
    
    outStatsDir = fullfile(outDir, '09_statistics');
    handles.outStatsDir  = outStatsDir;
    

function handles = createDirectories(handles)

    handles = getHandlesDirectories(handles);
    
    if ~exist(handles.outDir, 'dir'), mkdir(handles.outDir); end

    if ~exist(handles.outInputDir, 'dir'), mkdir(handles.outInputDir); end
  
    if ~exist(handles.outAdjDir, 'dir'), mkdir(handles.outAdjDir); end
    
    if ~exist(handles.outIlluDir, 'dir'), mkdir(handles.outIlluDir); end

    if ~exist(handles.outDenDir, 'dir'), mkdir(handles.outDenDir); end
    
    if ~exist(handles.outTubularityDir, 'dir'), mkdir(handles.outTubularityDir); end
    
    if ~exist(handles.outEdgesDir, 'dir'), mkdir(handles.outEdgesDir); end
    
    if ~exist(handles.outPreprocDir, 'dir'), mkdir(handles.outPreprocDir); end
    
    if ~exist(handles.outMaskBWDir, 'dir'), mkdir(handles.outMaskBWDir); end

    if ~exist(handles.outCleanDir, 'dir'), mkdir(handles.outCleanDir); end
    
    if ~exist(handles.outHullDir, 'dir'), mkdir(handles.outHullDir); end
    
    if ~exist(handles.outSkelDir, 'dir'), mkdir(handles.outSkelDir); end

    if ~exist(handles.outBranchesDir, 'dir'), mkdir(handles.outBranchesDir); end
    
    if ~exist(handles.outRingsDir, 'dir'), mkdir(handles.outRingsDir); end
    
    if ~exist(handles.outOverlayDir, 'dir'), mkdir(handles.outOverlayDir); end
    
    if ~exist(handles.outStatsDir, 'dir'), mkdir(handles.outStatsDir); end

%     adjusmentStr  = get(handles.popupmenuAdjustment, 'String');
%     adjustmentVal  = get(handles.popupmenuAdjustment, 'Value');
%     adjustmentType = lower(adjusmentStr{adjustmentVal});
%     adjustmentType = removeBlankSpaceStr(adjustmentType);
% 
%     doTubularity = get(handles.checkbox_Tubularity,'value');
%     valTubularity = str2double(get(handles.textEdit_Tubularity, 'string'));
%     tubuStr = 'noTubu';
%     if doTubularity
%         tubuStr = sprintf('wiTubu%02d', valTubularity);
%     end
%     
%     segmentStr  = get(handles.popupmenuThreshold, 'String');
%     segmentVal  = get(handles.popupmenuThreshold, 'Value');
%     segmentType = lower(segmentStr{segmentVal});
%     
%     %dateStr = datestr(datetime('today'));
%     dateStr = datestr(now,'dd.mm.yyyy-HH.MM.SS');
%     timeNow = '';
%     
%     handles.strParams.tubuStr = tubuStr;
%     handles.strParams.segmentType = segmentType;
% 
%     outDir      = fullfile(handles.fileFolder, 'output', sprintf('%s_%s_%s_%s%s', adjustmentType, tubuStr, segmentType, dateStr, timeNow));
%     if ~exist(outDir, 'dir'), mkdir(outDir); end
%     handles.outDir  = outDir;
%     
%     outInputDir = fullfile(outDir, '00_input');
%     if ~exist(outInputDir, 'dir'), mkdir(outInputDir); end
%     handles.outInputDir  = outInputDir;
%     
%     outAdjDir      = fullfile(outDir, '01_adjusted');
%     if ~exist(outAdjDir, 'dir'), mkdir(outAdjDir); end
%     handles.outAdjDir  = outAdjDir;
%     
%     outIlluDir     = fullfile(outDir, '02_illumination');
%     if ~exist(outIlluDir, 'dir'), mkdir(outIlluDir); end
%     handles.outIlluDir  = outIlluDir;
%     
%     outDenDir      = fullfile(outDir, '03_denoised');
%     if ~exist(outDenDir, 'dir'), mkdir(outDenDir); end
%     handles.outDenDir  = outDenDir;
%     
%     outTubularityDir   = fullfile(outDir, '04_tubularity');
%     if ~exist(outTubularityDir, 'dir'), mkdir(outTubularityDir); end
%     handles.outTubularityDir  = outTubularityDir;
%     
%     outEdgesDir    = fullfile(outDir, '04a_edges');
%     if ~exist(outEdgesDir, 'dir'), mkdir(outEdgesDir); end
%     handles.outEdgesDir  = outEdgesDir;
%     
%     outPreprocDir    = fullfile(outDir, '04b_preproc');
%     if ~exist(outPreprocDir, 'dir'), mkdir(outPreprocDir); end
%     handles.outPreprocDir  = outPreprocDir;
%     
%     outMaskBWDir   = fullfile(outDir, '05_thresholded');
%     if ~exist(outMaskBWDir, 'dir'), mkdir(outMaskBWDir); end
%     handles.outMaskBWDir  = outMaskBWDir;
% 
%     outCleanDir    = fullfile(outDir, '06a_cleaned');
%     if ~exist(outCleanDir, 'dir'), mkdir(outCleanDir); end
%     handles.outCleanDir  = outCleanDir;
%     
%     outHullDir    = fullfile(outDir, '06b_hull');
%     if ~exist(outHullDir, 'dir'), mkdir(outHullDir); end
%     handles.outHullDir  = outHullDir;
%     
%     outSkelDir     = fullfile(outDir, '07a_skel');
%     if ~exist(outSkelDir, 'dir'), mkdir(outSkelDir); end
%     handles.outSkelDir  = outSkelDir;
%     
%     outBranchesDir = fullfile(outDir, '07b_branches');
%     if ~exist(outBranchesDir, 'dir'), mkdir(outBranchesDir); end
%     handles.outBranchesDir  = outBranchesDir;
%     
%     outRingsDir = fullfile(outDir, '07c_rings');
%     if ~exist(outRingsDir, 'dir'), mkdir(outRingsDir); end
%     handles.outRingsDir  = outRingsDir;
%     
%     outOverlayDir = fullfile(outDir, '08_overlays');
%     if ~exist(outOverlayDir, 'dir'), mkdir(outOverlayDir); end
%     handles.outOverlayDir  = outOverlayDir;
%     
%     outStatsDir = fullfile(outDir, '09_statistics');
%     if ~exist(outStatsDir, 'dir'), mkdir(outStatsDir); end
%     handles.outStatsDir  = outStatsDir;
    
    
function handles = createSubDirectories(handles, imageName)
% createSubDirectories: this function creates subdirectories for each of
% the images generated after each of the processing steps.

    outDenDir = fullfile(handles.outDenDir);
    if ~exist(outDenDir, 'dir'), mkdir(outDenDir); end
    
    outSkelDir = fullfile(handles.outSkelDir, imageName);
    if ~exist(outSkelDir, 'dir'), mkdir(outSkelDir); end
    
    outBranchesDir = fullfile(handles.outBranchesDir, imageName);
    if ~exist(outBranchesDir, 'dir'), mkdir(outBranchesDir); end
     
    outOverlayDir = fullfile(handles.outOverlayDir, imageName);
    if ~exist(outOverlayDir, 'dir'), mkdir(outOverlayDir); end

    
function handles = getFileList(handles)
    
    opts.regexp   = '*.*';
    files     = dir(fullfile(handles.fileFolder, opts.regexp));
    dirIdx    = [files.isdir];
    fileNames = {files(~dirIdx).name}';
    validIndex = ~ismember(fileNames,{'.','..'}); 
    fileNames(~validIndex) = [];
    validFileNames = {};
    
    noFiles = length(fileNames);
   
    for ii=1:noFiles
        fileName = fileNames{ii};
        [~, ~, extension] = fileparts(fileName);
        extension = upper(extension);
        switch lower(extension)
            case {'.png', '.bmp', '.jpg', '.tif', '.tiff'}
                % Allow only PNG, TIF, JPG, or BMP images
                validFileNames = [validFileNames fileName];
            otherwise
        end
    end
    
    set(handles.listImageFiles, 'string', validFileNames);

    return

% --- Executes on button press in button_LoadIms.
function handles = button_LoadIms_Callback(hObject, eventdata, handles)
% hObject    handle to button_LoadIms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    returnValue = uigetdir(handles.fileFolder, 'Select folder');
    
    if returnValue ~= 0
        handles.fileFolder = returnValue;
		handles = getFileList(handles);
%         handles = createDirectories(handles);  
        handles = getHandlesDirectories(handles);
		guidata(hObject, handles);
    end
    
    return
    
    
function opts = getHandlesValues(handles) 
% getHandlesValues: this function retrieves the values of the selected
% options pressed by the user.

    opts.adjustment   = handles.adjustment.type;
    
    opts.illumination = [];
    opts.illumination.seSzeIllu = str2double(get(handles.textEdit_IlluminationDiskSze, 'string'));
	opts.params.doIlluCorrection = get(handles.checkbox_Illumination, 'value');
    
    opts.denoising = [];
    opts.params.doDenoising = get(handles.checkbox_Denoising, 'value');
    opts.denoising.type = handles.denoising.type;

    opts.tubularity = [];
    opts.tubularity.sigmaEnd = str2double(get(handles.textEdit_Tubularity, 'string'));
    opts.params.doTubularity = get(handles.checkbox_Tubularity,'value');
    
    opts.params.doEdges      = get(handles.checkbox_Edges,'value');

    opts.cleaning = [];
    opts.params.doCleanIsolatedRegions  = get(handles.checkbox_SmallRegions, 'value');
    opts.cleaning.minAreaPercent        = str2double(get(handles.textEdit_RegionAnalysis, 'string'))/100;
    if ~opts.params.doCleanIsolatedRegions
        opts.cleaning.minAreaPercent = 0;
    end
 
    opts.params.doFillHoles     = get(handles.checkbox_Holes, 'value');
    opts.cleaning.minHoleSize   = str2double(get(handles.textEdit_FillHoles, 'string'));
    if ~opts.params.doFillHoles
        opts.cleaning.minHoleSize = 0;
    end
    
    opts.params.doTubularity = get(handles.checkbox_Tubularity,'value');
    
    opts.params.doEdges      = get(handles.checkbox_Edges,'value');
    
    opts.params.doHull = get(handles.checkbox_convHull, 'value');

    opts.skeleton.cleanSpur         = get(handles.checkbox_Spurs,'value');
    opts.skeleton.spurLength        = handles.data.spurLength;
    
    opts.skeleton.cleanBranches     = get(handles.checkbox_SpatialBranches,'value');
    opts.skeleton.spatialBranches   = handles.data.spatialBranches;
    
    opts.thresh   = handles.segmentation.type;
    opts.minDist  = 5000000;
    
    return
    
function [fileNameSkelStrClean, fileNameSkelStrCleanOverlay] = getSkelFilename(opts)

    if opts.skeleton.cleanSpur == 1 &&  opts.skeleton.spurLength ~= 0
        fileNameSkelStrClean        = sprintf('_skelcut%03d', opts.skeleton.spurLength);
        fileNameSkelStrCleanOverlay = sprintf('_skelcut%03d', opts.skeleton.spurLength);
    elseif opts.skeleton.cleanSpur == 1 &&  opts.skeleton.spurLength == 0
        fileNameSkelStrClean        = sprintf('_skelcutall');
        fileNameSkelStrCleanOverlay = sprintf('_skelcutall');
    else
        fileNameSkelStrClean        = sprintf('_skelcutnone');
        fileNameSkelStrCleanOverlay = sprintf('_skelcutnone');
    end
    
    return
    
    
function [fileNameBranchStrClean, fileNameBranchStrCleanOverlay] = getBranchFilename(opts)
	if opts.skeleton.cleanBranches == 1 &&  opts.skeleton.spatialBranches ~= 0
        fileNameBranchStrClean         = sprintf('_branchcut%03d', opts.skeleton.spatialBranches);
        fileNameBranchStrCleanOverlay  = sprintf('_branchcut%03d', opts.skeleton.spatialBranches);
	else
        fileNameBranchStrClean         = sprintf('_branchcutnone');
        fileNameBranchStrCleanOverlay  = sprintf('_branchcutnone');
	end
    return
    
    
function handles = maskStackSzeControl(hObject, handles)

    if(handles.mask.currIdx > handles.mask.maskStackLimit)
        excedingIdx = handles.mask.currIdx - handles.mask.maskStackLimit;
        handles.mask.maskStack(1:excedingIdx) = [];
        handles.mask.currIdx = handles.mask.currIdx - excedingIdx;
    end
    
    guidata(hObject, handles);
    return 
  
    
function handles = addMask(hObject, handles)
% addMask: given the stack of edited masks, this function adds a new mask
% into the stack.
    
    handles.mask.currIdx = numel(handles.mask.maskStack)+1;
    handles = maskStackSzeControl(hObject, handles);
    handles.mask.maskStack{handles.mask.currIdx} = handles.mask.currMask;
    
    guidata(hObject, handles);
    return 


function handles = swapMask(hObject, handles)
    % This function replaces the current mask on the top of the stack
    % (does not increment the mask stack index). This is used in the
    % painting mode to avoid rapidly filling the stack with each paint
    % stroke. Instead, a watching timer is implemented so that the
    % stack index is only incremented after a certain amount of time
    % has past if painting.
    %handles.mask.maskStack{handles.mask.currIdx} = handles.mask.maskStack{max(handles.mask.currIdx-1,1)} | handles.mask.currMask;
    handles.mask.maskStack{handles.mask.currIdx} = handles.mask.currMask;


function pixLoc = fixPixLocation(handles, pixLoc)
% This function ensures that the location of the brushed region is integer and within the
% bounds of the image.
    
    % crop and swap x/y
    pixLoc = round(pixLoc);
    
    pixLoc(pixLoc < 1) = 1;
    
    [im_rows, im_cols, ~] = size(handles.images.im_clean);
    
    if(pixLoc(1) > im_rows)
      pixLoc(1) = im_rows;
    end
    
    if(pixLoc(2) > im_cols)
      pixLoc(2) = im_cols;
    end

    
function brushMask = drawBrush(handles, pixLoc)
% drawBrush: creates a binary image containg a circle at position 'pixLoc'
% and of radius: 'handles.data.brushRad'

    brushRad = handles.data.brushRad;
    brushSze = [brushRad, brushRad];
    
    pixLoc = fixPixLocation(handles, pixLoc);
    
    [im_rows, im_cols, ~] = size(handles.images.im_clean);
    
    yy = repmat([1:im_rows]',[1 im_cols]);
    xx = repmat([1:im_cols],[im_rows 1]);
    
    brushMask = ...
      (yy - pixLoc(1)).^2./brushSze(1).^2 + ...
      (xx - pixLoc(2)).^2./brushSze(2).^2 ...
      < 1;
    
    
function drawOverlayMask(hObject, eventdata, handles)
    brushMask = single(handles.mask.maskStack{handles.mask.currIdx});
    brushMaskOverlay = doOverlay(handles.images.im_in_color, brushMask, [0.9,0.1], handles.mask.alphaColor);
    handleImageName = 'imageDetection';
    
    set(eval(['handles.' handleImageName]), 'visible','on');
    imshow(brushMaskOverlay, 'InitialMagnification', 'fit', 'parent', eval(['handles.' handleImageName]));
    drawnow;
    

function maskColor = mask2Color(mask, rgb)
    [irows, icols, ~] = size(mask);
    maskColor = repmat( mask, [1 1 3] );
    maskColor = reshape(maskColor, irows * icols, 3) .* rgb;
    maskColor = reshape(maskColor, irows, icols, 3);


function drawPreviewMask(hObject, eventdata, handles)
    brushMask = handles.mask.previewMask;
    
    rgb_pix = handles.mask.maskColor_Overlay;
    if length(rgb_pix) == 1, rgb_pix = repmat(rgb_pix,1,3); end
    brushMaskColor = mask2Color(brushMask, rgb_pix);
    
    handleImageName = 'imageDetection';
    set(eval(['handles.' handleImageName]), 'visible','on');

    if ~isfield(handles.mask, 'currIdx')
        handles.mask.currIdx = 0;
    end

    if handles.mask.currIdx >= 1
        imMask = handles.mask.maskStack{handles.mask.currIdx};
    else
        if isfield(handles.images, 'im_clean')
            imMask = handles.images.im_clean;
        else
            imMask = zeros(size(handles.images.im_clean));
        end
    end
    
    imOverlay = doOverlay(handles.images.im_in_color, imMask, [0.9,0.1], handles.mask.alphaColor);
    imMaskOverlay = doOverlay(imOverlay, brushMask, [0.9,0.1], handles.mask.alphaColor);

    imshow(imMaskOverlay, 'InitialMagnification', 'fit', 'parent', eval(['handles.' handleImageName]));
    
    drawnow;

    
function handles = previewRegion(hObject, eventdata, handles, pixLoc)
    brushMask = drawBrush(handles, pixLoc);
    handles.mask.previewMask = zeros(size(handles.images.im_clean));
    handles.mask.previewMask(brushMask) = true;
    
    drawPreviewMask(hObject, eventdata, handles);
    
    
function handles = paintRegion(hObject, eventdata, handles, pixLoc)

    brushMask = drawBrush(handles, pixLoc);
    handles.mask.currMask = handles.mask.maskStack{handles.mask.currIdx};

    maskFillType = get(hObject,'SelectionType');

    if strcmpi(maskFillType, 'alt')
        handles.mask.currMask(brushMask) = false;
    else
        handles.mask.currMask(brushMask) = true;
    end

    handles = addMask(hObject, handles);

    drawOverlayMask(hObject, eventdata, handles);
    guidata(hObject, handles);
%     set(gcf,'UserData',{handles});
    return

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function mask_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to mainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	pixLoc = get(handles.imageDetection,'CurrentPoint');
    pixLoc = pixLoc(1, [2 1]); %# accessing mouse location and switching/flipping values.
    
    handles.mask.status = 'painting';
    
    if isfield(handles, 'images')
        handles = paintRegion(hObject, eventdata, handles, pixLoc);
    end
    
    guidata(hObject, handles);
    
    return
    
% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background. 
function startmovit(hObject, eventdata, handles)
% Unpack gui object

    pixLoc = get(handles.imageDetection,'CurrentPoint');
    pixLoc = pixLoc(1, [2 1]);

    if isfield(handles, 'images')
        handles.mask.status = 'painting';
        
        handles = paintRegion(hObject, eventdata, handles, pixLoc);

        % Set callbacks
        set(hObject,'WindowButtonMotionFcn',@movit);
        set(hObject,'WindowButtonUpFcn',@stopmovit);
        
        % Store gui object
        guidata(hObject, handles);
    end

    return
    
function flag=isMultipleCall()
    s = dbstack();
    % s(1) corresponds to isMultipleCall
    if numel(s)<=2, flag=false; return; end
    % compare all functions on stack to name of caller
    count = sum(strcmp(s(2).name,{s(:).name}));
    % is caller re-entrant?
    if count>1, flag=true; else flag=false; end
    
    
% --- Executes on mouse motion over figure - except title and menu.
function mask_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to mainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    pixLoc = get(handles.imageDetection, 'CurrentPoint');
    pixLoc = pixLoc(1, [2 1]); %# accessing mouse location and switching/flipping values.
    
    if isfield(handles, 'images')
        if(pixLoc(1) <= 0 || pixLoc(1)>size(handles.images.im_clean,1) || pixLoc(2) <= 0 || pixLoc(2)>size(handles.images.im_clean,2))
            return;
        end
        if strcmpi(handles.mask.status, 'painting')
            fprintf('# %s\n', handles.mask.status);
            handles = paintRegion(hObject, eventdata, handles, pixLoc);
            handles = previewRegion(hObject, eventdata, handles, pixLoc);
        elseif strcmpi(handles.mask.status, 'iddle') || strmcpi(handles.mask.status, 'ending')
            fprintf('# %s\n', handles.mask.status);
            handles = previewRegion(hObject, eventdata, handles, pixLoc);
        else
            error('# Unknown option: %s!\n', handles.mask.status);
        end
    else
        return;
    end
    guidata(hObject, handles);
    
    return
    
function movit_preview(hObject, eventdata, handles)
    pixLoc = get(handles.imageDetection, 'CurrentPoint');
    pixLoc = pixLoc(1, [2 1]);

    if isfield(handles, 'images')
        if(pixLoc(1) <= 0 || pixLoc(1)>size(handles.images.im_clean,1) || pixLoc(2) <= 0 || pixLoc(2)>size(handles.images.im_clean,2))
            return;
        end
        handles = previewRegion(hObject, eventdata, handles, pixLoc);
        guidata(hObject, handles);
    end
    return 
    
% --- Executes on mouse motion over figure - except title and menu.
% function movit(hObject, eventdata, handles)
function movit(hObject, eventdata)
% hObject    handle to mainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if isMultipleCall();  return;  end

    handles = guidata( ancestor(hObject, 'figure') ); %# last valid
    
    pixLoc = get(handles.imageDetection, 'CurrentPoint');
    pixLoc = pixLoc(1, [2 1]);

    if isfield(handles, 'images')
        if(pixLoc(1) <= 0 || pixLoc(1)>size(handles.images.im_clean,1) || pixLoc(2) <= 0 || pixLoc(2)>size(handles.images.im_clean,2))
            return;
        end
        if strcmpi(handles.mask.status, 'painting')
            fprintf('# %s\n', handles.mask.status);
            handles = paintRegion(hObject, eventdata, handles, pixLoc); %#last valid  
        elseif strcmpi(handles.mask.status, 'iddle')
            fprintf('# %s\n', handles.mask.status);
            handles = previewRegion(hObject, eventdata, handles, pixLoc);
        else
            error('# Unknown option: %s!\n', handles.mask.status);
        end
        guidata(hObject, handles);
    end
    return
    
    
% --- Executes on mouse release over figure.
function mask_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to mainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.mask.status = 'iddle';
    guidata(hObject, handles);
    return;
   
% --- Executes on mouse release over figure.
% function stopmovit(hObject, eventdata, handles)
function stopmovit(hObject, eventdata)
% hObject    handle to mainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Clean up the evidence ...
    set(hObject,'WindowButtonMotionFcn','');
    set(hObject,'WindowButtonUpFcn','');

    handles = guidata( ancestor(hObject, 'figure') );
    handles.mask.status = 'iddle';
    guidata(hObject, handles);

    return;

    
function handles = analyzeImage(hObject, imageFilename, handles)
% analyzeImage: main function for processing a given images.

    opts = getHandlesValues(handles);

    imChannelIdx = get(handles.popupmenuColorChannel, 'Value');
    
    branchOverlayColor = 'y';
    
    if opts.skeleton.spurLength == 0 && opts.skeleton.cleanSpur == 1
        fileNameGUISkelStr      = sprintf('%scutall', handles.data.fileName.skel);
    else
        fileNameGUISkelStr      = sprintf('%scut%03d', handles.data.fileName.skel, opts.skeleton.spurLength);
    end
    
    [imPath, imFileName, ~]  = fileparts(imageFilename);
    imFileName = strrep(imFileName, ' ', '_');
    
    handles = createSubDirectories(handles, imFileName);
    
    [fileNameSkelStrClean, fileNameSkelStrCleanOverlay] = getSkelFilename(opts);
    fileNameSkelStrOverlay      = sprintf('_skelcutnone');
    
    [fileNameBranchStrClean, fileNameBranchStrCleanOverlay] = getBranchFilename(opts);
    fileNameBranchStrOverlay       = sprintf('_branchcutnone');
    
    %# Output filenames
    fileNameInputImRGB = fullfile(handles.outInputDir, sprintf('%s%s_rgb%s', imFileName, handles.data.fileName.input, handles.data.fileName.extension));
    
    fileNameAdjusted = fullfile(handles.outAdjDir, sprintf('%s%s%s', imFileName, handles.data.fileName.adjusted, handles.data.fileName.extension));
    fileNameAdjustedRGB = fullfile(handles.outAdjDir, sprintf('%s%s_rgb%s', imFileName, handles.data.fileName.adjusted, handles.data.fileName.extension));
    
    fileNameIllumination    = fullfile(handles.outIlluDir, sprintf('%s%s%s', imFileName, handles.data.fileName.illumination, handles.data.fileName.extension));
    fileNameIlluminationRGB = fullfile(handles.outIlluDir, sprintf('%s%s_rgb%s', imFileName, handles.data.fileName.illumination, handles.data.fileName.extension));
    
    fileNameDenoised        = fullfile(handles.outDenDir, sprintf('%s%s%s', imFileName, handles.data.fileName.denoised, handles.data.fileName.extension));
    fileNameDenoisedRGB     = fullfile(handles.outDenDir, sprintf('%s%s_rgb%s', imFileName, handles.data.fileName.denoised, handles.data.fileName.extension));
    
    fileNameTubularity      = fullfile(handles.outTubularityDir, sprintf('%s%s%s', imFileName, handles.data.fileName.tubularity, handles.data.fileName.extension));
    fileNameTubularityRGB   = fullfile(handles.outTubularityDir, sprintf('%s%s_rgb%s', imFileName, handles.data.fileName.tubularity, handles.data.fileName.extension));
    
    fileNameEdges           = fullfile(handles.outEdgesDir, sprintf('%s%s.mat', imFileName, handles.data.fileName.edges));
    fileNameEdgesRGB        = fullfile(handles.outEdgesDir, sprintf('%s%s_rgb%s', imFileName, handles.data.fileName.edges, handles.data.fileName.extension));
    
    fileNamePreproc         = fullfile(handles.outPreprocDir, sprintf('%s%s%s', imFileName, handles.data.fileName.preproc, handles.data.fileName.extension));
    fileNamePreprocRGB      = fullfile(handles.outPreprocDir, sprintf('%s%s_rgb%s', imFileName, handles.data.fileName.preproc, handles.data.fileName.extension));
    
    fileNameMask_bw         = fullfile(handles.outMaskBWDir, sprintf('%s%s%s', imFileName, handles.data.fileName.mask_bw, handles.data.fileName.extension));
    
    fileNameHull            = fullfile(handles.outHullDir, sprintf('%s%s%s', imFileName, handles.data.fileName.hull, handles.data.fileName.extension));
    
    fileNameClean           = fullfile(handles.outCleanDir, sprintf('%s%s%s', imFileName, handles.data.fileName.clean, handles.data.fileName.extension));
    
    
    fileNameSkel            = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s%s', imFileName, sprintf('%scutnone', handles.data.fileName.skel), sprintf('%scutnone', handles.data.fileName.branch), handles.data.fileName.extension));
    fileNameBranch          = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s%s', imFileName, sprintf('%scutnone', handles.data.fileName.skel), sprintf('%scutnone', handles.data.fileName.branch), handles.data.fileName.extension));
    
    fileNameSkelClean           = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrClean, fileNameBranchStrClean, '_clean', handles.data.fileName.extension));
    fileNameSkelOverlay         = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrOverlay, fileNameBranchStrOverlay, '_overlay', handles.data.fileName.extension));
    fileNameSkelCleanOverlay    = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrCleanOverlay, fileNameBranchStrCleanOverlay, '_clean_overlay', handles.data.fileName.extension));
	fileNameSkelDelOverlayInputImRGB = fullfile(handles.outOverlayDir, imFileName, sprintf('%s%s%s%s_rgb_skeldel_input%s', imFileName, handles.data.fileName.skel, fileNameGUISkelStr, ...
                                                 sprintf('%scut%03d',handles.data.fileName.branch, opts.skeleton.spatialBranches), handles.data.fileName.extension));
    fileNameSkelOrigOverlayInputImRGB     = fullfile(handles.outOverlayDir, imFileName, sprintf('%s%s_rgb_skeldel_input%s', imFileName, sprintf('%s_orig', handles.data.fileName.skel), handles.data.fileName.extension));
    
    fileNameBranchClean         = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrClean, fileNameBranchStrClean, '_clean', handles.data.fileName.extension));
    fileNameBranchOverlay       = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrOverlay, fileNameBranchStrOverlay, '_overlay', handles.data.fileName.extension));
    fileNameBranchCleanOverlay  = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrCleanOverlay, fileNameBranchStrCleanOverlay, '_clean_overlay', handles.data.fileName.extension));
    fileNameGUIBranchesDel  = fullfile(handles.outOverlayDir, imFileName, sprintf('%s%s%s%s%s', imFileName, handles.data.fileName.branch, fileNameGUISkelStr, sprintf('%scut%03d',handles.data.fileName.branch, opts.skeleton.spatialBranches), handles.data.fileName.extension));
    fileNameGUIBranchesOrig = fullfile(handles.outOverlayDir, imFileName, sprintf('%s%s%s', imFileName, sprintf('%s_orig', handles.data.fileName.branch), handles.data.fileName.extension));
    
    fileNameSkelBranchOverlayOrig     = fullfile(handles.outOverlayDir, sprintf('%s_orig%s', imFileName, handles.data.fileName.extension));
    fileNameSkelBranchOverlayProc     = fullfile(handles.outOverlayDir, sprintf('%s_proc%s%s%s', imFileName, fileNameSkelStrClean, fileNameBranchStrClean, handles.data.fileName.extension));
    
    %hProgressBar = waitbar(0, 'Initializing waitbar...', 'Name', 'Progress Bar', 'WindowStyle', 'modal');
    hProgressBar = waitbar(0, 'Initializing waitbar...', 'Name', 'Progress Bar');

    noSteps = 4;
    
    if opts.params.doIlluCorrection == 1
        noSteps = noSteps + 1;
    end
    
    if opts.params.doDenoising
        noSteps = noSteps + 1;
    end
    
    if opts.params.doEdges
        noSteps = noSteps + 1;
    end
    
    if opts.params.doHull
        noSteps = noSteps + 1;
    end
    
    if opts.params.doTubularity
        noSteps = noSteps + 1;
    end
    
    progressFactor = 1/noSteps;
    progressMultiplier = 0;
    
    reProcess = false;
    
%     try
        handles.mask.maskStack = {};
        handles.mask.currIdx = 0; 
        
        handles.images.im_in = (im2uint8(imread(imageFilename)));
        
        [imrows, imcols, imch] = size(handles.images.im_in);
        
        imInputColor = handles.images.im_in;
        if imch == 1
            imInputColor = zeros(imrows, imcols, 3);
            imInputColor(:,:,imChannelIdx) = im2double(handles.images.im_in);
            imInputColor = im2uint8(imInputColor);
        end
        imwrite(imInputColor, fileNameInputImRGB);
        handles.images.im_in_color = imInputColor;
        
        %drawImageResultsSingle(handles.images, handles, 'im_in_color', 'text_inputImage', 'imageInput');
        drawImageResultsSingle(hObject, handles.images, handles, 'im_in_color', 'text_inputImage', 'imageInput');
        
        if imch > 1
            imHist = zeros(1, imch);
            for ch=1:imch
                imHist(1,ch) = sum(sum(handles.images.im_in(:,:,ch)~=0));
            end
            [~, imChannelIdx] = max(imHist);
            handles.images.im_in = handles.images.im_in(:,:,imChannelIdx);
        end
        
        if imChannelIdx == 2
            skelOverlayColor = 'r';
        else
            skelOverlayColor = 'g';
        end
        
        % Show Progress Bar
        progressMultiplier = progressMultiplier + 1;
        waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Adjusting Image Contrast ... %2.0f%% along',progressMultiplier*progressFactor*100));
        
        %# Performing Intensity-correction
        if exist(fileNameAdjusted, 'file') == 2
            handles.images.im_adj = imread(fileNameAdjusted);
        else
            handles.images.im_adj    = doAdjust(handles.images.im_in, opts);
            imwrite(handles.images.im_adj, fileNameAdjusted);
            imAdjustedColor = zeros(imrows, imcols, 3);
            imAdjustedColor(:,:,imChannelIdx) = im2double(handles.images.im_adj);
            imAdjustedColor = im2uint8(imAdjustedColor);
            imwrite(imAdjustedColor, fileNameAdjustedRGB);
            reProcess = true;
        end
        
        %# Performing Illumination-correction
        handles.images.im_illu   = handles.images.im_adj;
        if opts.params.doIlluCorrection == 1
            progressMultiplier = progressMultiplier + 1;
            waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Normalizing Illumination ... %2.0f%% along',progressMultiplier*progressFactor*100));
        
            if exist(fileNameIllumination, 'file') == 2 && ~reProcess
                handles.images.im_illu = imread(fileNameIllumination);
            else
                handles.images.im_illu   = doIlluCorrection(handles.images.im_adj, opts);
                imwrite(handles.images.im_illu, fileNameIllumination);
                imIluColor = zeros(imrows, imcols, 3);
                imIluColor(:,:,imChannelIdx) = im2double(handles.images.im_illu);
                imIluColor = im2uint8(imIluColor);
                imwrite(imIluColor, fileNameIlluminationRGB);
                reProcess = true;
            end
        end
        
        %# Performing Image-Denoising
        handles.images.im_den = handles.images.im_illu;
        if opts.params.doDenoising == 1
            
            progressMultiplier = progressMultiplier + 1;
            waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Denoising Image ... %2.0f%% along',progressMultiplier*progressFactor*100));
            
            if exist(fileNameDenoised, 'file') == 2  && ~reProcess
                handles.images.im_den = imread(fileNameDenoised);
                imDen       = doMaxminNorm(handles.images.im_den, 0, 1, 0, 255);
                imDenColor  = imread(fileNameDenoisedRGB);
            else
                handles.images.im_den = doDenoising(handles.images.im_illu, opts);
                imwrite(handles.images.im_den, fileNameDenoised);
                imDenColor = zeros(imrows, imcols, 3);
                imDenColor(:,:,imChannelIdx) = handles.images.im_den;
                imwrite((imDenColor), fileNameDenoisedRGB);
                reProcess = true;
                imDen  = doMaxminNorm(handles.images.im_den, 0, 1, 0, 1);
            end
        end
        
        %# Performing Tubularity
        handles.images.im_tubularity = handles.images.im_den;
        if opts.params.doTubularity == 1
            progressMultiplier = progressMultiplier + 1;
            waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Tubularity Extraction ... %2.0f%% along',progressMultiplier*progressFactor*100));
            
            if exist(fileNameTubularity, 'file') == 2  && ~reProcess
                handles.images.im_tubularity = imread(fileNameTubularity);
                imTubularityColor  = imread(fileNameTubularityRGB);
            else
                handles.images.im_tubularity = doFrangi(double(im2uint8(handles.images.im_den)), opts);
                imwrite(handles.images.im_tubularity, fileNameTubularity);
                imTubularityColor = zeros(imrows, imcols, 3);
                imTubularityColor(:,:,imChannelIdx) = handles.images.im_tubularity;
                imwrite((imTubularityColor), fileNameTubularityRGB);
                reProcess = true;
            end
        end
        
        %# Performing finer-tube detection
        handles.images.preProc  = handles.images.im_tubularity;
        if opts.params.doEdges == 1
            
            progressMultiplier = progressMultiplier + 1;
            waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Extracting Boundaries for finer tubes ... %2.0f%% along',progressMultiplier*progressFactor*100));
            
            if exist(fileNameEdges, 'file') == 2 && ~reProcess
                im_edges = load(fileNameEdges);
                im_edges = im_edges.im_edges;
                handles.images.im_edges = im_edges;
            else
                im_edges_t = doEdges(handles.images.im_den, handles.edge);
                im_edges = im_edges_t;
                handles.images.im_edges = im_edges;
                save(fileNameEdges, 'im_edges', '-v7.3');
                
                imEdgesColor = zeros(imrows, imcols, 3);
                imEdgesColor(:,:,imChannelIdx) = im2double(handles.images.im_edges);
                imEdgesColor = im2uint8(imEdgesColor);
                imwrite(imEdgesColor, fileNameEdgesRGB);
                
                reProcess = true;
            end
            
            imEdge = doMaxminNorm(handles.images.im_edges, 0, 1, 0, 1);
            imOut  = 0.5*(imEdge+imDen);
            sigma = 0.5; %0.25;
            hsze  = 2*ceil(3*sigma)+1;
            gf = fspecial('gaussian', hsze, sigma);
            handles.images.preProc   = imfilter(imOut, gf, 'same', 'replicate');
        end
        
        imPreprocColor = zeros(imrows, imcols, 3);
        imPreprocColor(:,:,imChannelIdx) = im2double(handles.images.preProc);
        imPreprocColor = im2uint8(imPreprocColor);
        imwrite(im2uint8(handles.images.preProc), fileNamePreproc);
        imwrite(imPreprocColor, fileNamePreprocRGB);
        
        % Show Progress Bar
        progressMultiplier = progressMultiplier + 1;
        waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Image Thresholding ... %2.0f%% along',progressMultiplier*progressFactor*100));
        
        %# Applying Threshold
        if exist(fileNameMask_bw, 'file') == 2 && ~reProcess
            handles.images.im_bw = imread(fileNameMask_bw);
        else
            [handles.images.im_bw, handles.images.im_mask_gray]  = doThresh(handles.images.preProc, opts);
            imwrite(handles.images.im_bw, fileNameMask_bw);
            
            reProcess  = true;
        end
        
        %# Cleaning Binary Mask
        if exist(fileNameClean, 'file') == 2 && ~reProcess
            handles.images.im_clean = imread(fileNameClean);
        else
            handles.images.im_clean = doCleaning(handles.images.im_bw, opts);
            imwrite(handles.images.im_clean, fileNameClean);
            
            reProcess  = true;
        end

        handles.mask.maskStack{1} = handles.images.im_clean;
        handles.mask.previewMask = zeros(size(handles.images.im_clean));
        handles.mask.currIdx = 1;            
        
        
        handles.images.im_cleanOverlay = doOverlay(handles.images.im_in_color, handles.images.im_clean, [0.8,0.2], handles.mask.alphaColor);
        
        drawImageResultsSingle(hObject, handles.images, handles, 'im_cleanOverlay', 'text_imageDetection', 'imageDetection');
        
        if ~exist(fileNameHull, 'file')
            handles.images.im_hull = doHull(handles.images.im_clean, opts);
            imwrite(handles.images.im_hull, fileNameHull);
        end
        
        % Show Progress Bar        
        progressMultiplier = progressMultiplier + 1;
        waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Extracting Skeleton & Branches ... %2.0f%% along',progressMultiplier*progressFactor*100));
        
        if (exist(fileNameSkelClean, 'file') == 2 && exist(fileNameBranchClean, 'file')) && ~reProcess
            handles.images.im_skel_clean     = imread(fileNameSkelClean);
            handles.images.im_branches_clean = imread(fileNameBranchClean);
            handles.images.im_skel           = imread(fileNameSkel);
            handles.images.im_branches       = imread(fileNameBranch);
        else
            [handles.images.im_skel, handles.images.im_branches, handles.images.im_skel_clean, handles.images.im_branches_clean] = doSkeleton(handles.images.im_clean, opts);

            % removing single isolated pixels from skeleton image using 8-neighbors
            CC = bwconncomp(handles.images.im_skel_clean, 8);
            L = labelmatrix(CC);
            S = regionprops(CC, 'Area');
            BW = ismember(L, find([S.Area] > 1));
            handles.images.im_skel_clean = BW;
            
            if opts.skeleton.cleanBranches == 1
                opts.branches.pruneDist  = get(handles.slider_spatialBranches, 'value');
                handles.images.im_branches_clean = doPruneBranches(handles.images.im_branches_clean, opts);
            end
            
            imwrite(handles.images.im_skel, fileNameSkel);
            imwrite(handles.images.im_branches, fileNameBranch);
            
            imwrite(handles.images.im_skel_clean, fileNameSkelClean);
            imwrite(handles.images.im_branches_clean, fileNameBranchClean);
            
            skel_overlay_del     = doOverlay(handles.images.im_clean, handles.images.im_skel_clean, [0.3,0.2], 'r');
            branch_overlay_del   = doOverlay(skel_overlay_del, imdilate(logical(handles.images.im_branches_clean),strel('disk',5)), [0.7,0.2], 'g');
            imwrite(skel_overlay_del, fileNameSkelCleanOverlay);
            imwrite(branch_overlay_del, fileNameBranchCleanOverlay);
           
            skel_overlay     = doOverlay(handles.images.im_clean, handles.images.im_skel, [0.3,0.2], 'r');
            branch_overlay   = doOverlay(skel_overlay, imdilate(logical(handles.images.im_branches),strel('disk',5)), [0.7,0.2], 'g');
            imwrite(skel_overlay, fileNameSkelOverlay);
            imwrite(branch_overlay, fileNameBranchOverlay);
            
            reProcess = true;
        end
        
        % Generating overlay for skeleton image.
        dilSzeSkel = 1;
        dilSzeBranch = 4;
        
        handles.images.skel_overlay = doOverlay(handles.images.preProc, bwmorph(handles.images.im_skel_clean, 'dilate', dilSzeSkel), [0.9,0.1], skelOverlayColor);
        overlaySkelProc = doOverlay(imInputColor, bwmorph(handles.images.im_skel_clean, 'dilate', dilSzeSkel), [0.9,0.1], skelOverlayColor);
        imwrite(overlaySkelProc, fileNameSkelDelOverlayInputImRGB);
        
        overlaySkelOrig = doOverlay(imInputColor, bwmorph(handles.images.im_skel, 'dilate', dilSzeSkel), [0.9,0.1], skelOverlayColor);
        imwrite(overlaySkelOrig, fileNameSkelOrigOverlayInputImRGB);
        
        handles.images.branches_overlay = doOverlay(handles.images.preProc, bwmorph(handles.images.im_branches_clean, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        overlayBranchClean = doOverlay(imInputColor, bwmorph(handles.images.im_branches_clean, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        imwrite(overlayBranchClean, fileNameGUIBranchesDel);
        overlayBranchOrig = doOverlay(imInputColor, bwmorph(handles.images.im_branches, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        imwrite(overlayBranchOrig, fileNameGUIBranchesOrig);
        
        % Saving images consisting of skeleton + branching overlap.
        orig_skel_branch_overlay = doOverlay(overlaySkelOrig, bwmorph(handles.images.im_branches, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        imwrite(orig_skel_branch_overlay, fileNameSkelBranchOverlayOrig);
        
        proc_skel_branch_overlay = doOverlay(overlaySkelProc, bwmorph(handles.images.im_branches_clean, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        handles.images.final_overlay = proc_skel_branch_overlay;
        imwrite(proc_skel_branch_overlay, fileNameSkelBranchOverlayProc);
        
        drawImageResultsSingle(hObject, handles.images, handles, 'final_overlay', 'text_imageOutput', 'imageOutput');
        
        close(hProgressBar);
        
% 	catch ME
% 		% Some error happened if you get here.
% 		errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
% 			ME.stack(1).name, ME.stack(1).line, ME.message);
% 		warnUser(errorMessage);
%     end
    
	set(gcf,'Pointer','arrow');
	drawnow;	% Cursor won't change right away unless you do this.
        
    return
    

function drawImageResultsSingle(hObject, resIms, handles, variableName, handleTextName, handleImageName, useLabel, labelStr)

    if nargin > 8
        error('myfuns:drawImageResultsSingle:TooManyInputs', ...
            'requires at most 8 optional inputs');
    end
    
    % Fill in unset optional values.
    switch nargin
        case 6
            useLabel = false;
            labelStr = 'Image';
        case 7
            labelStr = 'Image';
    end

    set(eval(['handles.' handleTextName]), 'visible','on');
    imshow(eval(['resIms.' variableName]), 'InitialMagnification', 'fit', 'parent', eval(['handles.' handleImageName]));
    
    
    % get current GUI size
    newPos = get(gcf, 'position');
    newWidth = newPos(3);
    newHeight = newPos(4);
    currWidth = handles.currPos(3);
    currHeight = handles.currPos(4);

    axes(eval(['handles.' handleImageName]));
    fontSize = 31;
    
    if(newWidth > currWidth || newHeight > currHeight)
        newSze = max(newWidth/currWidth, newHeight/currHeight);
        fontSize = fontSize + round(newSze);
    else
        newSze = max(currWidth/newWidth, currHeight/newHeight);
        fontSize = fontSize - round(newSze);
    end
    
    if useLabel
        xlabel(labelStr,'FontSize', fontSize, 'Units','normalized');
    end
    
    drawnow;
    
    handles.currPos = newPos;
    
    % Update handles structure
    guidata(hObject, handles);

    return


% --- Executes on button press in button_Analyze.
function button_Analyze_Callback(hObject, eventdata, handles)
% hObject    handle to button_Analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles = createDirectories(handles);
%     handles = getHandlesDirectories(handles);
    
	set(gcf,'Pointer','watch');
	drawnow;	% Cursor won't change right away unless you do this.
    
	% Get a list of those indexes that are selected, so we know which images to process.
    idxImageFiles   = get(handles.listImageFiles, 'value');
    
	% Then get list of all of the filenames in the list,
	% regardless of whether they are selected or not.
    imageNames = get(handles.listImageFiles, 'string');
    
    if ~isempty(imageNames)
        noFiles = length(idxImageFiles);

        for ii=1:noFiles
            idx = idxImageFiles(ii); 

            baseFilename   = strcat(cell2mat(imageNames(idx)));
            imageFilename  = fullfile(handles.fileFolder, baseFilename);

            handles = analyzeImage(hObject, imageFilename, handles);

            guidata(hObject, handles);
            
        end
    end
    
	set(gcf,'Pointer','arrow');
	drawnow;	% Cursor won't change right away unless you do this.

	guidata(hObject, handles);
    return;
    
    
function pixSze = getPixelSize(handles)
% this function computes the output pixel size in micrometers given the
% microscope parameters.
% handles    structure with handles and user data (see GUIDATA)

    dataMicro = handles.data.microscope;
    %pixSze = dataMicro.cameraPixSize * dataMicro.binning / ( dataMicro.objectiveMagnification * dataMicro.lenseMagnification * dataMicro.CMount );
    pixSze = dataMicro.cameraPixSize;
    
    return;


% --- Executes on button press in button_Statistics.
function button_Statistics_Callback(hObject, eventdata, handles)
% button_Statistics_Callback: this function generates an CSV file
% containing the topological measurements of the processed images.
%
% hObject    handle to button_Statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	set(gcf,'Pointer','watch');
	drawnow;	% Cursor won't change right away unless you do this.
    
    %handles = createDirectories(handles);
    handles = getHandlesDirectories(handles);
    
    imageNames = get(handles.listImageFiles, 'string');
    noFiles    = length(imageNames);
    
	opts    = [];
    opts.skeleton.cleanSpur         = get(handles.checkbox_Spurs,'value');
    opts.skeleton.spurLength        = handles.data.spurLength;
    
    opts.skeleton.cleanBranches     = get(handles.checkbox_SpatialBranches,'value');
    opts.skeleton.spatialBranches   = handles.data.spatialBranches;
    
    stats     = repmat(struct('fileName','', 'tubeArea',-1, 'hullArea', -1, 'skelLengthOrig',-1, 'skelLength',-1, 'noBranchesOrig',-1, 'noBranches',-1, 'tubeWidthOrig', -1, 'tubeWidth', -1), noFiles, 1);
    
    %hProgressBar = waitbar(0, 'Initializing waitbar...', 'Name', 'Progress Bar', 'WindowStyle', 'modal');
    hProgressBar = waitbar(0, 'Initializing waitbar...', 'Name', 'Progress Bar');
    noSteps = noFiles;
	progressFactor = 1/noSteps;
	progressMultiplier = 0;
    
    pixSzeUM = getPixelSize(handles);

    for ii=1:noFiles
        
        [~, imFileName, ~] = fileparts(imageNames{ii});
        imFileName         = strrep(imFileName, ' ', '_');
        
        progressMultiplier = progressMultiplier + 1;
        imFileName4Bar          = strrep(imFileName, '_', '\_');
        waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Obtaining statistics for image %s ... %2.0f%% along', imFileName4Bar, progressMultiplier*progressFactor*100));
        
        [fileNameSkelStrClean, ~]   = getSkelFilename(opts);
        [fileNameBranchStrClean, ~] = getBranchFilename(opts);

        fileNameTubesClean   = fullfile(handles.outCleanDir, sprintf('%s%s%s', imFileName, '_clean', handles.data.fileName.extension));
        fileNameHull         = fullfile(handles.outHullDir, sprintf('%s%s%s', imFileName, '_hull', handles.data.fileName.extension));
        fileNameSkelClean    = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrClean, fileNameBranchStrClean, '_clean', handles.data.fileName.extension));
        fileNameBranchClean  = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrClean, fileNameBranchStrClean, '_clean', handles.data.fileName.extension));        
        fileNameSkelOrig     = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s', imFileName, '_skelcutnone_branchescutnone', handles.data.fileName.extension));
        fileNameBranchOrig   = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s', imFileName, '_skelcutnone_branchescutnone', handles.data.fileName.extension));
        
        stats(ii).fileName          = imFileName;
        
        if (exist(fileNameTubesClean, 'file') == 2) && (exist(fileNameSkelClean, 'file') == 2) && (exist(fileNameBranchClean, 'file') == 2) && (exist(fileNameSkelOrig, 'file') == 2) && (exist(fileNameBranchOrig, 'file') == 2)
            imTubes                     = imread(fileNameTubesClean);
            
            tubeArea                    = numel(find(imTubes));
            stats(ii).tubeArea          = tubeArea;
            
            
            imHull = imread(fileNameHull);
            hullArea = tubeArea./numel(find(imHull));
            stats(ii).hullArea          = hullArea;

            imSkelOrig                  = imread(fileNameSkelOrig);
            skelLengthOrig              = numel(find(imSkelOrig));
            stats(ii).skelLengthOrig    = skelLengthOrig;

            imSkel                      = imread(fileNameSkelClean);
            skelLength                  = numel(find(imSkel));
            stats(ii).skelLength        = skelLength;

            imBranchOrig                = imread(fileNameBranchOrig);
            noBranchesOrig              = numel(find(imBranchOrig));
            stats(ii).noBranchesOrig    = noBranchesOrig;

            imBranch                    = imread(fileNameBranchClean);
            noBranches                  = numel(find(imBranch));
            stats(ii).noBranches        = noBranches;
            
            stats(ii).tubeWidthOrig    =  stats(ii).tubeArea/stats(ii).skelLengthOrig;
            stats(ii).tubeWidth        =  stats(ii).tubeArea/stats(ii).skelLength;
        else
            stats(ii).tubeArea = -1;
            stats(ii).hullArea = -1;
            stats(ii).skelLengthOrig = -1;
            stats(ii).skelLength = -1;
            stats(ii).noBranchesOrig = -1;
            stats(ii).noBranches = -1;
            stats(ii).tubeWidthOrig = -1;
            stats(ii).tubeWidth = -1;
        end
    end
    
    if noFiles
        
%         cameraPixSize = handles.data.microscope.cameraPixSize;
        
        header ={'Image ID', 'Tube Area(pix)', 'Hull Area(pix)', 'Skeleton Length Orig(pix)', 'Skeleton Length(pix)', 'No. of Branches Orig', 'No. of Branches', 'Tube Width Orig(pix)', 'Tube Width(pix)', 'Tube Area(um)', 'Hull Area(um)', 'Skeleton Length Orig(um)', 'Skeleton Length(um)', 'No. of Branches Orig', 'No. of Branches', 'Tube Width Orig(um)', 'Tube Width(um)'};
        fid = fopen(fullfile(handles.outStatsDir, sprintf('stats_%s_%s_%s%s_hull.csv', handles.strParams.tubuStr, handles.strParams.segmentType, fileNameSkelStrClean, fileNameBranchStrClean)), 'w');
        fprintf(fid, '%s,', header{:});
        fprintf(fid, '\n');
        for nf=1:noFiles
%             fprintf(fid, '%s, %d, %2.4f, %d, %d, %d, %d, %2.4f, %2.4f, %d, %2.4f, %d, %d, %d, %d, %2.4f, %2.4f\n', stats(nf).fileName, stats(nf).tubeArea, stats(nf).hullArea, stats(nf).skelLengthOrig, stats(nf).skelLength, stats(nf).noBranchesOrig, stats(nf).noBranches, stats(nf).tubeWidthOrig, stats(nf).tubeWidth, ...
%                 stats(nf).tubeArea*pixSzeUM, stats(nf).hullArea*pixSzeUM, stats(nf).skelLengthOrig*pixSzeUM, stats(nf).skelLength*pixSzeUM, stats(nf).noBranchesOrig, stats(nf).noBranches, stats(nf).tubeWidthOrig*pixSzeUM, stats(nf).tubeWidth*pixSzeUM);
%             fprintf(fid, '%s, %d, %2.4f, %d, %d, %d, %d, %2.4f, %2.4f, %d, %2.4f, %d, %d, %d, %d, %2.4f, %2.4f\n', stats(nf).fileName, stats(nf).tubeArea, stats(nf).hullArea, stats(nf).skelLengthOrig, stats(nf).skelLength, stats(nf).noBranchesOrig, stats(nf).noBranches, stats(nf).tubeWidthOrig, stats(nf).tubeWidth, ...
%                 stats(nf).tubeArea*(cameraPixSize^2), stats(nf).hullArea*pixSzeUM, stats(nf).skelLengthOrig*pixSzeUM, stats(nf).skelLength*pixSzeUM, stats(nf).noBranchesOrig, stats(nf).noBranches, stats(nf).tubeWidthOrig*pixSzeUM, stats(nf).tubeWidth*pixSzeUM);
            fprintf(fid, '%s, %d, %2.4f, %d, %d, %d, %d, %2.4f, %2.4f, %d, %2.4f, %d, %d, %d, %d, %2.4f, %2.4f\n', stats(nf).fileName, stats(nf).tubeArea, stats(nf).hullArea, stats(nf).skelLengthOrig, stats(nf).skelLength, stats(nf).noBranchesOrig, stats(nf).noBranches, stats(nf).tubeWidthOrig, stats(nf).tubeWidth, ...
                stats(nf).tubeArea*(pixSzeUM^2), stats(nf).hullArea*(pixSzeUM^2), stats(nf).skelLengthOrig*pixSzeUM, stats(nf).skelLength*pixSzeUM, stats(nf).noBranchesOrig, stats(nf).noBranches, stats(nf).tubeWidthOrig*pixSzeUM, stats(nf).tubeWidth*pixSzeUM);
        end
        fclose(fid);
    end

    close(hProgressBar);
    
	set(gcf,'Pointer','arrow');
	drawnow;	% Cursor won't change right away unless you do this.

% 	guidata(hObject, handles);
    return;

    
 % % Update Analyze button and tooltip string depending on how many files in the listbox were selected.
function updateAnalyzeButtonCaption(handles)

    selIms = get(handles.listImageFiles, 'value');
    if length(selIms) > 1 
        buttonCaption = {'Analyze '};   % MATLAB quirk - needs to be cell array to keep trailing spaces.
        buttonCaption = strcat(buttonCaption, num2str(length(selIms)));
        buttonCaption = strcat(buttonCaption, ' images');
        set(handles.button_Analyze, 'string', buttonCaption);
		set(handles.button_Analyze, 'Tooltipstring', 'Analyze the selected image(s)');
    elseif length(selIms) == 1
        set(handles.button_Analyze, 'string', 'Analyze 1 image');
		set(handles.button_Analyze, 'Tooltipstring', 'Analyze the selected image');	
    else
        set(handles.button_Analyze, 'string', 'No images selected');
		set(handles.button_Analyze, 'Tooltipstring', 'Please select the input image(s)');
    end
    
    return;

% --- Executes on selection change in listImageFiles.
function listImageFiles_Callback(hObject, eventdata, handles)
% hObject    handle to listImageFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listImageFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listImageFiles

	set(gcf,'Pointer','watch');
	drawnow;	% Cursor won't change right away unless you do this.
    
    % Update the number of images in the Analyze button caption.
    updateAnalyzeButtonCaption(handles);
    
	% Get image name
    selectedIm = get(handles.listImageFiles, 'value');
	if isempty(selectedIm)
		% Bail out if nothing was selected.
		% Change mouse pointer (cursor) to an arrow.
		set(gcf,'Pointer','arrow');
		drawnow;	% Cursor won't change right away unless you do this.
		return;
	end
    if length(selectedIm) > 1 
		% Change mouse pointer (cursor) to an arrow.
		set(gcf,'Pointer','arrow')
		drawnow;	% Cursor won't change right away unless you do this.
        return;
    end
    
    % If only one is selected, display it.
    listImageNames = get(handles.listImageFiles, 'string');
    imageFilename = strcat(cell2mat(listImageNames(selectedIm)));
    fullImageFilename = [handles.fileFolder '/' imageFilename];	% Prepend folder.

%     set(handles.imageFilename, 'string', imageFilename);
    fprintf('# %s\n', fullImageFilename);
    
    %# ADD HERE CODE TO VISUALIZE IMAGE SELECTED ON MOUSE CLICK.
	opts = getHandlesValues(handles);
    
    [fileNameSkelStrClean, ~] = getSkelFilename(opts);
    [fileNameBranchStrClean, ~] = getBranchFilename(opts);
    
    [~,imFileName,~] = fileparts(fullImageFilename);
    imFileName = strrep(imFileName, ' ', '_');

    fileNameInputImRGB = fullfile(handles.outInputDir, sprintf('%s%s_rgb%s', imFileName, handles.data.fileName.input, handles.data.fileName.extension));
    if exist(fileNameInputImRGB, 'file')
        imInputColor = imread(fileNameInputImRGB);
        handles.images.im_in_color = imInputColor;
        %drawImageResultsSingle(handles.images, handles, 'im_in_color', 'text_inputImage', 'imageInput');
        drawImageResultsSingle(hObject, handles.images, handles, 'im_in_color', 'text_inputImage', 'imageInput');
    else
        set(handles.imageInput, 'visible','off');
        set(get(handles.imageInput, 'children'), 'visible', 'off');
    end

    fileNameClean           = fullfile(handles.outCleanDir, sprintf('%s%s%s', imFileName, handles.data.fileName.clean, handles.data.fileName.extension));
    if exist(fileNameClean, 'file')
        handles.images.im_clean = imread(fileNameClean);
        
        handles.mask.maskStack{1} = handles.images.im_clean;
        [irows, icols, ~] = size(handles.images.im_clean);
        handles.mask.previewMask = zeros(size(irows, icols));
        handles.mask.currIdx = 1;
        
        handles.images.im_cleanOverlay = doOverlay(handles.images.im_in_color, handles.images.im_clean, [0.8,0.2], handles.mask.alphaColor);
        drawImageResultsSingle(hObject, handles.images, handles, 'im_cleanOverlay', 'text_imageDetection', 'imageDetection');
    else
        set(handles.imageDetection, 'visible','off');
        set(get(handles.imageDetection, 'children'), 'visible', 'off');
    end

    fileNameSkelBranchOverlayProc     = fullfile(handles.outOverlayDir, sprintf('%s_proc%s%s%s', imFileName, fileNameSkelStrClean, fileNameBranchStrClean, handles.data.fileName.extension));
    if exist(fileNameSkelBranchOverlayProc, 'file')
        handles.images.final_overlay = imread(fileNameSkelBranchOverlayProc);
        drawImageResultsSingle(hObject, handles.images, handles, 'final_overlay', 'text_imageOutput', 'imageOutput');
    else
        set(handles.imageOutput, 'visible', 'off');
        set(get(handles.imageOutput, 'children'), 'visible', 'off');
    end
    %# FINISH HERE CODE TO VISUALIZE IMAGE SELECTED ON MOUSE CLICK.
    
    % Change mouse pointer (cursor) to an arrow.
	set(gcf,'Pointer','arrow');
	drawnow;	% Cursor won't change right away unless you do this.
    
    %set(hListBox, 'String', listString);
    guidata(hObject, handles);
    return;


% --- Executes during object creation, after setting all properties.
function listImageFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listImageFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in checkbox_Illumination.
function checkbox_Illumination_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Illumination (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Illumination


% --- Executes on button press in checkbox_Denoising.
function checkbox_Denoising_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Denoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Denoising



function textEdit_denoisingBM3D_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_denoisingBM3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_denoisingBM3D as text
%        str2double(get(hObject,'String')) returns contents of textEdit_denoisingBM3D as a double


% --- Executes during object creation, after setting all properties.
function textEdit_denoisingBM3D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_denoisingBM3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on selection change in popupmenuAdjustment.
function popupmenuAdjustment_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuAdjustment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuAdjustment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuAdjustment

    adjustmentStr = get(hObject, 'String');
    adjustmentVal = get(hObject,'Value');
    handles.adjustment.type = lower(adjustmentStr{adjustmentVal});
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenuAdjustment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuAdjustment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on selection change in popupmenuDenoise.
function popupmenuDenoise_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuDenoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuDenoise contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuDenoise

    denoiseStr = get(hObject, 'String');
    denoiseVal = get(hObject, 'Value');
    handles.denoising.type = lower(denoiseStr{denoiseVal});
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenuDenoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuDenoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in checkbox_Adjustment.
function checkbox_Adjustment_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Adjustment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Adjustment



function textEdit_IlluminationDiskSze_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_IlluminationDiskSze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_IlluminationDiskSze as text
%        str2double(get(hObject,'String')) returns contents of textEdit_IlluminationDiskSze as a double


% --- Executes during object creation, after setting all properties.
function textEdit_IlluminationDiskSze_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_IlluminationDiskSze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in checkbox_Edges.
function checkbox_Edges_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Edges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Edges


% --- Executes on button press in checkbox_Tubularity.
function checkbox_Tubularity_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Tubularity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Tubularity



function textEdit_Tubularity_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_Tubularity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_Tubularity as text
%        str2double(get(hObject,'String')) returns contents of textEdit_Tubularity as a double


% --- Executes during object creation, after setting all properties.
function textEdit_Tubularity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_Tubularity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in checkbox_SmallRegions.
function checkbox_SmallRegions_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_SmallRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_SmallRegions



function textEdit_RegionAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_RegionAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_RegionAnalysis as text
%        str2double(get(hObject,'String')) returns contents of textEdit_RegionAnalysis as a double


% --- Executes during object creation, after setting all properties.
function textEdit_RegionAnalysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_RegionAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on selection change in popupmenuThreshold.
function popupmenuThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuThreshold contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuThreshold

    segmentStr = get(hObject, 'String');
    segmentVal = get(hObject, 'Value');
    handles.segmentation.type = lower(segmentStr{segmentVal});
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenuThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in checkbox_Holes.
function checkbox_Holes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Holes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Holes



function textEdit_FillHoles_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_FillHoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_FillHoles as text
%        str2double(get(hObject,'String')) returns contents of textEdit_FillHoles as a double


% --- Executes during object creation, after setting all properties.
function textEdit_FillHoles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_FillHoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in checkbox_Spurs.
function checkbox_Spurs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Spurs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Spurs


% --- Executes on button press in checkbox_SpatialBranches.
function checkbox_SpatialBranches_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_SpatialBranches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_SpatialBranches


% --- Executes on slider movement.
function slider_spurLength_Callback(hObject, eventdata, handles)
% hObject    handle to slider_spurLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    handles.data.spurLength  = round(get(hObject,'Value'));
    set(handles.text_spurLength,'String', sprintf('Spur Length: %d', handles.data.spurLength));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_spurLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_spurLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end


% --- Executes on slider movement.
function slider_spatialBranches_Callback(hObject, eventdata, handles)
% hObject    handle to slider_spatialBranches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    handles.data.spatialBranches  = round(get(hObject,'Value'));
    set(handles.text_spatialBranches,'String', sprintf('Spatial Distance: %d', handles.data.spatialBranches));
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider_spatialBranches_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_spatialBranches (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end


% --- Executes on selection change in popupmenuColorChannel.
function popupmenuColorChannel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuColorChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuColorChannel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuColorChannel


% --- Executes during object creation, after setting all properties.
function popupmenuColorChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuColorChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in checkbox_convHull.
function checkbox_convHull_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_convHull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_convHull


% --- Executes on slider movement.
function slider_brushSize_Callback(hObject, eventdata, handles)
% hObject    handle to slider_brushSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    handles.data.brushRad  = round(get(hObject,'Value'));
    set(handles.text_brushSize,'String', sprintf('Brush size: %d', handles.data.brushRad));
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_brushSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_brushSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end


% --- Executes on button press in button_RedoMask.
function button_RedoMask_Callback(hObject, eventdata, handles)
% hObject    handle to button_RedoMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.mask.currIdx = handles.mask.currIdx + 1;
    if handles.mask.currIdx >= length(handles.mask.maskStack)
        handles.mask.currIdx  = length(handles.mask.maskStack);
    end
    
    drawOverlayMask(hObject, eventdata, handles);
    guidata(hObject,handles);

% --- Executes on button press in button_UndoMask.
function button_UndoMask_Callback(hObject, eventdata, handles)
% hObject    handle to button_UndoMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.mask.currIdx = handles.mask.currIdx - 1;
    if handles.mask.currIdx <= 1
        handles.mask.currIdx  = 1;
    end
    
    drawOverlayMask(hObject, eventdata, handles);
    guidata(hObject,handles);

% --- Executes on button press in button_RedoSkel.
function button_RedoSkel_Callback(hObject, eventdata, handles)
% button_RedoSkel_Callback: this function recomputes the skeleton, if the
% used desires it so.
%
% hObject    handle to button_RedoSkel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    opts = getHandlesValues(handles);
    
    imFileNames = get(handles.listImageFiles,'String');
    
	set(gcf,'Pointer','watch');
	drawnow;	% Cursor won't change right away unless you do this.
    
    %hProgressBar = waitbar(0, 'New skeleton and branching points ...', 'Name', 'Progress Bar', 'WindowStyle', 'modal');
    hProgressBar = waitbar(0, 'New skeleton and branching points ...', 'Name', 'Progress Bar');
	noSteps = 4;
	progressFactor = 1/noSteps;
	progressMultiplier = 0;
    
    if isfield(handles.mask, 'maskStack')
        
        imChannelIdx = get(handles.popupmenuColorChannel, 'Value');

        if imChannelIdx == 2
            skelOverlayColor = 'r';
        else
            skelOverlayColor = 'g';
        end
        branchOverlayColor = 'y';

        if opts.skeleton.spurLength == 0 && opts.skeleton.cleanSpur == 1
            fileNameGUISkelStr      = sprintf('%scutall', handles.data.fileName.skel);
        else
            fileNameGUISkelStr      = sprintf('%scut%03d', handles.data.fileName.skel, opts.skeleton.spurLength);
        end
        
        selectedIdx = get(handles.listImageFiles,'Value');
        imageFileName = imFileNames{selectedIdx};
        
        [fileNameSkelStrClean, fileNameSkelStrCleanOverlay] = getSkelFilename(opts);
        fileNameSkelStrOverlay      = sprintf('_skelcutnone');

        [fileNameBranchStrClean, fileNameBranchStrCleanOverlay] = getBranchFilename(opts);
        fileNameBranchStrOverlay       = sprintf('_branchcutnone');

        [~,imFileName,~] = fileparts(imageFileName);
        imFileName = strrep(imFileName, ' ', '_');
        
        fileNamePreproc         = fullfile(handles.outPreprocDir, sprintf('%s%s%s', imFileName, handles.data.fileName.preproc, handles.data.fileName.extension));
        fileNameInputImRGB = fullfile(handles.outInputDir, sprintf('%s%s_rgb%s', imFileName, handles.data.fileName.input, handles.data.fileName.extension));

        fileNameSkel            = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s%s', imFileName, sprintf('%scutnone', handles.data.fileName.skel), sprintf('%scutnone', handles.data.fileName.branch), handles.data.fileName.extension));
        fileNameBranch          = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s%s', imFileName, sprintf('%scutnone', handles.data.fileName.skel), sprintf('%scutnone', handles.data.fileName.branch), handles.data.fileName.extension));
        fileNameSkelClean    = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrClean, fileNameBranchStrClean, '_clean', handles.data.fileName.extension));
        fileNameBranchClean  = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrClean, fileNameBranchStrClean, '_clean', handles.data.fileName.extension)); 
        fileNameSkelCleanOverlay    = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrCleanOverlay, fileNameBranchStrCleanOverlay, '_clean_overlay', handles.data.fileName.extension));
        fileNameBranchCleanOverlay  = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrCleanOverlay, fileNameBranchStrCleanOverlay, '_clean_overlay', handles.data.fileName.extension));
        fileNameSkelOverlay         = fullfile(handles.outSkelDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrOverlay, fileNameBranchStrOverlay, '_overlay', handles.data.fileName.extension));
        fileNameBranchOverlay       = fullfile(handles.outBranchesDir, imFileName, sprintf('%s%s%s%s%s', imFileName, fileNameSkelStrOverlay, fileNameBranchStrOverlay, '_overlay', handles.data.fileName.extension));
        
        fileNameClean           = fullfile(handles.outCleanDir, sprintf('%s%s%s', imFileName, handles.data.fileName.clean, handles.data.fileName.extension));
		fileNameSkelDelOverlayInputImRGB = fullfile(handles.outOverlayDir, imFileName, sprintf('%s%s%s%s_rgb_skeldel_input%s', imFileName, handles.data.fileName.skel, fileNameGUISkelStr, ...
                                                 sprintf('%scut%03d',handles.data.fileName.branch, opts.skeleton.spatialBranches), handles.data.fileName.extension));
		fileNameSkelOrigOverlayInputImRGB     = fullfile(handles.outOverlayDir, imFileName, sprintf('%s%s_rgb_skeldel_input%s', imFileName, sprintf('%s_orig', handles.data.fileName.skel), handles.data.fileName.extension));
		fileNameGUIBranchesDel  = fullfile(handles.outOverlayDir, imFileName, sprintf('%s%s%s%s%s', imFileName, handles.data.fileName.branch, fileNameGUISkelStr, sprintf('%scut%03d',handles.data.fileName.branch, opts.skeleton.spatialBranches), handles.data.fileName.extension));
		fileNameGUIBranchesOrig = fullfile(handles.outOverlayDir, imFileName, sprintf('%s%s%s', imFileName, sprintf('%s_orig', handles.data.fileName.branch), handles.data.fileName.extension));
		fileNameSkelBranchOverlayOrig     = fullfile(handles.outOverlayDir, sprintf('%s_orig%s', imFileName, handles.data.fileName.extension));
		fileNameSkelBranchOverlayProc     = fullfile(handles.outOverlayDir, sprintf('%s_proc%s%s%s', imFileName, fileNameSkelStrClean, fileNameBranchStrClean, handles.data.fileName.extension));
    
        %fprintf('# curr-id:%d\n', handles.mask.currIdx);
        imClean = handles.mask.maskStack{handles.mask.currIdx};

        imwrite(imClean, fileNameClean);
        
        progressMultiplier = progressMultiplier + 1;
        waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Recomputing skeleton and branching points ... %2.0f%% along', progressMultiplier*progressFactor*100));
        
        [handles.images.im_skel, handles.images.im_branches, handles.images.im_skel_clean, handles.images.im_branches_clean] = doSkeleton(imClean, opts);
        
        % removing single isolated pixels from skeleton image using 8-neighbors
        CC = bwconncomp(handles.images.im_skel_clean, 8);
        L = labelmatrix(CC);
        S = regionprops(CC, 'Area');
        BW = ismember(L, find([S.Area] > 1));
        handles.images.im_skel_clean = BW;

        if opts.skeleton.cleanBranches == 1
            
            progressMultiplier = progressMultiplier + 1;
            waitbar(progressMultiplier*progressFactor, hProgressBar, sprintf('Recomputing skeleton and branching points ... %2.0f%% along', progressMultiplier*progressFactor*100));
        
            opts.branches.pruneDist  = get(handles.slider_spatialBranches, 'value');
            handles.images.im_branches_clean = doPruneBranches(handles.images.im_branches_clean, opts);
        end

        imwrite(handles.images.im_skel, fileNameSkel);
        imwrite(handles.images.im_branches, fileNameBranch);

        imwrite(handles.images.im_skel_clean, fileNameSkelClean);
        imwrite(handles.images.im_branches_clean, fileNameBranchClean);

        skel_overlay_del     = doOverlay(handles.images.im_clean, handles.images.im_skel_clean, [0.3,0.2], 'r');
        branch_overlay_del   = doOverlay(skel_overlay_del, imdilate(logical(handles.images.im_branches_clean),strel('disk',5)), [0.7,0.2], 'g');
        imwrite(skel_overlay_del, fileNameSkelCleanOverlay);
        imwrite(branch_overlay_del, fileNameBranchCleanOverlay);

        skel_overlay     = doOverlay(handles.images.im_clean, handles.images.im_skel, [0.3,0.2], 'r');
        branch_overlay   = doOverlay(skel_overlay, imdilate(logical(handles.images.im_branches),strel('disk',5)), [0.7,0.2], 'g');
        imwrite(skel_overlay, fileNameSkelOverlay);
        imwrite(branch_overlay, fileNameBranchOverlay);
        
        % Generating overlay for skeleton image.
        dilSzeSkel = 1;
        dilSzeBranch = 4;
        
        handles.images.preProc = imread(fileNamePreproc);
        imInputColor = imread(fileNameInputImRGB);
        handles.images.skel_overlay = doOverlay(handles.images.preProc, bwmorph(handles.images.im_skel_clean, 'dilate', dilSzeSkel), [0.9,0.1], skelOverlayColor);
        overlaySkelProc = doOverlay(imInputColor, bwmorph(handles.images.im_skel_clean, 'dilate', dilSzeSkel), [0.9,0.1], skelOverlayColor);
        imwrite(overlaySkelProc, fileNameSkelDelOverlayInputImRGB);
        
        overlaySkelOrig = doOverlay(imInputColor, bwmorph(handles.images.im_skel, 'dilate', dilSzeSkel), [0.9,0.1], skelOverlayColor);
        imwrite(overlaySkelOrig, fileNameSkelOrigOverlayInputImRGB);
        
        handles.images.branches_overlay = doOverlay(handles.images.preProc, bwmorph(handles.images.im_branches_clean, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        overlayBranchClean = doOverlay(imInputColor, bwmorph(handles.images.im_branches_clean, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        imwrite(overlayBranchClean, fileNameGUIBranchesDel);
        overlayBranchOrig = doOverlay(imInputColor, bwmorph(handles.images.im_branches, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        imwrite(overlayBranchOrig, fileNameGUIBranchesOrig);
        
        % Saving images consisting of skeleton + branching overlap.
        orig_skel_branch_overlay = doOverlay(overlaySkelOrig, bwmorph(handles.images.im_branches, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        imwrite(orig_skel_branch_overlay, fileNameSkelBranchOverlayOrig);
        
        proc_skel_branch_overlay = doOverlay(overlaySkelProc, bwmorph(handles.images.im_branches_clean, 'dilate', dilSzeBranch), [0.9,0.1], branchOverlayColor);
        handles.images.final_overlay = proc_skel_branch_overlay;
        imwrite(proc_skel_branch_overlay, fileNameSkelBranchOverlayProc);
        
        %drawImageResultsSingle(handles.images, handles, 'final_overlay', 'text_imageOutput', 'imageOutput');
        drawImageResultsSingle(hObject, handles.images, handles, 'final_overlay', 'text_imageOutput', 'imageOutput');
    end
    
    close(hProgressBar);
    
	set(gcf,'Pointer','arrow');
	drawnow;



function textEdit_ObjectiveMagnification_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_ObjectiveMagnification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_ObjectiveMagnification as text
%        str2double(get(hObject,'String')) returns contents of textEdit_ObjectiveMagnification as a double


% --- Executes during object creation, after setting all properties.
function textEdit_ObjectiveMagnification_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_ObjectiveMagnification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function textEdit_LenseMagnification_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_LenseMagnification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_LenseMagnification as text
%        str2double(get(hObject,'String')) returns contents of textEdit_LenseMagnification as a double


% --- Executes during object creation, after setting all properties.
function textEdit_LenseMagnification_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_LenseMagnification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function textEdit_CMount_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_CMount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_CMount as text
%        str2double(get(hObject,'String')) returns contents of textEdit_CMount as a double


% --- Executes during object creation, after setting all properties.
function textEdit_CMount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_CMount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function textEdit_CameraPixelSize_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_CameraPixelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_CameraPixelSize as text
%        str2double(get(hObject,'String')) returns contents of textEdit_CameraPixelSize as a double


% --- Executes during object creation, after setting all properties.
function textEdit_CameraPixelSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_CameraPixelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function textEdit_Binning_Callback(hObject, eventdata, handles)
% hObject    handle to textEdit_Binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textEdit_Binning as text
%        str2double(get(hObject,'String')) returns contents of textEdit_Binning as a double


% --- Executes during object creation, after setting all properties.
function textEdit_Binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textEdit_Binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
