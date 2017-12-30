function [ data, wavelengths, params ] = loadSPE( filename )
%% LOADSPE load LightField SPE v3.0 or WinSpec SPE v2.x data from a file
% [data, wavelengths, params] = loadSPE( filename )
%
%% *data*
% If the file has only one ROI and one frame then data is a simple (2D) array.
% If the file has one ROI but several frames then data is 3D array of
% dimensions ( height, width, number of frames). 
%
% If the file has several ROIs then the data variable will look different
% for v2.x and v3.x spe files due the difference in how the data is stored.
% 
% For v3.x if the file has one frame and several ROIs data is a struct 
% data.ROI{numOfROIs}( height, width). If the file has several frames and 
% several ROIs in each frame then data is  an array of structs 
% data(numOfFrames).ROI{numOfROIs}( height, width).
%
% For v2.x the data will be either a 2D or 3D array, depending on the
% number of frames. Then params.ROI will contain a struct with info on
% sizes of ROIs and location on the full frame, which then can be
% individually extracted. This may be added in the future releases
%
%% *params*
% params.version contains version info of the file. For SPE v3.x the whole
% XML footer is parsed as a struct in params.SpeFormat. 
% SPE v3.0 file and XML specification is at <a href="matlab:web('ftp://ftp.princetoninstruments.com/public/Manuals/Princeton%20Instruments/SPE%203.0%20File%20Format%20Specification.pdf')">the PI website</a>.
%
% For the SPE v2.x all parameters were stored in file header some of the
% key paramters are stored directly in the 'root' of the params struct:
% params.ROI
% params.xdim
% params.ydim
% params.xlabel
% params.ylabel
% params.dlabel
% params.SpecGrooves
% params.ExperimentTimeLocal
% params.date
% params.exp_sec
%
% Full set of parameters are stored in params.full except for x and y
% calibrations which are stored separately as params.xcalib and
% params.ycalib. Full specification of the SPE v2.x header is avaliable at
% <a href="matlab:web('ftp://ftp.piacton.com/Public/Manuals/Princeton%20Instruments/WinSpec%202.6%20Spectroscopy%20Software%20User%20Manual.pdf')">the PI website</a>.
% 
%% *wavelengths*
% For SPE v2.x is a single 1D array defined by either x-axis or y-axis
% calibration, whichever is present in the file using the polynomial
% method.
% For SPE v3.x wavelegths is a 1D array is thereis a single ROI or a cell
% array of 1D arrays corresponding to the different ROIs
%
%% Examples:
% - 1 frame and 1 ROI per frame
%   [ d, w, ~] = loadSPE( 'example.spe' );
%   plot( w, d );
%
% - 1 frame and several ROIs per frame. Some ROIs are 1D spectra, some are
%   2D images of CCD
%   [ d, w, ~] = loadSPE( 'example2.spe' );
%   for i = 1 : length( w )
%       [ x, y ] = d.ROI{i}';
%       if y == 1 %1D spectra
%           plot( w{i}, d.ROI{i} );
%       else %2D image
%           imagesc( w{i}, 1:y, d.ROI{i}' );
%           % or could simply plot pixels imagesc( d.ROI{i}');
%       end
%   end
%
%% Copyright 
% XML parsing is implemented using functions based on <a href="matlab:web('https://mathworks.com/matlabcentral/fileexchange/28518-xml2struct')">xml2struct</a>,
% written by W. Falkena, modified by A. Wanner, I. Smirnov, X. Mo.
% xml2struct (C) 2012, W. Falkena
%
% (C) 2017, M. Sich, The University of Sheffield
% v2.0 30-12-2017

%% Main code

fid = fopen( filename, 'rb', 'ieee-le' );

% get SPE file version first
fseek( fid, 1992, 'bof' );
params.version = fread( fid, 1, 'single' );

if params.version >= 3
    % SPE v3.x uses XML in the end of the file instead of the header to
    % store information. Although some info might be dublicated in the
    % header section, extract all key parameters from the XML footer.
    fseek( fid, 678, 'bof' );
    footerOffset = fread( fid, 1, 'uint64' );
    fseek( fid, footerOffset, 'bof' );
    XMLFooter = join( string( fread( fid, Inf, 'int8=>char' )), '');
    
    % try to parse the xml string
    % more info https://undocumentedmatlab.com/blog/parsing-xml-strings
    try
        % The following avoids the need for file I/O:
        inputObject = java.io.StringBufferInputStream(XMLFooter);  % or: org.xml.sax.InputSource(java.io.StringReader(xmlString))
        try
            % Parse the input data directly using xmlread's core functionality
            parserFactory = javaMethod('newInstance','javax.xml.parsers.DocumentBuilderFactory');
            p = javaMethod('newDocumentBuilder',parserFactory);
            xmlTreeObject = p.parse(inputObject);
        catch
            % Use xmlread's semi-documented inputObject input feature
            xmlTreeObject = xmlread(inputObject);
        end
    catch
        % Fallback to standard xmlread usage, using a temporary XML file:
        % Store the XML data in a temp *.xml file
        filename = [tempname '.xml'];
        fid = fopen(filename,'Wt');
        fwrite(fid,xmlString);
        fclose(fid);
        % Read the file into an XML model object
        xmlTreeObject = xmlread(filename);
        % Delete the temp file
        delete(filename);
    end
    temp = parseChildNodes(xmlTreeObject);
    params.SpeFormat = temp.SpeFormat;
    
    % to save some typing
    db = params.SpeFormat.DataFormat.DataBlock;
    
    % check pixel format
    switch db.Attributes.pixelFormat
        case 'MonochromeUnsigned16'
            strDataType = 'uint16';
        case 'MonochromeUnsigned32'
            strDataType = 'uint32';
        case 'MonochromeFloating32'
            strDataType = 'single';
        otherwise
            fclose( fid );
            error( ['Unsupported pixel format: ' db.Attributes.pixelFormat]);
    end
    
    % number of frames
    nFrames = str2num(db.Attributes.count);
    % get number of ROIs in each frame
    nROI = length( db.DataBlock );
    fseek( fid, 4100, 'bof' );
    if nROI == 1
        if nFrames == 1
            height = str2num( db.DataBlock.Attributes.height );
            width = str2num( db.DataBlock.Attributes.width );
            data = fread( fid, [ height, width ], strDataType );
        else
            for j = 1 : nFrames
                height = str2num( db.DataBlock.Attributes.height );
                width = str2num( db.DataBlock.Attributes.width );
                data( :, :, j ) = fread( fid, [ height, width ], strDataType );
            end
        end
    else
        if nFrames == 1
            for i = 1 : nROI
                height = str2num( db.DataBlock{1,i}.Attributes.height );
                width = str2num( db.DataBlock{1,i}.Attributes.width );
                data.ROI{i}(:, :) = fread( fid, [ height, width ], strDataType );
            end
        else
            for i = 1 : nROI
                for j = 1 : nFrames
                    height = str2num( db.DataBlock{1,i}.Attributes.height );
                    width = str2num( db.DataBlock{1,i}.Attributes.width );
                    data(j).ROI{i}(:, :) = fread( fid, [ height, width ], strDataType );
                end
            end
        end
    end
    
    % export wavelengths data to a separate parameter
    w = str2num(params.SpeFormat.Calibrations.WavelengthMapping.Wavelength.Text);
    if nROI == 1
        % if there is only one ROI then wavelengths should be a single 1D
        % array
        x1 = str2num(params.SpeFormat.Calibrations.SensorMapping.Attributes.x);
        x2 = x1 + str2num(params.SpeFormat.Calibrations.SensorMapping.Attributes.width);
        wavelengths = w( x1:x2 );
    else
        % if there are several ROIs in the file then create cell array of
        % 1D arrays
        for i = 1 : nROI
            x1 = str2num(params.SpeFormat.Calibrations.SensorMapping{1,i}.Attributes.x);
            x2 = x1 - 1 + str2num(params.SpeFormat.Calibrations.SensorMapping{1,i}.Attributes.width);
            wavelengths{i} = w( x1:x2 );
        end
    end
    
elseif params.version >= 2
    % SPE v2.5
    % All info is contained in the header.
    % Start of Header Information (0 - 2996)
    fseek( fid, 0, 'bof' );
    p.ControllerVersion    = fread( fid, 1, 'int16' );
    p.LogicOutput          = fread( fid, 1, 'int16' );
    p.AmpHiCapLowNoise     = fread( fid, 1, 'uint16' );
    p.xDimDet              = fread( fid, 1, 'uint16' );
    p.mode                 = fread( fid, 1, 'int16' );
    p.exp_sec              = fread( fid, 1, 'single' );
    p.VChipXdim            = fread( fid, 1, 'int16' );
    p.VChipYdim            = fread( fid, 1, 'int16' );
    p.yDimDet              = fread( fid, 1, 'uint16' ); %y dimension of CCD or detector
    p.date                 = fread( fid, 10, 'char' ); %date
    p.date                 = join(string(char(p.date)), '');
    p.VirtualChipFlag      = fread( fid, 1, 'int16' ); %On/Off
    fread( fid, 2, 'char' );    %skipping spare chars
    p.noscan               = fread( fid, 1, 'int16' ); %Old number of scans - should always be -1
    p.DetTemperature       = fread( fid, 1, 'single' ); %Detector Temperature Set
    p.DetType              = fread( fid, 1, 'int16' ); %CCD/DiodeArray type
    p.xdim                 = fread( fid, 1, 'uint16' ); %actual # of pixels on x axis
    p.stdiode              = fread( fid, 1, 'int16' ); %trigger diode
    p.DelayTime            = fread( fid, 1, 'single' ); %46 Used with Async Mode
    p.ShutterControl       = fread( fid, 1, 'uint16' ); %50 Normal, Disabled Open, Disabled Closed
    p.AbsorbLive           = fread( fid, 1, 'int16' ); %52 On/Off
    p.AbsorbMode           = fread( fid, 1, 'uint16' ); %54 Reference Strip or File
    p.CanDoVirtualChipFlag = fread( fid, 1, 'int16' ); %56 T/F Cont/Chip able to do Virtual Chip
    p.ThresholdMinLive     = fread( fid, 1, 'int16' ); %58 On/Off
    p.ThresholdMinVal      = fread( fid, 1, 'single' ); %60 Threshold Minimum Value
    p.ThresholdMaxLive     = fread( fid, 1, 'int16' ); %64 On/Off
    p.ThresholdMaxVal      = fread( fid, 1, 'single' ); %66 Threshold Maximum Value
    p.SpecAutoSpectroMode  = fread( fid, 1, 'int16' ); %70 T/F Spectrograph Used
    p.SpecCenterWlNm       = fread( fid, 1, 'single' ); %72 Center Wavelength in Nm
    p.SpecGlueFlag         = fread( fid, 1, 'int16' ); %T/F File is Glued
    p.SpecGlueStartWlNm    = fread( fid, 1, 'single' ); %78 Starting Wavelength in Nm
    p.SpecGlueEndWlNm      = fread( fid, 1, 'single' ); %82 Starting Wavelength in Nm
    p.SpecGlueMinOvrlpNm   = fread( fid, 1, 'single' ); %86 Minimum Overlap in Nm
    p.SpecGlueFinalResNm   = fread( fid, 1, 'single' ); %90 Final Resolution in Nm
    p.PulserType           = fread( fid, 1, 'int16' ); %94 0=None, PG200=1, PTG=2, DG535=3
    p.CustomChipFlag       = fread( fid, 1, 'int16' ); %96 T/F Custom Chip Used
    p.XPrePixels           = fread( fid, 1, 'int16' ); %98 Pre Pixels in X direction
    p.XPostPixels          = fread( fid, 1, 'int16' ); %100 Post Pixels in X direction
    p.YPrePixels           = fread( fid, 1, 'int16' ); %102 Pre Pixels in Y direction
    p.YPostPixels          = fread( fid, 1, 'int16' ); %104 Post Pixels in Y direction
    p.asynen               = fread( fid, 1, 'int16' ); %106 asynchron enable flag 0 = off
    p.datatype             = fread( fid, 1, 'int16' ); %108 experiment datatype
    switch p.datatype
        case 0
            p.datatype = 'single';
        case 1
            p.datatype = 'int32';
        case 2
            p.datatype = 'int16';
        case 3
            p.datatype = 'uint16';
        otherwise
            temp = p.datatype;
            p.datatype = [ 'Unknown :' num2str(temp)];
            fclose( fid );
            error( ['Unknown .SPE v2.x data type: ' num2str(temp)] );
    end
    p.PulserMode           = fread( fid, 1, 'int16' ); %110 Repetitive/Sequential
    p.PulserOnChipAccums   = fread( fid, 1, 'uint16' ); %112 Num PTG On-Chip Accums
    p.PulserRepeatExp      = fread( fid, 1, 'uint32' ); %114 Num Exp Repeats (Pulser SW Accum)
    p.PulseRepWidth        = fread( fid, 1, 'single' ); %118 Width Value for Repetitive pulse (usec)
    p.PulseRepDelay        = fread( fid, 1, 'single' ); %122 Width Value for Repetitive pulse (usec)
    p.PulseSeqStartWidth   = fread( fid, 1, 'single' ); %126 Start Width for Sequential pulse (usec)
    p.PulseSeqEndWidth     = fread( fid, 1, 'single' ); %130 End Width for Sequential pulse (usec)
    p.PulseSeqStartDelay   = fread( fid, 1, 'single' ); %134 Start Delay for Sequential pulse (usec)
    p.PulseSeqEndDelay     = fread( fid, 1, 'single' ); %138 End Delay for Sequential pulse (usec)
    p.PulseSeqIncMode      = fread( fid, 1, 'int16' ); %142 Increments: 1=Fixed, 2=Exponential
    p.PImaxUsed            = fread( fid, 1, 'int16' ); %144 PI-Max type controller flag
    p.PImaxMode            = fread( fid, 1, 'int16' ); %146 PI-Max mode
    p.PImaxGain            = fread( fid, 1, 'int16' ); %148 PI-Max Gain
    p.BackGrndApplied      = fread( fid, 1, 'int16' ); %150 1 if background subtraction done
    p.PImax2nsBrdUsed      = fread( fid, 1, 'int16' ); %152 T/F PI-Max 2ns Board Used
    p.minblk               = fread( fid, 1, 'uint16' ); %154 min. # of strips per skips
    p.numminblk            = fread( fid, 1, 'uint16' ); %156 # of min-blocks before geo skps
    p.SpecMirrorLocation(:)= fread( fid, 2, 'int16' ); % 158 Spectro Mirror Location, 0=Not Present
    p.SpecSlitLocation(:)  = fread( fid, 4, 'int16' ); % 162 Spectro Slit Location, 0=Not Present
    p.CustomTimingFlag     = fread( fid, 1, 'int16' ); % 170 T/F Custom Timing Used
    p.ExperimentTimeLocal  = fread( fid, 7, 'char' );%[TIMEMAX]172 Experiment Local Time as hhmmss\0
    p.ExperimentTimeLocal  = join(string(char(p.ExperimentTimeLocal)), '');
    p.ExperimentTimeUTC    = fread( fid, 7, 'char' );%[TIMEMAX] 179 Experiment UTC Time as hhmmss\0
    p.ExperimentTimeUTC    = join(string(char(p.ExperimentTimeUTC)), '');
    p.ExposUnits           = fread( fid, 1, 'int16' ); %186 User Units for Exposure
    p.ADCoffset            = fread( fid, 1, 'uint16' ); %188 ADC offset
    p.ADCrate              = fread( fid, 1, 'uint16' ); %190 ADC rate
    p.ADCtype              = fread( fid, 1, 'uint16' ); %192 ADC type
    p.ADCresolution        = fread( fid, 1, 'uint16' ); %194 ADC resolution
    p.ADCbitAdjust         = fread( fid, 1, 'uint16' ); %196 ADC bit adjust
    p.gain                 = fread( fid, 1, 'uint16' ); %198 gain
    p.Comments             = fread( fid, 80, 'char' ); %200 File Comments
    p.Comments             = join(string(char(p.Comments)), '');
    fseek( fid, 600, 'bof' );
    p.geometric            = fread( fid, 1, 'uint16' ); % 600 geometric ops: 
    % rotate 0x01,
    % reverse 0x02,
    % flip 0x04
    p.xlabel               = fread( fid, 16, 'char' ); % 602 intensity display string
    p.xlabel               = join(string(char(p.xlabel)), '');
    p.cleans               = fread( fid, 1, 'uint16' ); %618 cleans
    p.NumSkpPerCln         = fread( fid, 1, 'uint16' ); %620 number of skips per clean.
    p.SpecMirrorPos(:)     = fread( fid, 2, 'int16' ); %622 Spectrograph Mirror Positions
    p.SpecSlitPos(:)       = fread( fid, 4, 'single' ); %626 Spectrograph Slit Positions
    p.AutoCleansActive     = fread( fid, 1, 'int16' ); %642 T/F
    p.UseContCleansInst    = fread( fid, 1, 'int16' ); %644 T/F
    p.AbsorbStripNum       = fread( fid, 1, 'int16' ); %646 Absorbance Strip Number
    p.SpecSlitPosUnits     = fread( fid, 1, 'int16' ); %648 Spectrograph Slit Position Units
    p.SpecGrooves          = fread( fid, 1, 'single' ); %50 Spectrograph Grating Grooves
    p.srccmp               = fread( fid, 1, 'int16' ); %654 number of source comp.diodes
    p.ydim                 = fread( fid, 1, 'uint16' ); %656 y dimension of raw data.
    p.scramble             = fread( fid, 1, 'int16' ); %658 0=scrambled,1=unscrambled
    p.ContinuousCleansFlag = fread( fid, 1, 'int16' ); %660 T/F Continuous Cleans Timing Option
    p.ExternalTriggerFlag  = fread( fid, 1, 'int16' ); %662 T/F External Trigger Timing Option
    p.lnoscan              = fread( fid, 1, 'int32' ); %664 Number of scans (Early WinX)
    p.lavgexp              = fread( fid, 1, 'int32' ); %668 Number of Accumulations
    p.ReadoutTime          = fread( fid, 1, 'single' ); %672 Experiment readout time
    p.TriggeredModeFlag    = fread( fid, 1, 'int16' ); %676 T/F Triggered Timing Option
    fread( fid, 10, 'char' );    %skipping spare chars
    p.sw_version           = fread( fid, 16, 'char' ); %[FILEVERMAX] 688 Version of SW creating this file
    p.sw_version           = join(string(char(p.sw_version)), '');    
    temp                   = fread( fid, 1, 'int16' ); %704 
    switch temp
        case 1
            p.type = 'new120 (Type II)';
        case 2
            p.type = 'old120 (Type I)';
        case 3
            p.type = 'ST130';
        case 4
            p.type = 'ST121';
        case 5
            p.type = 'ST138';
        case 6
            p.type = 'DC131 (PentaMax)';
        case 7
            p.type = 'ST133 (MicroMax/SpectroMax)';
        case 8
            p.type = 'ST135 (GPIB)';
        case 9
            p.type = 'VICCD';
        case 10
            p.type = 'ST116 (GPIB)';
        case 11
            p.type = 'OMA3 (GPIB)';
        case 12
            p.type = 'OMA4';
        otherwise
            p.type = ['Unknown: ' num2str(temp)];
    end
    p.flatFieldApplied      = fread( fid, 1, 'int16' ); %706 1 if flat field was applied.
    fread( fid, 16, 'char' );    %skipping spare chars
    p.kin_trig_mode         = fread( fid, 1, 'int16' ); %724 Kinetics Trigger Mode
    p.dlabel                = fread( fid, 16, 'char' ); %[LABELMAX] 726 Data label.
    p.dlabel                = join(string(char(p.dlabel)), '');    
    fread( fid, 436, 'char' );    %skipping spare chars Spare_4[436] 742
    % HDRNAMEMAX = 120
    fseek( fid, 1178, 'bof' );
    p.PulseFileName         = fread( fid, 120, 'char' );%[HDRNAMEMAX] 1178 Name of Pulser File with Pulse Widths/Delays (for Z-Slice)
    p.PulseFileName         = join(string(char(p.PulseFileName)), '');
    p.AbsorbFileName        = fread( fid, 120, 'char' );%[HDRNAMEMAX] 1298 Name of Absorbance File (if File Mode)
    p.AbsorbFileName        = join(string(char(p.AbsorbFileName)), '');
    p.NumExpRepeats         = fread( fid, 1, 'uint32' ); %1418 Number of Times experiment repeated
    p.NumExpAccums          = fread( fid, 1, 'uint32' ); %1422 Number of Time experiment accumulated
    p.YT_Flag               = fread( fid, 1, 'int16' ); %1426 Set to 1 if this file contains YT data
    p.clkspd_us             = fread( fid, 1, 'single' ); %1428 Vert Clock Speed in micro-sec
    p.HWaccumFlag           = fread( fid, 1, 'int16' ); %1432 set to 1 if accum done by Hardware.
    p.StoreSync             = fread( fid, 1, 'int16' ); %1434 set to 1 if store sync used
    p.BlemishApplied        = fread( fid, 1, 'int16' ); %1436 set to 1 if blemish removal applied
    p.CosmicApplied         = fread( fid, 1, 'int16' ); %1438 set to 1 if cosmic ray removal applied
    p.CosmicType            = fread( fid, 1, 'int16' ); %1440 if cosmic ray applied, this is type
    p.CosmicThreshold       = fread( fid, 1, 'single' ); %1442 Threshold of cosmic ray removal.
    p.NumFrames             = fread( fid, 1, 'int32' ); %1446 number of frames in file.
    p.MaxIntensity          = fread( fid, 1, 'single' ); % MaxIntensity 1450 max intensity of data (future)
    p.MinIntensity          = fread( fid, 1, 'single' ); % MinIntensity 1454 min intensity of data future)
    p.ylabel                = fread( fid, 16, 'char' ); %[LABELMAX] 1458 y axis label.
    p.ylabel                = join(string(char(p.ylabel)), '');
    p.ShutterType           = fread( fid, 1, 'uint16' ); %1474 shutter type.
    p.shutterComp           = fread( fid, 1, 'single' ); %1476 shutter compensation time.
    p.readoutMode           = fread( fid, 1, 'uint16' ); %1480 readout mode, full,kinetics, etc
    p.WindowSize            = fread( fid, 1, 'uint16' ); %1482 window size for kinetics only.
    p.clkspd                = fread( fid, 1, 'uint16' ); %1484 clock speed for kinetics & frame transfer
    p.interface_type        = fread( fid, 1, 'uint16' ); %1486 computer interface (isa-taxi, pci, eisa,etc.)
    p.NumROIsInExperiment   = fread( fid, 1, 'int16' ); %1488 May be more than the 10 allowed in this header (if 0, assume 1)
    if p.NumROIsInExperiment == 0
        p.NumROIsInExperiment = 1;
    end
    fread( fid, 16, 'char' ); %Spare_5[16] 1490
    p.controllerNum         = fread( fid, 1, 'uint16' ); % 1506 if multiple controller system will have controller number data came from. This is a future item.
    p.SWmade                = fread( fid, 1, 'uint16' ); %1508 Which software package created this file
    p.NumROI                = fread( fid, 1, 'int16' ); %1510 number of ROIs used. if 0 assume 1.
    n = p.NumROI;
    if p.NumROI == 0
        n = 1;
        p.NumROI = 1;
    end
    if p.NumROI > 10
        n = 10;
    end
    % Make ROI struct
    ROIOffsets = [ 1512 1524 1536 1548 1560 1572 1584 1596 1608 1620 ];
    for i = 1 : n
        fseek( fid, ROIOffsets(i), 'bof' );
        p.ROI{i}.startx     = fread( fid, 1, 'uint16' ); % left x start value.
        p.ROI{i}.endx       = fread( fid, 1, 'uint16' ); %right x value.
        p.ROI{i}.groupx     = fread( fid, 1, 'uint16' ); %amount x is binned/grouped in hw.
        p.ROI{i}.starty     = fread( fid, 1, 'uint16' ); %top y start value.
        p.ROI{i}.endy       = fread( fid, 1, 'uint16' ); %bottom y value.
        p.ROI{i}.groupy     = fread( fid, 1, 'uint16' ); %amount y is binned/grouped in hw.
    end
    fseek( fid, 1632, 'bof' );
    p.FlatField             = fread( fid, 120, 'char' );%[HDRNAMEMAX] 1632 Flat field file name.
    p.FlatField             = join(string(char(p.FlatField)), '');
    p.background            = fread( fid, 120, 'char' );%[HDRNAMEMAX] 1752 background sub. file name.
    p.background            = join(string(char(p.background)), '');
    p.blemish               = fread( fid, 120, 'char' );%[HDRNAMEMAX] 1872 blemish file name.
    p.blemish               = join(string(char(p.blemish)), '');
    p.file_header_ver       = fread( fid, 1, 'single' ); %1992 version of this file header
    p.YT_Info               = fread( fid, 1000, 'char' );%[1000] 1996-2995 Reserved for YT information
    p.YT_Info               = join(string(char(p.YT_Info)), '');
    p.WinView_id            = fread( fid, 1, 'int32' );%2996 == 0x01234567L if file created by WinX
    
    % Pick some key parameters
    params.ROI = p.ROI;
    params.xdim = p.xdim;
    params.ydim = p.ydim;
    params.xlabel = p.xlabel;
    params.ylabel = p.ylabel;
    params.dlabel = p.dlabel;
    params.SpecGrooves = p.SpecGrooves;
    params.ExperimentTimeLocal = p.ExperimentTimeLocal;
    params.date = p.date; 
    params.exp_sec = p.exp_sec;
    
    % Load calibrations
    % X Calib
    fseek( fid, 3000, 'bof' );
    xc.offset           = fread( fid, 1, 'double' ); %3000 offset for absolute data scaling
    xc.factor           = fread( fid, 1, 'double' ); %3008 factor for absolute data scaling
    xc.current_unit     = fread( fid, 1, 'char' ); % 3016 selected scaling unit
    fseek( fid, 3018, 'bof' ); % char reserved1 3017 reserved
    temp                = fread( fid, 40, 'char' );% special string for scaling
    xc.special_string   = join(string(char(temp)), '');
    % char reserved2[40] 3058 reserved
    fseek( fid, 3098, 'bof' );
    xc.calib_valid      = fread( fid, 1, 'char' ); %3098 flag if calibration is valid
    xc.input_unit       = fread( fid, 1, 'char' ); %3099 current input units for "calib_value"
    xc.polynom_unit     = fread( fid, 1, 'char' ); %3100 linear UNIT and used in the "polynom_coeff"
    xc.polynom_order    = fread( fid, 1, 'char' ); %3101 ORDER of calibration POLYNOM
    xc.calib_count      = fread( fid, 1, 'char' ); %3102 valid calibration data pairs
    xc.pixel_position(:)= fread( fid, 10, 'double' ); % 3103 pixel pos. of calibration data
    xc.calib_value(:)   = fread( fid, 10, 'double' ); % 3183 calibration VALUE at above pos
    xc.polynom_coeff(:) = fread( fid, 6, 'double' ); % 3263 polynom COEFFICIENTS
    xc.laser_position   = fread( fid, 1, 'double' ); %3311 laser wavenumber for relative WN
    % char reserved3 3319 reserved
    fseek( fid, 3320, 'bof' );
    xc.new_calib_flag   = fread( fid, 1, 'uint8' ); % 3320 If set to 200, valid label below
    xc.calib_label(:)   = fread( fid, 81, 'char' ); % 3321 Calibration label (NULL term'd)
    xc.calib_label      = join(string(char(xc.calib_label)), '');
    xc.expansion(:)     = fread( fid, 87, 'char' ); %3402 Calibration Expansion area
    xc.expansion        = join(string(char(xc.expansion)), '');
    params.xcalib = xc;
   
    % Y Calib
    fseek( fid, 3498, 'bof' );
    yc.offset           = fread( fid, 1, 'double' ); %3489 offset for absolute data scaling
    yc.factor           = fread( fid, 1, 'double' ); %3497 factor for absolute data scaling
    yc.current_unit     = fread( fid, 1, 'char' ); % 3505 selected scaling unit
    fseek( fid, 3507, 'bof' ); % char reserved1 3506 reserved
    temp                = fread( fid, 40, 'char' );% 3507 special string for scaling
    yc.special_string   = join(string(char(temp)), '');
    fseek( fid, 3587, 'bof' ); %char reserved2[40] 3547 reserved
    yc.calib_valid      = fread( fid, 1, 'char' ); %3587 flag if calibration is valid
    yc.input_unit       = fread( fid, 1, 'char' ); %3588 current input units for "calib_value"
    yc.polynom_unit     = fread( fid, 1, 'char' ); %3589 linear UNIT and used in the "polynom_coeff"
    yc.polynom_order    = fread( fid, 1, 'char' ); %3590 ORDER of calibration POLYNOM
    yc.calib_count      = fread( fid, 1, 'char' ); %3591 valid calibration data pairs
    yc.pixel_position(:)= fread( fid, 10, 'double' ); % 3592 pixel pos. of calibration data
    yc.calib_value(:)   = fread( fid, 10, 'double' ); % 3672 calibration VALUE at above pos
    yc.polynom_coeff(:) = fread( fid, 6, 'double' ); % 3752 polynom COEFFICIENTS
    yc.laser_position   = fread( fid, 1, 'double' ); %3800 laser wavenumber for relative WN
    fseek( fid, 3809, 'bof' ); %3808 reserved
    yc.new_calib_flag   = fread( fid, 1, 'uint8' ); %3809 If set to 200, valid label below
    yc.calib_label(:)   = fread( fid, 81, 'char' ); %3812 Calibration label (NULL term'd)
    yc.calib_label      = join(string(char(yc.calib_label)), '');
    yc.expansion(:)     = fread( fid, 87, 'char' ); %3891 Calibration Expansion area
    yc.expansion        = join(string(char(yc.expansion)), '');
    params.ycalib = yc;
    
    % Wavelengths calibrations
    w = 1 : p.xdim;
    if xc.calib_valid == 1
        wavelengths = zeros( 1, p.xdim);
        for i = 1 : xc.polynom_order
            wavelengths = xc.polynom_coeff(i) * w.^(i-1) + wavelengths;
        end
        if yc.calib_valid == 1
            w = 1 : p.ydim;
            wavelengths = zeros( 1, p.ydim);
            for i = 1 : yc.polynom_order
                wavelengths = yc.polynom_coeff(i) * w.^(i-1) + wavelengths;
            end
            params.wavelengthsY = wavelengthsY;
        end
    else
        if yc.calib_valid == 1
            warning( ['Applying y-axis calibrations to wavelengths, no x-axis calibration found in ' filename ]);
            w = 1 : p.ydim;
            wavelengths = zeros( 1, p.ydim);
            for i = 1 : yc.polynom_order
                wavelengths = yc.polynom_coeff(i) * w.^(i-1) + wavelengths;
            end
        else
            warning( ['No valid x- or y-axis wavelengths calibrations found in ' filename...
                '. Reverting to pixels.' ]);
            wavelengths = w;
        end
    end

    % End of Calibration Structures (3978-4098)
    fseek( fid, 3978, 'bof' );
    p.Istring(:)        = fread( fid, 40, 'char' ); %3978 special intensity scaling string
    p.Istring           = join(string(char(p.Istring)), '');
    fseek(fid, 4043, 'bof' );% char Spare_6[25] 4018
    p.SpecType          = fread( fid, 1, 'uint8' ); %4043 spectrometer type (acton, spex, etc.)
    p.SpecModel         = fread( fid, 1, 'uint8' ); %4044 spectrometer model (type dependent)
    p.PulseBurstUsed    = fread( fid, 1, 'uint8' ); %4045 pulser burst mode on/off
    p.PulseBurstCount   = fread( fid, 1, 'uint32' ); %4046 pulser triggers per burst
    p.PulseBurstPeriod  = fread( fid, 1, 'double' ); %4050 pulser burst period (in usec)
    p.PulseBracketUsed  = fread( fid, 1, 'uint8' ); %4058 pulser bracket pulsing on/off
    p.PulseBracketType  = fread( fid, 1, 'uint8' ); %4059 pulser bracket pulsing type
    p.PulseTimeConstFast = fread( fid, 1, 'double' ); %4060 pulser slow exponential time constant (in usec)
    p.PulseAmplitudeFast = fread( fid, 1, 'double' ); %4068 pulser fast exponential amplitude constant
    p.PulseTimeConstSlow = fread( fid, 1, 'double' ); % 4076 pulser slow exponential time constant (in usec)
    p.PulseAmplitudeSlow = fread( fid, 1, 'double' ); %4084 pulser slow exponential amplitude constant
    p.AnalogGain        = fread( fid, 1, 'int16' ); % 4092 analog gain
    p.AvGainUsed        = fread( fid, 1, 'int16' ); %4094 avalanche gain was used
    p.AvGain            = fread( fid, 1, 'int16' ); %4096 avalanche gain value
    p.lastvalue         = fread( fid, 1, 'int16' ); %4098 Always the LAST value in the header
    params.full = p;
    
    % Load data
    fseek( fid, 4100, 'bof' );
    if p.NumFrames == 1 || p.NumFrames == 0
        data = fread( fid, [p.xdim, p.ydim], p.datatype );
    else
        for i = 1 : p.NumFrames
            data( :, :, i ) = fread( fid, [p.xdim, p.ydim], p.datatype );
        end
    end
else
    fclose( fid );
    error( ['Unsupported version of SPE file: ' num2str(params.version) '!']);
end
fclose(fid);

%% Subfunctions to parse Matlabs XML DOM to a structure
% ----- Subfunction parseChildNodes -----
    function [children,ptext,textflag] = parseChildNodes(theNode)
        % Recurse over node children.
        children = struct;
        ptext = struct; textflag = 'Text';
        if hasChildNodes(theNode)
            childNodes = getChildNodes(theNode);
            numChildNodes = getLength(childNodes);
            
            for count = 1:numChildNodes
                theChild = item(childNodes,count-1);
                [text,name,attr,childs,textflag] = getNodeData(theChild);
                
                if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
                    %XML allows the same elements to be defined multiple times,
                    %put each in a different cell
                    if (isfield(children,name))
                        if (~iscell(children.(name)))
                            %put existsing element into cell format
                            children.(name) = {children.(name)};
                        end
                        index = length(children.(name))+1;
                        %add new element
                        children.(name){index} = childs;
                        if(~isempty(fieldnames(text)))
                            children.(name){index} = text;
                        end
                        if(~isempty(attr))
                            children.(name){index}.('Attributes') = attr;
                        end
                    else
                        %add previously unknown (new) element to the structure
                        children.(name) = childs;
                        if(~isempty(text) && ~isempty(fieldnames(text)))
                            children.(name) = text;
                        end
                        if(~isempty(attr))
                            children.(name).('Attributes') = attr;
                        end
                    end
                else
                    ptextflag = 'Text';
                    if (strcmp(name, '#cdata_dash_section'))
                        ptextflag = 'CDATA';
                    elseif (strcmp(name, '#comment'))
                        ptextflag = 'Comment';
                    end
                    
                    %this is the text in an element (i.e., the parentNode)
                    if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                        if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                            ptext.(ptextflag) = text.(textflag);
                        else
                            %what to do when element data is as follows:
                            %<element>Text <!--Comment--> More text</element>
                            
                            %put the text in different cells:
                            % if (~iscell(ptext)) ptext = {ptext}; end
                            % ptext{length(ptext)+1} = text;
                            
                            %just append the text
                            ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                        end
                    end
                end
                
            end
        end
    end

% ----- Subfunction getNodeData -----
    function [text,name,attr,childs,textflag] = getNodeData(theNode)
        % Create structure of node info.
        
        %make sure name is allowed as structure name
        name = toCharArray(getNodeName(theNode))';
        name = strrep(name, '-', '_dash_');
        name = strrep(name, ':', '_colon_');
        name = strrep(name, '.', '_dot_');
        
        attr = parseAttributes(theNode);
        if (isempty(fieldnames(attr)))
            attr = [];
        end
        
        %parse child nodes
        [childs,text,textflag] = parseChildNodes(theNode);
        
        if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
            %get the data of any childless nodes
            % faster than if any(strcmp(methods(theNode), 'getData'))
            % no need to try-catch (?)
            % faster than text = char(getData(theNode));
            text.(textflag) = toCharArray(getTextContent(theNode))';
        end
        
    end

% ----- Subfunction parseAttributes -----
    function attributes = parseAttributes(theNode)
        % Create attributes structure.
        
        attributes = struct;
        if hasAttributes(theNode)
            theAttributes = getAttributes(theNode);
            numAttributes = getLength(theAttributes);
            
            for count = 1:numAttributes
                %attrib = item(theAttributes,count-1);
                %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
                %attributes.(attr_name) = char(getValue(attrib));
                
                %Suggestion of Adrian Wanner
                str = toCharArray(toString(item(theAttributes,count-1)))';
                k = strfind(str,'=');
                attr_name = str(1:(k(1)-1));
                attr_name = strrep(attr_name, '-', '_dash_');
                attr_name = strrep(attr_name, ':', '_colon_');
                attr_name = strrep(attr_name, '.', '_dot_');
                attributes.(attr_name) = str((k(1)+2):(end-1));
            end
        end
    end
end