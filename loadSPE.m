function [ data, wavelengths, params ] = loadSPE( filename )
%% LOADSPE load LightField SPE v3.0 or WinSpec SPE v2.x data from a file
% [data, wavelengths, params] = loadSPE( filename )
% If the file has only one ROI and one frame then data is a simple (2D) array.
% If the file has one ROI but several frames then data is 3D array of
% dimensions ( height, width, number of frames). If the file has one frame
% and several ROIs data is a struct data.ROI{numOfROIs}( height, width). If
% the file has several frames and several ROIs in each frame then data is
% an array of structs data(numOfFrames).ROI{numOfROIs}( height, width).
%
% params.version contains version info of the file. For SPE v2.x the only
% other parameter passed is params.exposureTime. For SPE v3.x the whole XML
% footer is parsed as a struct in params.SpeFormat. 
% SPE v3.0 file and XML specification is at <a href="matlab:web('ftp://ftp.princetoninstruments.com/public/Manuals/Princeton%20Instruments/SPE%203.0%20File%20Format%20Specification.pdf')">the PI website</a>.
%
% Examples:
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
% This code is based on load_winspec function for SPE v2.x files written by
% R. Bradley, modified by P. Walker, L. Tinkler.
%
% XML parsing is implemented using functions based on <a href="matlab:web('https://mathworks.com/matlabcentral/fileexchange/28518-xml2struct')">xml2struct</a>,
% written by W. Falkena, modified by A. Wanner, I. Smirnov, X. Mo.
% xml2struct (C) 2012, W. Falkena
%
% (C) 2017, M. Sich, The University of Sheffield
% v1.0 20-12-2017

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
    % SPE v2.x
    % All info is contained in the header.
    
    % load header
    fseek( fid, 0, 'bof' );
    [ header, ~ ] = fread( fid, 4100, 'uint8' );
    
    % determine data type
    dataType = header( 109 );
    dataSize = max( (dataType<2)+1 ) * 2;
    strDataType = { 'float32', 'int32', 'int16', 'uint16' };
    if dataType > 3
        fclose( fid );
        error( ['Unknown WinSpec data type: ' num2str(dataType)] );
    end
    
    % data size
    xdim = header(43) + header(44)*2^8;
    ydim = header(657) + header(658)*2^8;
    frames = header(1447) + header(1448)*2^8 +header(1449)*2^16 + header(1450)*2^24;
    
    %Work out number of 16MB reads needed to load all data (for some reason
    %fread only works for 16MB at a time...)
    size_per_col = dataSize*xdim;
    max_cols_per_fread = floor(16*1024*1024/size_per_col);
    num_reads = 1 + floor(ydim/max_cols_per_fread);
    
    %Allocate storage for data
    data = zeros( ydim, xdim, frames );
    myType=strcat('*',strDataType{dataType+1});
    for framenum = 1 : frames
        dataxy = zeros( xdim, ydim );
        for a=1:num_reads
            startY = 1 + max_cols_per_fread*(a-1);
            endY = startY + max_cols_per_fread - 1;
            if endY > ydim
                endY = ydim;
            end
            curr_y_size = endY-startY+1;
            dataTmp=zeros(xdim,curr_y_size);
            seekstatus = fseek(fid,4100+dataSize*xdim*(startY-1) +...
                dataSize*xdim*ydim*(framenum-1),'bof');
            [ dataTmp, ~ ]=fread(fid,[xdim curr_y_size],myType);
            dataxy(:,startY:endY) = dataTmp;
        end
        %Rotate xy and store in main data array
        dataxy = dataxy.';
        data(:,:,framenum) = dataxy;
    end
    
    % load exposure time
    fseek( fid, 10, 'bof');
    params.exposureTime = fread( fid, 1, 'single' );
    
    % load wavelengths
    order = header( 3102 );
    fseek( fid, 3263, 'bof' );
    [coeffs, ~ ] = fread( fid, order+1, 'double' );
    wavelengths = 1:xdim;
    powers = (1:length(coeffs))-1;
    wavelengths = sum( (coeffs*ones(1,length(wavelengths))).*...
        ((ones(length(coeffs),1)*wavelengths).^...
        ((ones(length(wavelengths),1)*powers).')));
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