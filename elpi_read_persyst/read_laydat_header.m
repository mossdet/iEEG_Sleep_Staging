function ft_hdr = read_laydat_header(lay_file_name)
    % input - file name of a .lay file that corresponds to a lay-dat pair
    % output:   header - information from .lay file
    %           record - EEG data from .dat file
    
    
    
    %-------------------
    %     lay file
    %-------------------
    [lay_file_dir, ~, ~] = fileparts(lay_file_name);

    %use inifile to read .lay file as .ini file
    data = laydat_get_info(lay_file_name);
    %change "empty" cells to '...'
    emptyCells = cellfun(@isempty,data);
    data(find(emptyCells)) = cellstr('');
    
    %fileinfo
    %find the fileinfo section of .lay file
    fileInfoArray = data(find(strcmp('fileinfo',data(:,1))),:);
    for ii=1:size(fileInfoArray,1)
        %get field name 
        field_name = char(fileInfoArray(ii,3));
        %make it syntax-valid (remove spaces, etc)
        valid_field_name = matlab.lang.makeValidName(field_name);
        %add field to struct with corresponding string data
        rawhdr.fileinfo.(valid_field_name)=char(fileInfoArray(ii,4));
    end
    
    %patient
    patientArray = data(find(strcmp('patient',data(:,1))),:);
    for ii=1:size(patientArray,1)
        %same as above, but more compact
        rawhdr.patient.(matlab.lang.makeValidName(char(patientArray(ii,3))))=char(patientArray(ii,4));
    end
    
    %archive
    archiveArray = data(find(strcmp('archive',data(:,1))),:);
    for ii=1:size(archiveArray,1)
        %same as above, but more compact
        rawhdr.archive.(matlab.lang.makeValidName(char(archiveArray(ii,3))))=char(archiveArray(ii,4));
    end
    %assert(strcmp(rawhdr.archive.studyduration, rawhdr.archive.studyspan), "Study duration and study span are different")

    if strcmp(rawhdr.fileinfo.datatype, '7')
        precision_num_bits = 32;
    else
        precision_num_bits = 16;
    end
    s = dir(strcat(lay_file_dir, '\', rawhdr.fileinfo.file));
    filesize_bits = s.bytes*8;
    if isfield(rawhdr, 'archive')
        file_tot_bits = (str2double(rawhdr.archive.datafilesize)*8);
        assert(filesize_bits==file_tot_bits, "Could not determine number of samples")
    end
    nSamples = filesize_bits / precision_num_bits / str2double(rawhdr.fileinfo.waveformcount);
    
    %montage
    montageArray = data(find(strcmp('montage',data(:,1))),:);
    for ii=1:size(montageArray,1)
        %storing a 2D vector of info on specific montage, rather than a string
        montage_data = data(find(strcmp(montageArray(ii,3),data(:,1))),3:4);
        rawhdr.montage.(matlab.lang.makeValidName(char(montageArray(ii,3)))) = montage_data;
    end
    
    %sampletimes
    sampleTimesArray = data(find(strcmp('sampletimes',data(:,1))),:);
    %sampletimes is a cell array, rather than a struct
    rawhdr.sampletimes = {};
    for ii=1:size(sampleTimesArray,1)
        %store it as a string like in .lay file
        %sampletimes_data = strcat(char(sampleTimesArray(ii,3)),'=',char(sampleTimesArray(ii,4)));
        rawhdr.sampletimes{ii}.sample = str2double(char(sampleTimesArray(ii,3)));
        rawhdr.sampletimes{ii}.time = str2double(char(sampleTimesArray(ii,4)));
    end
    
    %channelmap
    channelMapArray = data(find(strcmp('channelmap',data(:,1))),:);
    %channelmap is a cell array of channel names
    rawhdr.channelmap = {};
    for ii=1:size(channelMapArray,1)
        rawhdr.channelmap{ii} = char(channelMapArray(ii,3));
    end
    
    %-------------------
    %      header
    %-------------------
    %move some info from raw header to header
    if isfield(rawhdr,'fileinfo')
        %checking individual fields exist before moving them
        if isfield(rawhdr.fileinfo,'file')
            header.datafile = rawhdr.fileinfo.file;
        end
        if isfield(rawhdr.fileinfo,'samplingrate')
            header.samplingrate = str2double(rawhdr.fileinfo.samplingrate);
        end
        if isfield(rawhdr.fileinfo,'waveformcount')
            header.waveformcount = str2double(rawhdr.fileinfo.waveformcount);
        end
    end
    %making start time into one form
    date = strrep(rawhdr.patient.testdate,'.','/');
    time = strrep(rawhdr.patient.testtime,'.',':');
    dn = datenum(strcat(date, ',', time));
    header.starttime = datetime(dn,'ConvertFrom','datenum');
    header.patient = rawhdr.patient;
    
    %-------------------
    %     comments
    %-------------------
    lay_file_ID = fopen(lay_file_name);
    %comments need to be extracted manually
    header.annotations = {};
    rawhdr.comments = {};
    comments = 0;
    cnum = 1;
    tline = fgets(lay_file_ID);
    while ischar(tline)
        if (1==comments)
            contents = strsplit(tline,',');
            if numel(contents) < 5
                %this means there are no more comments
                break;
            elseif numel(contents) > 5
                %rejoin comments that have a comma in the text
                contents(5) = {strjoin(contents(5:end),',')};
            end
            %raw header contains just the original lines
            rawhdr.comments{cnum} = tline;
            samplenum = str2double(char(contents(1))) * str2double(char(rawhdr.fileinfo.samplingrate));
            %this calculates sample time
            i=1;
            while i<numel(rawhdr.sampletimes) && samplenum > rawhdr.sampletimes{i+1}.sample
                i=i+1;
            end
            samplenum = samplenum - rawhdr.sampletimes{i}.sample;
            samplesec = samplenum / str2double(char(rawhdr.fileinfo.samplingrate));
            timesec = samplesec + rawhdr.sampletimes{i}.time;
            commenttime= datestr(timesec/86400, 'HH:MM:SS');
            %use date calculated earlier
            dn = datenum(strcat(date, ',', commenttime));
            %put all that into a struct in the header
            header.annotations{cnum}.time = datetime(dn,'ConvertFrom','datenum');
            header.annotations{cnum}.duration = str2double(char(contents(2)));
            header.annotations{cnum}.text = char(contents(5));
            cnum=cnum+1;
        elseif strncmp('[Comments]',tline,9)
            %read until get to comments
            comments = 1;
        end
        tline=fgets(lay_file_ID);
    end
    fclose(lay_file_ID);
    
    %put raw header in header
    header.rawheader = rawhdr;
    
    % Build header structure compatible with fieldtrip
    ft_hdr.Fs = header.samplingrate;
    ft_hdr.nChans = str2double(rawhdr.fileinfo.waveformcount);
    
    ft_hdr.label = header.rawheader.channelmap';
    ft_hdr.label = strrep(ft_hdr.label, '-ref', '');
    ft_hdr.label = upper(ft_hdr.label);
    %ft_hdr.nSamples = size(record,2);
    ft_hdr.nSamples = nSamples;
    ft_hdr.nSamplesPre = 0;
    ft_hdr.nTrials = 1;
    ft_hdr.orig = 0;
    ft_hdr.chantype = repmat(' ', size(ft_hdr.label,1), size(ft_hdr.label,2));
    ft_hdr.chanunit = repmat('uV', size(ft_hdr.label,1), size(ft_hdr.label,2));
    ft_hdr.physmin = repmat(-2^(precision_num_bits-1), size(ft_hdr.label,1), size(ft_hdr.label,2));
    ft_hdr.physmax = repmat(2^(precision_num_bits-1), size(ft_hdr.label,1), size(ft_hdr.label,2));
    ft_hdr.digimin = repmat(-2^(precision_num_bits-1), size(ft_hdr.label,1), size(ft_hdr.label,2));
    ft_hdr.digimax = repmat(2^(precision_num_bits-1), size(ft_hdr.label,1), size(ft_hdr.label,2));

    % Archived files can be an extract of another file and thus have a
    % separate studydate field that could be different to the testdate from
    % the original file
    if isfield(rawhdr, 'archive')
        study_start_datetime_str = strcat(rawhdr.archive.studydate, " 0:0:0");
    else
        study_start_datetime_str = strcat(rawhdr.patient.testdate, " 0:0:0");
    end

    study_start_datetime = datetime(datevec(study_start_datetime_str, 'yyyy.mm.dd HH:MM:SS')) +seconds(round(rawhdr.sampletimes{1, 1}.time));
    ft_hdr.study_starttime = datevec(study_start_datetime);

    for i = 1:length(rawhdr.sampletimes)
        datetime_info = datetime(datevec(study_start_datetime_str, 'yyyy.mm.dd HH:MM:SS'))+seconds(round(rawhdr.sampletimes{1, i}.time));
        rawhdr.sampletimes{1, i}.study_datetime = datevec(datetime_info);
    end
    ft_hdr.sampletimes = rawhdr.sampletimes;

    ft_hdr.calibration = str2double(rawhdr.fileinfo.calibration);
    ft_hdr.datatype = rawhdr.fileinfo.datatype;
    ft_hdr.annotations = rawhdr.comments;

end