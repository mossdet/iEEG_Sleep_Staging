function [keys, sections, subsections] = laydat_get_info(fileName)

    [keys,sections,subsections] = readallkeys(fileName);
    
    return    


function [keys,sections,subsections] = readallkeys(fileName)
% Reads all the keys out as well as the sections and subsections

keys = [];
sections = [];
subsections = [];
% Read the whole file's contents out
try
    dataout = textread(fileName,'%s','delimiter','\n');
catch
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
nLines = size(dataout,1);

% Go through all the lines and construct the keys variable
keys = cell(nLines,4);
sections = cell(nLines,1);
subsections = cell(nLines,2);
keyN = 0;
secN = 0;
subsecN = 0;
secStr = '';
subsecStr = '';
for ii=1:nLines
    [status,value,key] = processiniline(dataout{ii});
    if status == 1
        secN = secN + 1;
        secStr = value;
        sections(secN) = {secStr};
    elseif status == 2
        subsecN = subsecN + 1;
        subsecStr = value;
        subsections(subsecN,:) = {secStr,subsecStr};
    elseif status == 3
        keyN = keyN + 1;
        keys(keyN,:) = {secStr,subsecStr,key,value};
    end
end
keys(keyN+1:end,:) = [];
sections(secN+1:end,:) = [];
subsections(subsecN+1:end,:) = [];
%------------------------------------



%------------------------------------
function [status,value,key] = processiniline(line)
% Processes a line read from the ini file and
% returns the following values:
%   - status:  -1   => unknown string found
%               0   => empty line found
%               1   => section found
%               2   => subsection found
%               3   => key-value pair found
%               4   => comment line found (starting with ;)
%   - value:    value-string of a key, section, subsection, comment, or unknown string
%   - key:      key as string

status = 0;
value = [];
key = [];
line = strim(line);                         % removes any leading and trailing spaces
if isempty(line)                            % empty line
    return
end
if strcmpi(line(1),';')                     % comment found
    status = 4;
    value = line(2:end);
elseif (line(1) == '[') & (line(end) == ']') & (length(line) >= 3)  % section found
    value = lower(line(2:end-1));
    status = 1;
elseif (line(1) == '{') &...                % subsection found
       (line(end) == '}') & (length(line) >= 3)
    value = lower(line(2:end-1));
    status = 2;
else                                        % either key-value pair or unknown string
    pos = findstr(line,'=');
    if ~isempty(pos)                        % key-value pair found
        status = 3;
        key = lower(line(1:pos-1));
        value = line(pos+1:end);
        key = strim(key);                   % removes any leading and trailing spaces
        value = strim(value);               % removes any leading and trailing spaces
        if isempty(key)                     % empty keys are not allowed
            status = 0;
            key = [];
            value = [];
        end
    else                                    % unknown string found
        status = -1;
        value = line;
    end
end


%------------------------------------
function outstr = strim(str)
% Removes leading and trailing spaces (spaces, tabs, endlines,...)
% from the str string.
if isnumeric(str);
    outstr = str;
    return
end
ind = find( ~isspace(str) );        % indices of the non-space characters in the str    
if isempty(ind)
    outstr = [];        
else
    outstr = str( ind(1):ind(end) );
end
