%
% This is a MATLAB rewrite of a generic E-Prime TXT file
% reader that creates a ePrime stucture as output
% rewrite of Perl eprime2csv file
% author: Eugene Tseytlin (University of Pittsburgh)
%

function  eprime = read_eprime(eprime_file)
% read in the input file
fid = fopen(eprime_file);
tline = fgets(fid);
while ischar(tline)
    tline = strtrim(tline);
    
    % look at cues for start and stop of frame and header
    hd = regexp(tline,'\*\*\* ([A-Za-z]+) (Start|End) \*\*\*','tokens');
    kv = regexp(tline,'([\w\.]+): (.*)','tokens');
    
    if ~isempty(hd) && strcmp(hd{1}{2},'Start')
        % start new block
        block = strtrim(hd{1}{1});
        if ~strcmp(block,'Header')
            try
                eprime.(block).(level).N = 1 + eprime.(block).(level).N;
            catch
                eprime.(block).(level).N = 1;
            end
        end
    elseif ~isempty(hd) && strcmp(hd{1}{2},'End') && ~strcmp(block,'Header')
        % add blank values to key/val that has
        % not been mentioned
        f = fieldnames(eprime.(block).(level));
        for i=1:length(f)
            x = length(eprime.(block).(level).(f{i}));
            y = eprime.(block).(level).N;
            if ~strcmp(f{i},'N') && x < y
                for j=x:y-1
                    eprime.(block).(level).(f{i}){end+1} = '';
                end
            end
        end
    elseif ~isempty(kv) && strcmp(kv{1}{1},'Level')
        % read Level: 1
        level = ['Level' strtrim(kv{1}{2})];
    elseif ~isempty(kv)
        % while in block, add to current map
        key = strrep(strtrim(kv{1}{1}),'.','_');
        val = strtrim(kv{1}{2});
        % do merge with previous values
        % write header just a tad different
        if strcmp(block,'Header')
            try
                eprime.(block).(key);
            catch
                eprime.(block).(key) = {};
            end
            eprime.(block).(key){end+1} = val;
        else
            try
                eprime.(block).(level).(key);
            catch
                % prepand a buffer if it started as blank
                eprime.(block).(level).(key) = {};
                if eprime.(block).(level).N > 1
                    for j=1:eprime.(block).(level).N-1
                        eprime.(block).(level).(key){end+1} = '';
                    end
                end
            end
            eprime.(block).(level).(key){end+1} = val;
        end
    end
    tline = fgets(fid);
end
fclose(fid);
