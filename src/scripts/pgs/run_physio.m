% prep physio data before doing peak detection
[id pulse_scan_start resp_scan_start] = prep_physio_data();

% figure out order sequence
order = 0;
f = dir('../fmri/cry/FinalCry_*.txt');
if ~isempty(f)
    m = regexp(char(f(1).name),'FinalCry_(\d)\-.*\.txt','tokens');
    if ~isempty(m)
        order = str2double(char(m{1}));
    end
end    
if order == 0
   order = input(['Please enter cry order sequence for ' id ' (1-3): ']); 
end
    

% run peak detection for all of the combinations
% Which interval, own=1, other=2, garbled=3, rest=4, questions=5? 
cry = {'own','other','garbled','rest','questions'};

% write header
fid = fopen([id '_cry_data.csv'],'wt');
fprintf(fid, 'id,cry,block,AvgHR,Hi,Lo\n');

% iterate over cry types
for i=1:4
	fprintf('Subject: %s, Cry Interval: %s, Cry Order: %d\n',id,char(cry(i)),order);
	[avgHR hi lo] = peak_detect_pulseox(id,i,order);
	for j=1:length(avgHR)	
		fprintf(fid,'%s,%s,%d,%d,%d,%d\n',id,char(cry(i)),j,avgHR(j),hi(j),lo(j));
	end;
end
fclose(fid);
