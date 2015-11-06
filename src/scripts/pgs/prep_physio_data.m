function [id pulse_scan_start resp_scan_start] = prep_physio_data()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This program pulls out and computes the onset time of the physio data from the physio files
%  These files should start with 'wpc' and end with the extensions .puls or .resp 
%  This program also takes the data from the first dicom file from the cry data and computes the
%  onset time of the cry scan from it. THIS DATA IS CURRENTLY INPUT FROM THE USER 
%
%  The files 'clear_stats.sh' must be located in the same directory as this program in order
%  for it to run properly.
%
%  Written by Rebecca McNamee
%  January 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND RENAME INPUT
%
% id = five digit subject id
% pulse_onset_time.txt and resp_onset_time.txt are computed from the 'clear_stats.sh' program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
id = '';
m = regexp(pwd(),'.*/\d+\.(\d+)/physio','tokens');
if ~isempty(m)
	id = char(m{1});
else
	m = regexp(pwd(),'.*/(\d+)','tokens');
	if ~isempty(m)
		id = char(m{1});
	end 
end
% prompt if can't figure out
if isempty(id)
	id = input('Please enter the subject ID: ','s');
end

%this should be taken care by calling script
%in = sprintf('! source clear_stats.sh %s',id); 
%eval(in)
%
load pulse_onset_time.txt
load resp_onset_time.txt
%
tpul=pulse_onset_time/1000;
tresp=resp_onset_time/1000;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER ENTERS THE HOURS, MINUTES, and SECONDS from the dicom file
% Data is converted to seconds since midnight, as this is the format of the physio data
%
% hsm = hours since midnight
% msm = minutes since midnight
% ssm = seconds since midnight
% tssm = total seconds since midnight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
mr = dir('MR.*');
hr='';
min='';
sec='';

if ~isempty(mr)
	mr = char(mr(1).name);
	m = regexp(mr,'MR.*\.(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})(\d)\d+','tokens');
	if ~isempty(m)
		year=char(m{1}{1});
		mo=char(m{1}{2});
		day=char(m{1}{3});
		hr=char(m{1}{4});
		min=char(m{1}{5});
		sec=char(m{1}{6});
		ml=char(m{1}{7});
		if str2num(ml) > 5
			sec = num2str(str2num(sec)+1);
		end
	end
end

if isempty(hr)
	hr=input('Please enter the hours from the first cry dicom file: ','s')
	min=input('Please enter the minutes from the first cry dicom file: ','s')
	sec=input('Please enter the seconds from the first cry dicom file: ','s')
end 

%
hr1=str2num(hr);
min1=str2num(min);
ssm=str2num(sec);
hsm=hr1*60*60;
msm=min1*60;
tssm=hsm+msm+ssm;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the difference between the scan start time and the physio start time 
% (physio should start first)
%
% dsm_pulse = difference since midnight between scan onset and pulse-ox onset
% dsm_resp = difference since midnight between the scan onset and repiration onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
dsm_pulse=tssm-tpul;
dsm_resp=tssm-tresp;
% 
pulse_scan_start=dsm_pulse*50;
resp_scan_start=dsm_resp.*50;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the output files to use in the peak_detect programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
save pulse_scan_start.txt pulse_scan_start -ascii;
save resp_scan_start.txt resp_scan_start -ascii;
