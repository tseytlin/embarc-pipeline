function [avgHR hi lo] = peak_detect_pulseox(id,loop2,order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This program lines up the fMRI data from the cry paradigm with the pulse-ox data, 
%  finds the peaks in the pulse ox data to compute the heart rate, and finds the
%  average of the heart rate over the each cry interval type (own, other, etc)
%
%  Run 'prep_physio_data.m' prior to using this program
%
%  Written by Rebecca McNamee
%  January 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% LOAD INPUT DATA AND START TIMES
%
% id = subject five digit id number (without the date)
% loop = cry interval type 
% new_pulseox, k = raw pulse-ox data 
% new_resp, r = raw respiratory data
% tpul = time (from midnight) that the pulse ox data was turned on
% tresp = time (from midnight) that the respiration data was turned on 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
format longg 
%id = input('Please enter the subject ID: ','s')
%
%loop= input('Which interval, own=1, other=2, garbled=3, rest=4, questions=5? ','s') 
%loop2=str2num(loop);
%
load new_pulseox.txt
load new_resp.txt
%
k=new_pulseox;
r=new_resp;
%
load pulse_onset_time.txt
load resp_onset_time.txt 
%
tpul=pulse_onset_time/1000;
tresp=resp_onset_time/1000;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE SCANNER PULSES AND RECORD PULSE ONSETS (spikes where the signal=5000);
%
% c, p = loop counters
% realk_pre = pulse-ox data with spikes removed
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
p=1;
c=1;
for j=1:length(k)
 if k(j)==5000;
  on(c)=j;
  c=c+1;
  j=j+1;
 elseif k(j)==5003;
  j=j+1;
 else
  realk_pre(p)=k(j);  
  p=p+1;
 end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Low-pass filter the data while preserving the correct phase (Butterworth filter)
%
% f1, f2 = filter parameters
% realk = realk_pre after filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
[f1 f2]=butter(3,0.1);
realk=filtfilt(f1,f2,realk_pre);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Load the start time of the scan (taken from the first dicom file from the cry scan data 
% and saved to the output file 'pulse_scan_start.txt' upon running prep_physio_data.m)
%
% a = rounded start time of the scan
% cry_length_time = length of cry paradigm in time
% cry_length_samples = length of cry paradigm in samples (based on 50Hz sampling rate)
% b = rounded end time of the scan
% realk_cry = realk truncated to match the time of the cry scan
% t = time of each physio data point computed from the 50Hz sampling rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
load pulse_scan_start.txt
a=round(pulse_scan_start);
cry_length_time=600; %% in seconds
cry_length_samples=600*50;
b=round(pulse_scan_start+cry_length_samples);
realk_cry=realk(a:b);
t=0:0.02:(length(realk_cry)/50); %%% Time in seconds;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Compute the physio stats based on the segment of cry data of interest (own baby cry,
% other baby cry, garbled, etc) using a switch loop
%
% onset = onset of each cry data interval
% dur = duration of the interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if order == 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
	% Cry order times - ORDER 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%

	switch loop2
	 case 1
	 onset=50.*[25 259 301 535];
	 dur=36*50;
	 case 2
	 onset=50.*[109 193 409 493];
	 dur = 36*50;
	 case 3
	 onset=50.*[67 151 343 451];
	 dur=36*50;
	 case 4
	 onset=50.*[1 235 385 577];
	 dur=36*50;
	 case 5
	 onset=50.*[61 103 145 187]; % 229 295 337 379 445 529 571]; % Can only evaulate 4 intervals
	 dur=6*50;
	end
elseif order == 2
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
	% Cry order times - ORDER 2
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	switch loop2
	 case 1
	 onset=50.*[151 259 409 451];
	 dur=36*50;
	 case 2
	 onset=50.*[25 193 343 535];
	 dur = 36*50;
	 case 3
	 onset=50.*[67 109 301 493];
	 dur=36*50;
	 case 4
	 onset=50.*[1 235 385 577];
	 dur=36*50;
	 case 5
	 onset=50.*[61 103 145 187]; % 229 295 337 379 445 529 571];
	 dur=6*50;
	end
else
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
	% Cry order times - ORDER 3
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	switch loop2
	 case 1
	 onset=50.*[67 301 493 535];
	 dur=36*50;
	 case 2
	 onset=50.*[109 151 343 451];
	 dur = 36*50;
	 case 3
	 onset=50.*[25 193 259 409];
	 dur=36*50;
	 case 4
	 onset=50.*[1 235 385 577];
	 dur=36*50;
	 case 5
	 onset=50.*[61 103 145 187]; % 229 295 337 379 445 529 571]; % Can only evaulate 4 intervals
	 dur=6*50;
	end
end 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Find the pulse-ox peaks to compute the heart rate using 5 criteria
%
% After choosing the above case, the program will loop through each of the 4 intervals
%
% tick = counter for each of the onsets intervals (case dependent using the above switch) 
% cryint = onset time 
% fix = duration of interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
tick=1;
clf;
for cryint=onset;
 fix=cryint+dur;
 realk_cry=realk(cryint:fix);
%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CRITERIA 1: Peak must occur at an inflection point
 % 
 % giff = first derivative of realk_cry
 % lol = first marker of potential peaks
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear lol;
 clear gdiff;
 gdiff=diff(realk_cry);
 p=1;
 for j=1:length(gdiff)-1
 if gdiff(j)>=0 & gdiff(j+1)<0
  %tmpk(p)=realk_cry(j+1); %%% Value of pulse-ox at pulse
  lol(p)=j+1;	    %%% Value of indice at pulse
  p=p+1;
  end
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CRITERIA 2: Majority of points around peak must be increasing and decreasing consistently
 %             This again uses the first derivative data
 %
 % tmp =  holds the lol point for evaulation
 % p = counter for peaks of lol as they are evaluated
 % d, d2 = evalutaion criterias (if 3 out of 4 points are consistently decreasing or increasing, 
 %         the point is considered as a potential peak 
 % lol2 = second marker of potential peaks (passed criteria 2)
 % 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear q;
 clear lol2;
 clear tmp;
 p=1;
 tmp=lol(p);
 while tmp<=4
  q(p)=lol(p);
  p=p+1;
  tmp=lol(p);
 end
 while p<length(lol)-4;
  d=[gdiff(tmp-4)>=0 gdiff(tmp-3)>=0 gdiff(tmp-2)>=0 gdiff(tmp-1)>=0];
  d2=[gdiff(tmp+1)<=0 gdiff(tmp+2)<=0 gdiff(tmp+3)<=0 gdiff(tmp+4)<=0];
  if sum(d)>=3 & sum(d2)>=3 
   q(p)=tmp;  %%%% Value of indice at pulse
  end
  p=p+1;
  tmp=lol(p);
 end
 while p<=length(lol);
  q(p)=lol(p);
  p=p+1;
 end
 lol2=nonzeros(q)';
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CRITERIA 3: Amplitude of peak must be greater than the average of the surrounding signal
 %             (Plus and minus 10 points)
 %
 % w = holds the point for evaluation
 % p = counter for peaks of lol2 as they are evaluated
 % lol3 = third marker of potential peaks (passed criteria 3) 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 clear w;
 clear lol3;
 clear val;
 p=1;
 val=lol2(p);
 while val<=10
  w(p)=lol2(p);
  p=p+1;
  val=lol2(p);
 end
 while val<max(lol2)-10
  if realk_cry(val)>mean(realk_cry(val-10:val)) & realk_cry(val)>mean(realk_cry(val:val+10)) 
   w(p)=val;
  end
  p=p+1;
  val=lol2(p);
 end
 while p<=length(lol2);
  w(p)=lol2(p);
  p=p+1;
 end
 lol3=(nonzeros(w))'; 
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CRITERIA 4: TIME CRITERIA (max heart rate can only be 240 beats per min)
 %
 % c = interval between one peak to the next (in samples)
 % x = holds point for evaluation
 % p = counter for peaks of lol3 as they are evaluated
 % lol4 = fourth marker of potential peaks (passed criteria 4)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 p=1;
 clear lol4;
 clear x;
 while p<length(lol3)-1
  c=lol3(p+1)-lol3(p);
  if c>=12.5;
   x(p)=lol3(p);
  end
  p=p+1;
 end
 while p<=length(lol3)
  x(p)=lol3(p);
  p=p+1;
 end
 lol4=nonzeros(x)';
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CRITERIA 5: Amplitude of the peak must again be greater than average of the surrounding signal
 %             (Plus and minus 30 points)
 % val2 = holds point for evaluation
 % p = counter for peaks of lol4 as they are evaluated
 % finaltmp = detected peak (passed criteria 5)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 p=1;
 clear val2;
 clear finaltmp;
 val2=lol4(p);
 while val2<=30
  finaltmp(p)=lol4(p);
  p=p+1;
  val2=lol4(p);
 end
 while val2< max(lol4)-30
  if realk_cry(val2)>mean(realk_cry((val2-30):(val2+30)))
   finaltmp(p)=lol4(p);
  end
  p=p+1;
  val2=lol4(p);
 end
 while p<=length(lol4)
  finaltmp(p)=lol4(p);
  p=p+1;
 end
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Compute final statistics and plot the signal showing detected peaks
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 h=nonzeros(finaltmp); %%%%% Indice of detected heart beat
 hbd=0.02.*h;          %%%%% Time of detected heart beat
 peak=realk_cry(h);    %%%%% Amplitude of detected heart beat
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 subplot(4,1,tick)
 plot(t(1:length(realk_cry)),realk_cry(1:length(realk_cry)))
 hold on
 plot(hbd(1:length(hbd)),peak(1:length(hbd)),'*')
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Compute the RR intervals 
 %
 % Point_int = number of sample points between detected peaks
 % RR_int = time between detected peaks
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 p=1;
 clear RRint;
 while p<length(h-1)
  Point_int(p)=h(p+1)-h(p);
  p=p+1;
 end
 p=1;
 while p<length(hbd-1)
  RRint(p)=hbd(p+1)-hbd(p);
  p=p+1;
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Replace outliers with average of surrounding intervals
 %
 % If the interval is found to be greater than 1.5 times or less than 0.5 times the surrounding 
 % intervals,replace with the average of the surrounding signal
 %
 % out_high = counter of number of outlying points that were greater than 1.5 times the previous interval
 % out_low = counter of number of outlying points that were less than 0.5 times the previous interval
 % v = holds surrounding interval values for evaluation
 % RRint_final = final peak to peak values with the outliers removed
 % RRint_vec = final peak to peak values over the specific cry interval 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 p=1;
 clear RRint_final;
 while p<=2
  RRint_final(p)=RRint(p);
  p=p+1;
 end
 out_high=0;
 out_low=0;
 while p<(length(RRint)-2)
  if RRint(p)>1.5*mean(RRint((p-2):(p+2)))
   v=[RRint(p-2) RRint(p+2) RRint(p-1) RRint(p+1)];
   RRint_final(p)=mean(v);
   p=p+1;
   out_low=out_low+1;
  elseif RRint(p)<0.5*mean(RRint((p-2):(p+2)))
   v=[RRint(p-2) RRint(p+2) RRint(p-1) RRint(p+1)];
   RRint_final(p)=mean(v);
   p=p+1;
   out_high=out_high+1;
  else
   RRint_final(p)=RRint(p);
   p=p+1;
  end
 end
 %%%%%
 if tick==1
  RRint_vec1=RRint_final';
 elseif tick==2
  RRint_vec2=RRint_final';
 elseif tick==3
  RRint_vec3=RRint_final';
 elseif tick==4
  RRint_vec4=RRint_final';
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Compute the final statistics
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 avgRR(1,tick)=mean(RRint_final);             %%%%% Average peak to peak value over the interval
 avgHR(1,tick)=((1/avgRR(1,tick))*60);        %%%%% Heart rate over the interval
 hi(1,tick)=out_high/length(RRint_final)*100; %%%%% Percentage of replaced high outliers over the interval
 lo(1,tick)=out_low/length(RRint_final)*100;  %%%%% Percentage of replaced low outliers over the interval
 tick=tick+1;
end;
fprintf('The Average Heart Rate over this entire interval was: %f bpm\n',avgHR);
fprintf('The percentage of replaced outliers that were too high was: %f \n', hi);
fprintf('The percentage of replaced outliers that were too low was: %f \n', lo);
fprintf('\n'); 

