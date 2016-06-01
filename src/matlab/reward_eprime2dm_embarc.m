function reward_eprime2dm_embsarc(eprime_file)
	% extract subject if possible
	directory = fileparts(eprime_file);
	participant = 'subject';
	remote_site = '';
	pt = regexp(eprime_file,'(.*)embarc_CU_([A-Z]{2}\d+)_\dR\d_mri_fmriraw_\d+/.*','tokens');
	if ~isempty(pt)
		maindir = char(pt{1}{1});
		participant = char(pt{1}{2});
		% figure out remote site
	    pt = regexp(maindir,'.*/([A-Z]{2})/','tokens');
	    if ~isempty(pt)
	    	remote_site = char(pt{1});
	    end
	end

	v1 = 1; %over 13
	v2 = 1; %over 71


	[st data]= system(['eprime2csv ', eprime_file, ' 1']);
	[head content] = parse_csv(data);
	[st2 data2]= system(['eprime2csv ', eprime_file, ' 2']);
	[head2 content2] = parse_csv(data2);
	[st3 data3]= system(['eprime2csv ', eprime_file, ' 3']);
	[head3 content3] = parse_csv(data3);
	[st4 data4]= system(['eprime2csv ', eprime_file, ' 4']);
	[head4 content4] = parse_csv(data4);

	anticipation_onset=0;
	outcome_trial = 0;
	outcome_onset=0;
	Value_R = 0;
	PE_R = 0;

	GamStim_RESP_c1 = 1;
	GamStim_RESP_c2 = 2;
	if strcmp(remote_site ,'SB')
		GamStim_RESP_c1 = 5;
		GamStim_RESP_c2 = 6;
	end

	resp1 = 0;
	response_onset=0;
	response = zeros(24,1)-1;
	ant_win = 0;
	ant_win_onset=0;
	ant_loss = 0;
	ant_loss_onset=0;
	out_winwin = 0;
	out_winwin_onset=0;
	out_winloss = 0;
	out_winloss_onset=0;
	out_losswin = 0;
	out_losswin_onset=0;
	out_lossloss =0;
	out_lossloss_onset=0;
	baseline = 0;
	baseline_onset=0;
	error_trial = 0;
	error_trial_onset=0;

	% offset
	Condition =         find(ismember(head3{1}, 'Condition')==1);	       %1
	GamStim_OnsetTime = find(ismember(head3{1}, 'GamStim.OnsetTime')==1);  %27 
	GamStim_RESP  =     find(ismember(head3{1}, 'GamStim.RESP')==1);       %28 
	GamStim_RT    =     find(ismember(head3{1}, 'GamStim.RT')==1);         %29
	RewardShuffleImage_OnsetTime = find(ismember(head3{1}, 'RewardShuffleImage.OnsetTime')==1);  %46
	FeedbackR_OnsetTime          = find(ismember(head3{1}, 'FeedbackR.OnsetTime')==1);     %22 %16
	FeedbackRN_OnsetTime         = find(ismember(head3{1}, 'FeedbackRN.OnsetTime')==1);          %19
	LossShuffleImage_OnsetTime   = find(ismember(head3{1}, 'LossShuffleImage.OnsetTime')==1);    %40 
	FeedbackL_OnsetTime          = find(ismember(head3{1}, 'FeedbackL.OnsetTime')==1);           %10
	FeedbackLN_OnsetTime         = find(ismember(head3{1}, 'FeedbackLN.OnsetTime')==1);          %13


	first_scan = str2num(char(content3{1,1}(GamStim_OnsetTime)));


	for i=1:24
	   % trial_scan_time = str2num(char(content3{1,1}(GamStim_OnsetTime)));
	    
	    flag = 0;
	    if isempty(str2num(char(content3{i}(GamStim_RESP)))) == 1
	        
	        flag = 1;
	        error_trial = error_trial +1;
	        error_trial_onset(error_trial) = str2num(char(content3{i}(GamStim_OnsetTime)))  - first_scan;
	        baseline = baseline + 1;
	        baseline_onset(baseline) = 17000 + str2num(char(content3{i}(GamStim_OnsetTime)))  - first_scan;
	        
	    end
	    
	    if str2num(char(content3{i,1}(GamStim_RESP))) == GamStim_RESP_c1
	        %response(i) =0;
	        resp1 = resp1 + 1;
	        response_onset(resp1) = str2num(char(content3{i}(GamStim_OnsetTime))) - first_scan;
	        baseline = baseline + 1;
	        baseline_onset(baseline) = 17000 + str2num(char(content3{i}(GamStim_OnsetTime))) - first_scan;
	        RT(resp1) = str2num(char(content3{i}(GamStim_RT)));
	    end
	    if str2num(char(content3{i}(GamStim_RESP))) == GamStim_RESP_c2
	        %response(i) =1;
	        resp1 = resp1 + 1;
	        response_onset(resp1) = str2num(char(content3{i}(GamStim_OnsetTime))) - first_scan;
	        baseline = baseline + 1;
	        baseline_onset(baseline) = 17000 + str2num(char(content3{i}(GamStim_OnsetTime))) - first_scan;
	        RT(resp1) = str2num(char(content3{i}(GamStim_RT)));
	    end
	    
	    
	    if flag == 0;
	        if strcmp(	content3{i}(Condition), 'RewardWin') == 1
	            %     stimulus(i) = 1;
	            %     outcome(i) = 1;
	            %
	            %     if response(i) == 1
	            %     correct(i) = 1;
	            %     end
	            %     if response(i) == 0
	            %     correct(i) = 0;
	            %     end
	      		 Value_R(resp1) = 0.5;
	            PE_R(resp1) = 0.5;
	            ant_win = ant_win + 1;
	            out_winwin = out_winwin+1;
	            ant_win_onset(ant_win) = str2num(char(content3{i}(RewardShuffleImage_OnsetTime))) - first_scan;
		       	out_winwin_onset(out_winwin) = str2num(char(content3{i}(FeedbackR_OnsetTime ))) - first_scan;
	       		anticipation_onset(resp1) = str2num(char(content3{i}(RewardShuffleImage_OnsetTime))) - first_scan;
	            outcome_onset(resp1) = str2num(char(content3{i}(FeedbackR_OnsetTime ))) - first_scan;
			 end
	        
	        if strcmp(content3{i}(Condition), 'RewardNeu') == 1
	            %     stimulus(i) = 1;
	            %     outcome(i) = 0;
	            %
	            %     if response(i) == 0
	            %     correct(i) = 1;
	            %     end
	            %     if response(i) == 1
	            %     correct(i) = 0;
	            %     end
	            Value_R(resp1) = 0.5;
	            PE_R(resp1) = -0.5;
	            ant_win = ant_win + 1;
	            out_winloss = out_winloss+1;
	            ant_win_onset(ant_win) = str2num(char(content3{i}(RewardShuffleImage_OnsetTime))) - first_scan;
	            out_winloss_onset(out_winloss) = str2num(char(content3{i}(FeedbackRN_OnsetTime  ))) - first_scan;
	            anticipation_onset(resp1) = str2num(char(content3{i}(RewardShuffleImage_OnsetTime))) - first_scan;
	            outcome_onset(resp1) = str2num(char(content3{i}(FeedbackRN_OnsetTime ))) - first_scan;
	            
	        end
	        
	        if strcmp(content3{i}(Condition), 'LossNeu') == 1
	            %     stimulus(i) = 0;
	            %     outcome(i) = 0;
	            %
	            %     if response(i) == 1
	            %     correct(i) = 1;
	            %     end
	            %     if response(i) == 0
	            %     correct(i) = 0;
	            %     end
	            Value_R(resp1) = -0.5;
	            PE_R(resp1) = 0.25;
	            ant_loss = ant_loss + 1;
	            out_losswin = out_losswin+1;
	            ant_loss_onset(ant_loss) = str2num(char(content3{i}(LossShuffleImage_OnsetTime))) - first_scan;
	            out_losswin_onset(out_losswin) = str2num(char(content3{i}(FeedbackLN_OnsetTime)))- first_scan;
	            anticipation_onset(resp1) = str2num(char(content3{i}(LossShuffleImage_OnsetTime))) - first_scan;
	            outcome_onset(resp1) = str2num(char(content3{i}(FeedbackLN_OnsetTime ))) - first_scan;
	        end
	        
	        if strcmp(content3{i}(Condition), 'LossLose') == 1
	            %     stimulus(i) = 0;
	            %     outcome(i) = -1;
	            %
	            %     if response(i) == 0
	            %     correct(i) = 1;
	            %     end
	            %     if response(i) == 1
	            %     correct(i) = 0;
	            %     end
	            %
	            Value_R(resp1) = -0.5;
	            PE_R(resp1) = -0.25;
	            ant_loss = ant_loss + 1;
	            out_lossloss = out_lossloss+1;
	            ant_loss_onset(ant_loss) = str2num(char(content3{i}(LossShuffleImage_OnsetTime))) - first_scan;
	            out_lossloss_onset(out_lossloss) = str2num(char(content3{i}(FeedbackL_OnsetTime)))- first_scan;
	            anticipation_onset(resp1) = str2num(char(content3{i}(LossShuffleImage_OnsetTime))) - first_scan;
	            outcome_onset(resp1) = str2num(char(content3{i}(FeedbackL_OnsetTime ))) - first_scan;
	        end
	    end
	end

	names=cell(5,1);
	onsets=cell(5,1);
	durations=cell(5,1);


	names{1}='response';
	names{2}='anticipation';
	names{3}='outcome';
	names{4}='baseline';
	names{5}='error';


	onsets{1} = response_onset/1000;
	onsets{2} = anticipation_onset/1000;
	onsets{3} = outcome_onset/1000;
	onsets{4} = (response_onset/1000)+17;

	if error_trial > 0
		onsets{5} = error_trial_onset/1000;
	else
		onsets{5} = [500];
	end


	pmod = struct('name', {' '}, 'param', {}, 'poly',{});
	pmod(2).name{1} = 'anti';
	pmod(2).param{1} = Value_R;
	pmod(2).poly{1} = 1;

	pmod(3).name{1} = 'signedPE';
	pmod(3).param{1} = PE_R;
	pmod(3).poly{1} = 1;

	%pmod(3).name{2} = 'outcome_PE';
	%pmod(3).param{2} = PE_O;
	%pmod(3).poly{2} = 1;


	durations{1}=[4];
	durations{2}=[6];
	durations{3}=[1];
	durations{4}=[3];
	if error_trial > 0
		durations{5}=[17];
	else
		durations{5}=[0];
	end

	save([directory '/nDM_',participant,'_reward.mat'],'names','onsets','durations','pmod');

	behaviour(1) = median(RT);
	behaviour(2) = error_trial;
	standard(1) = std(RT);

	key(1:2,1) = cellstr('RewardTask');
	key(1,2) = cellstr('medianRT');
	key(2,2) = cellstr('TotalErrors');
	key(1,3) = cellstr('Milliseconds');
	key(2,3) = cellstr('Errors');


	outputfiles = ([directory,'/',participant,'_Rewardbehaviour.csv']);
	fid = fopen(outputfiles, 'wt');


	fprintf(fid, '%s,', char([cellstr(key(1,1))]));
	fprintf(fid, '%s,', char([cellstr(key(1,2))]));
	fprintf(fid, '%d,', behaviour(1));
	fprintf(fid, '%d,', standard(1));
	fprintf(fid, '%s\n', char([cellstr(key(1,3))]));

	fprintf(fid, '%s,', char([cellstr(key(2,1))]));
	fprintf(fid, '%s,', char([cellstr(key(2,2))]));
	fprintf(fid, '%d,', behaviour(2));
	fprintf(fid, 'N/A,');
	fprintf(fid, '%s\n', char([cellstr(key(2,3))]));

	fclose(fid);

