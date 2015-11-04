function ert_eprime2dm(eprime_file)
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
	
	[st data]= system(['eprime2csv ', eprime_file, ' 1']);
	[head content] = parse_csv(data);
	[st2 data2]= system(['eprime2csv ', eprime_file, ' 2']);
	[head2 content2] = parse_csv(data2);



	% figure out offsets
	procedure =  find(ismember(head2{1}, 'Procedure')==1); % 10; %'Procedure'
	congruency = find(ismember(head2{1}, 'congruency')==1); % 12; %'congruency'
	stim_RT =    find(ismember(head2{1}, 'stim.RT')==1); % 27; % 'stim.RT'
	fixation_RT =    find(ismember(head2{1}, 'fixation.RT')==1); % 19 % 'fixation.RT'
	stim_OnsetTime = find(ismember(head2{1}, 'stim.OnsetTime')==1); % 25;  %'stim.OnsetTime'
	pictureOffs =        find(ismember(head2{1}, 'Picture')==1); % 9; % 'Picture'
	stim_RESP =      find(ismember(head2{1}, 'stim.RESP')==1); % 26; %'stim.RESP'
	fixation_RESP =      find(ismember(head2{1}, 'fixation.RESP')==1); % 18; %'fixation.RESP'
	tasktitle_OnsetTime= find(ismember(head{1}, 'tasktitle.OnsetTime')==1); % 38/39

	% some constants that can change based on site
	start_time_offset = 0;
	H_resp = 1;
	F_resp = 2;

	% if stonybrook
	if strcmp('SB',remote_site )
		start_time_offset = 6000;
		H_resp = 5;
		F_resp = 6;
	end

	H_count = 0;
	F_count = 0;
	iC_count = 0;
	iI_count = 0;
	cI_count = 0;
	cC_count = 0;
	misc_count = 0;

	H_onset = 0;
	F_onset = 0;
	iC_onset = 0;
	iI_onset = 0;
	cI_onset = 0;
	cC_onset = 0;
	misc_onset = 0;
	error_onset = 0;
	post_error_onset = 0;

	error_count = 0;
	miss_count = 0;
	post_error_count = 0;
	post_error = 0;
	miss = 0;
	err = 0;
	start_time = (str2num(char(content{1,1}(tasktitle_OnsetTime)))+start_time_offset); 


	iC_error = 0;
	iI_error = 0;
	cI_error = 0;
	cC_error = 0;


	for j = 1:length(content2)-1
	    
	    post_error = 0;
	    
	    if err == 1 && (strcmp(content2{j,1}(congruency), '') ~= 1) && (strcmp(cellstr(content2{j,1}(procedure)), 'paradigm') == 1)
	        post_error = 1;
	        post_error_count = post_error_count + 1;
	    end
	    
	    
	    err = 0;
	    miss = 0;
	    
	    stimRT = str2num(char(content2{j,1}(stim_RT)));
	    fixRT = str2num(char(content2{j,1}(fixation_RT)));
	    stim_onset = str2num(char(content2{j,1}(stim_OnsetTime)));
	    
	    if strcmp(char(content2{j,1}(procedure)), 'rest') == 0
	        if strcmp(char(content2{j,1}(procedure)), 'paradigm') == 1
	            picture = char((content2{j,1}(pictureOffs)));
	            l_pic = length(picture);
	            
	            if strcmp(picture(l_pic-2), 'H') == 1
	                happy(j) = 1;
	                
	                if isempty(char(content2{j,1}(stim_RESP))) == 0
	                    if str2num(char(content2{j,1}(stim_RESP))) == H_resp
	                        error_count = error_count + 1;
	                        err = 1;
	                    end
	                end
	                
	                if (isempty(char(content2{j,1}(fixation_RESP))) == 0) && (isempty(char(content2{j,1}(stim_RESP))) == 1)
	                    if str2num(char(content2{j,1}(fixation_RESP))) == H_resp
	                        error_count = error_count + 1;
	                        err = 1;
	                    end
	                end
	                
	                if (isempty(char(content2{j,1}(fixation_RESP))) == 1) && (isempty(char(content2{j,1}(stim_RESP))) == 1)
	                    miss_count = miss_count + 1;
	                    miss = 1;
	                    error_count = error_count + 1;
	                    err = 1;
	                end
	                
	                
	                
	            end
	            if strcmp(picture(l_pic-2), 'F') == 1
	                happy(j) = -1;
	                
	                
	                if isempty(char(content2{j,1}(stim_RESP))) == 0
	                    if str2num(char(content2{j,1}(stim_RESP))) == F_resp
	                        error_count = error_count + 1;
	                        err = 1;
	                    end
	                end
	                
	                if (isempty(char(content2{j,1}(fixation_RESP))) == 0) && (isempty(char(content2{j,1}(stim_RESP))) == 1)
	                    if str2num(char(content2{j,1}(fixation_RESP))) == F_resp
	                        error_count = error_count + 1;
	                        err = 1;
	                    end
	                end
	                
	                if (isempty(char(content2{j,1}(fixation_RESP))) == 1) && (isempty(char(content2{j,1}(stim_RESP))) == 1)
	                    miss_count = miss_count + 1;
	                    miss = 1;
	                    error_count = error_count + 1;
	                    err = 1;
	                end
	                
	            end
	            
	            
	            if (err == 1) && (post_error == 1)
	                
	                post_error = 0;
	                post_error_count = post_error_count - 1;
	                
	            end
	            
	            if post_error == 1
	                post_error_onset(post_error_count) = stim_onset - start_time;
	            end
	            
	        end
	        
			if strcmp(remote_site ,'SB')
				if happy(j) == 1
			       if (err == 0) && (miss == 0) && (post_error == 0)
			            H_count = H_count + 1;
			            H_onset(H_count) = stim_onset - start_time;
			       end 
			        
			    end
			    if happy(j) == -1
			       if (err == 0) && (miss == 0) && (post_error == 0)
			            F_count = F_count + 1;
			            F_onset(F_count) = stim_onset - start_time;
			       end 
			        
			    end
			end
	        if strcmp(char(content2{j,1}(congruency)), 'ci') == 1
	            trialtype(j) = 2;
	            
	            if (err == 0) && (miss == 0) && (post_error == 0)
	                cI_count = cI_count + 1;
	                cI_onset(cI_count) = stim_onset - start_time;
	                
	                if stimRT == 0
	                    cI_RT(cI_count) = fixRT + 1000;
	                end
	                if stimRT > 30
	                    cI_RT(cI_count) = stimRT;
	                end
	                
	                if (0 >stimRT >= 30)
	                    cI_count = cI_count - 1;
	                end
	                
	            end
	            
	            
	            if err == 1
	                cI_error = cI_error+1;
	            end
	            %if miss == 1
	            %cI_error = cI_error+1;
	            %end
	        end
	        
	        
	        
	        if strcmp(char(content2{j,1}(congruency)), 'cc') == 1
	            trialtype(j) = 1;
	            
	            if (err == 0) && (miss == 0) && (post_error == 0)
	                cC_count = cC_count + 1;
	                cC_onset(cC_count) = stim_onset - start_time;
	                
	                if stimRT ==0
	                    cC_RT(cC_count) = fixRT + 1000;
	                end
	                if stimRT > 30
	                    cC_RT(cC_count) = stimRT;
	                end
	                if (0 >stimRT >= 30)
	                    cC_count = cC_count - 1;
	                end
	            end
	            
	            if err == 1
	                cC_error = cC_error+1;
	            end
	            %if miss == 1
	            %cC_error = cC_error+1;
	            %end
	        end
	        
	        
	        if strcmp(char(content2{j,1}(congruency)), 'ic') == 1
	            trialtype(j) = 1;
	            
	            if (err == 0) && (miss == 0) && (post_error == 0)
	                iC_count = iC_count + 1;
	                iC_onset(iC_count) = stim_onset - start_time;
	                
	                if stimRT == 0
	                    iC_RT(iC_count) = fixRT + 1000;
	                end
	                if stimRT > 30
	                    iC_RT(iC_count) = stimRT;
	                end
	                
	                if (0 >stimRT >= 30)
	                    iC_count = iC_count - 1;
	                end
	            end
	            
	            if err == 1
	                iC_error = iC_error+1;
	            end
	            %if miss == 1
	            %iC_error = iC_error+1;
	            %end
	            
	            
	            
	        end
	        
	        
	        
	        if strcmp(char(content2{j,1}(congruency)), 'ii') == 1
	            trialtype(j) = 3;
	            
	            if (err == 0) && (miss == 0) && (post_error == 0)
	                iI_count = iI_count + 1;
	                iI_onset(iI_count) = stim_onset - start_time;
	                
	                if stimRT == 0
	                    iI_RT(iI_count) = fixRT + 1000;
	                end
	                if stimRT > 30
	                    iI_RT(iI_count) = stimRT;
	                end
	                
	                if (0 >stimRT >= 30)
	                    iI_count = iI_count - 1;
	                end
	                
	                
	            end
	            
	            if err == 1
	                iI_error = iI_error+1;
	            end
	            %if miss == 1
	            %iI_error = iI_error+1;
	            %end
	        end
	        
	        if err == 1
	            error_onset(error_count) = stim_onset - start_time;
	        end
	        
	        
	        if (strcmp(content2{j,1}(congruency), '') == 1) && (strcmp(cellstr(content2{j,1}(procedure)), 'paradigm') == 1)
	            trialtype(j) = -1;
	            
	            if err == 0
	                if miss == 0
	                    misc_count = misc_count + 1;
	                    misc_onset(misc_count) = stim_onset - start_time;
	                end
	            end
	        end
	    end
	end


	names=cell(7,1);
	onsets=cell(7,1);
	durations=cell(7,1);
	%
	%
	names{1}='cC';
	names{2}='cI';
	names{3}='iC';
	names{4}='iI';
	names{5}='error';
	names{6}='posterror';
	names{7}='misc';

	onsets{1} = cC_onset/1000;
	onsets{2} = cI_onset/1000;
	onsets{3} = iC_onset/1000;
	onsets{4} = iI_onset/1000;
	onsets{5} = error_onset/1000;
	onsets{6} = post_error_onset/1000;
	onsets{7} = misc_onset/1000;


	%
	durations{1}=[0];
	durations{2}=[0];
	durations{3}=[0];
	durations{4}=[0];
	durations{5}=[0];
	durations{6}=[0];
	durations{7}=[0];

	% save design matrix
	save([directory '/nDM_',participant,'_ert.mat'],'names','onsets','durations');


	names=cell(4,1);
	onsets=cell(4,1);
	durations=cell(4,1);
	
	names{1}='Happy';
	names{2}='Fear';
	names{3}='error';
	names{4}='posterror';
	%names{5}='misc';
	 
	onsets{1} = H_onset/1000;
	onsets{2} = F_onset/1000;
	onsets{3} = error_onset/1000;
	onsets{4} = post_error_onset/1000;
	%onsets{5} = misc_onset/1000;
	
	durations{1}=[0];
	durations{2}=[0];
	durations{3}=[0];
	durations{4}=[0];
	%durations{5}=[0];

	save([directory '/nDM_',participant,'_ert2.mat'],'names','onsets','durations');


	RTm(1) = mean(cC_RT);
	RTm(2) = mean(cI_RT);
	RTm(3) = mean(iC_RT);
	RTm(4) = mean(iI_RT);
	RTs(1) = std(cC_RT);
	RTs(2) = std(cI_RT);
	RTs(3) = std(iC_RT);
	RTs(4) = std(iI_RT);
	RThigh(1) = RTm(1) + (RTs(1)*2);
	RThigh(2) = RTm(2) + (RTs(2)*2);
	RThigh(3) = RTm(3) + (RTs(3)*2);
	RThigh(4) = RTm(4) + (RTs(4)*2);
	RTlow(1) = RTm(1) - (RTs(1)*2);
	RTlow(2) = RTm(2) - (RTs(2)*2);
	RTlow(3) = RTm(3) - (RTs(3)*2);
	RTlow(4) = RTm(4) - (RTs(4)*2);

	cC_count_final = 0;
	cI_count_final = 0;
	iC_count_final = 0;
	iI_count_final = 0;

	cC_RT_final = 0;
	cI_RT_final = 0;
	iC_RT_final = 0;
	iI_RT_final = 0;

	for zz = 1:cC_count
	    
	    if (cC_RT(zz) < RThigh(1)) && (cC_RT(zz) > RTlow(1))
	        cC_count_final = cC_count_final + 1;
	        cC_RT_final(cC_count_final) = cC_RT(zz);
	    end
	    
	end

	for zz = 1:cI_count
	    
	    if (cI_RT(zz) < RThigh(2)) && (cI_RT(zz) > RTlow(2))
	        cI_count_final = cI_count_final + 1;
	        cI_RT_final(cI_count_final) = cI_RT(zz);
	    end
	    
	end

	for zz = 1:iC_count
	    
	    if (iC_RT(zz) < RThigh(3)) && (iC_RT(zz) > RTlow(3))
	        iC_count_final = iC_count_final + 1;
	        iC_RT_final(iC_count_final) = iC_RT(zz);
	    end
	    
	end

	for zz = 1:iI_count
	    
	    if (iI_RT(zz) < RThigh(4)) && (iI_RT(zz) > RTlow(4))
	        iI_count_final = iI_count_final + 1;
	        iI_RT_final(iI_count_final) = iI_RT(zz);
	    end
	    
	    
	end

	RTm_final(1) = mean(cC_RT_final);
	RTm_final(2) = mean(cI_RT_final);
	RTm_final(3) = mean(iC_RT_final);
	RTm_final(4) = mean(iI_RT_final);

	behaviour(1) = iI_error - cI_error;
	behaviour(2) = (iI_error + cI_error) - (iC_error + cC_error);
	behaviour(3) = RTm_final(4) - RTm_final(2);
	behaviour(4) = (RTm_final(4) + RTm_final(2)) - (RTm_final(1) + RTm_final(3));
	standard(3) = (std(iI_RT_final)+std(cI_RT_final))/2;
	standard(4) = (std(cC_RT_final)+std(cI_RT_final)+std(iC_RT_final)+std(iI_RT_final))/4;
	key(1:4,1) = cellstr('EmotionRecognitionTask');
	key(1,2) = cellstr('iI_minus_cI_Error_difference');
	key(2,2) = cellstr('I_minus_C_Error_difference');
	key(3,2) = cellstr('iI_minus_cI_RTdifference');
	key(4,2) = cellstr('I_minus_C_RTdifference');
	key(1,3) = cellstr('Errors');
	key(2,3) = cellstr('Errors');
	key(3,3) = cellstr('Milliseconds');
	key(4,3) = cellstr('Milliseconds');


	full_behaviour(1) = (148 - (cC_error + cI_error + iC_error + iI_error))/148;
	full_behaviour(2) = (37 - iI_error)/37;
	full_behaviour(3) = (36 - iC_error)/36;
	full_behaviour(4) = (36 - cI_error)/36;
	full_behaviour(5) = (35 - cC_error)/35;
	full_behaviour(6) = behaviour(3);
	full_behaviour(7) = (RTm_final(2) + RTm_final(4))/2;
	full_behaviour(8) = (RTm_final(1) + RTm_final(3))/2;
	full_behaviour(9) = RTm_final(4);
	full_behaviour(10) = RTm_final(3);
	full_behaviour(11) = RTm_final(2);
	full_behaviour(12) = RTm_final(1);
	full_behaviour(13) = iI_count_final;
	full_behaviour(14) = iC_count_final;
	full_behaviour(15) = cI_count_final;
	full_behaviour(16) = cC_count_final;

	outputfiles = ([directory,'/',participant,'_ERT_behaviour.csv']);
	fid = fopen(outputfiles, 'wt');
	
	for i = 1:2
		fprintf(fid, '%s,', char([cellstr(key(i,1))]));
	  	fprintf(fid, '%s,', char([cellstr(key(i,2))]));
	 	fprintf(fid, '%d,', behaviour(i));
	 	fprintf(fid, 'N/A,');
	 	fprintf(fid, '%s\n', char([cellstr(key(i,3))]));
	 end

	 for i = 3:4
		 fprintf(fid, '%s,', char([cellstr(key(i,1))]));
		 fprintf(fid, '%s,', char([cellstr(key(i,2))]));
		 fprintf(fid, '%d,', behaviour(i));
		 fprintf(fid, '%d,', standard(i));
		 fprintf(fid, '%s\n', char([cellstr(key(i,3))]));
	 end

	 fclose(fid);

	outputfiles = ([directory,'/',participant,'_ERT_behaviour_full.csv']);
	fid = fopen(outputfiles, 'wt');
	
	for i = 1:16
		fprintf(fid, '%d\n,', full_behaviour(i));
	 	
	 end


	 fclose(fid);

