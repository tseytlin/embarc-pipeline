%
% run MATLAB specific ExploreDTI code
%
function dti_sequence(exploreDTI_dir, subject)
	% define template location
	template_filename = [exploreDRI_dir 'Templates/AAL_Template.nii'];	
	labels_filename = [exploreDRI_dir 'Templates/AAL_Label_Names.txt'];

	% Resample ${EXPLOREDTI}/T1_1mm_RAS/${SUBJECT}_T1_masked_RAS.nii to [2 2 2] mm^3:
	% script: home/versacea/bin/ExploreDTI/Source/Plugins/E_DTI_Resample_nii.p 
	% (or ???  E_DTI_Resample_nii_ex.p  or E_DTI_Resample_nii_ex_me.p or 
	% E_DTI_resample_nii_file.p  ??? NOT SURE WHICH ONE IS THE CORRECT ONE) 
	% input: ${EXPLOREDTI}/T1_1mm_RAS/${SUBJECT}_T1_masked_RAS.nii 
	% destination/output: ${EXPLOREDTI}/REORDERED/${SUBJECT}_T1.nii
	% voxel size: [2 2 2]
	in_t1  = [ exploreDTI_dir '/T1_1mm_RAS/' subject '_T1_masked_RAS.nii'];
	out_t1 =  [ exploreDTI_dir '/REORDERED/' subject '_T1.nii']; 
	vox_sz = [2,2,2];
	interp = 2; % Interpolator: 0=nearest neighbor, 1= trilinear, 2 = cubic B-spline	
	E_DTI_Resample_nii_ex(int_t1,out_t1,vox_sz,interp);


	% Generate DWI.mat with customized parameters
	% Convet bvec and bval in B-matrix 
	% script: home/versacea/bin/ExploreDTI/Source/Plugins/E_DTI_convert_nii_dic_2_txt.p
	% input and output destination folder: ${EXPLOREDTI}/REORDERED/.
	bvec = [exploreDTI_dir '/REORDERED/' subject '.bvec';
	bval = [exploreDTI_dir '/REORDERED/' subject '.bval'];
	txt  = [exploreDTI_dir '/REORDERED/' subject '.txt'];
	E_DTI_convert_nii_dic_2_txt_exe(bvec,bval,txt);


	% Calculate DTI*.mat file --> convert ${EXPLOREDTI}/REORDERED/${SUBJECT}.nii to 
	% ${EXPLOREDTI}/REORDERED/${SUBJECT}.mat file with the follwoing parameters  
	% [permute (y x z); flip (-x y z) B value (1000); voxel size (2 2 2); N no-DWI (7);
	% N DWI (61); matrix size (128 128 64)] 
	% script: E_DTI_convert_nii_dic_2_txt.p
	% input and output destination folder: ${EXPLOREDTI}/REORDERED/.
	% ?? Script to Convert folder of DICOMs to an ExploreDTI .mat file:
	% script syntax: E_DTI_Script_Get_DTI_folders(source_dir,target_dir)
	% source_dir = /data/Phillips2/projects/dtistudy/BIOS/data/exploreDTI/REORDERED
	% target_dir = /data/Phillips2/projects/dtistudy/BIOS/data/exploreDTI/REORDERED
	% NOT SURE HOW YOU CAN ENTER THE FOLLOWING PARAMETERS (NEXT SLIDE)?
	source_dir = [exploreDTI_dir '/REORDERED/'];
	target_dir = [exploreDTI_dir '/REORDERED/'];
	E_DTI_Script_Get_DTI_folders(source_dir,target_dir);


	% Setting the parameters/ Correcting the DWIs for Subject Motion (SM), Eddy Current (EC) and EPI distortions
	% script: E_DTI_SMECEPI_Main(parameter_filename);
	% destination folder: /data/Phillips2/projects/dtistudy/BIOS/data/exploreDTI/MC_D_EPI/ 
	parameter_file = ['']; % hard-coded parameter file that I got for everyone
	E_DTI_SMECEPI_Main(parameter_filename);
	
	% copy results to different folder
	dos([ 'mv ' exploreDTI_dir '/MC_D_EPI/*trafo.mat  ' exploreDTI_dir '/MC_D_EPI/trafo/']);
	dos([ 'mv ' exploreDTI_dir '/MC_D_EPI/*native.mat  ' exploreDTI_dir '/MC_D_EPI/native/']);
	

	% export DWI.mat to bvec, bval
	% figure out parameters
	% DWIs with B0(s) ('_DWIs.nii')
	% Gradient directions ('_grads.txt')	
	E_DTI_Export_mat_2_analyze();

	% Here is an example of the script needed for each subject:
	% All the subjects will need a dmrirc file. 
	% I usually, just sed ‘s/20120303.1953/${SUBEJECT}/g’
	
	% Going back to ExploreDTI (in parallel) run whole-brain tractography (DTI)
	% With the following parameters:
	paramaters.SeedPointRes = [3 3 3];
	paramaters.StepSize = 1;
	paramaters.FAThresh = 0.2000;
	paramaters.AngleThresh = 45;
	paramaters.FiberLengthRange = [50 500];
	filename_in = [exploreDTI_dir '/MD_C_EPI/trafo/' ];
	filename_out = [exploreDTI_dir '/MD_C_EPI/trafo/'];
	
	% run determenistic tractography
	WholeBrainTrackingDTI(filename_in, filename_out, parameters);

	% create connectivity matricies
	dti_filename = [exploreDTI_dir '/MD_C_EPI/trafo/*MD_C_trafo.mat'];
	tracts_filename = [exploreDTI_dir '/MD_C_EPI/trafo/*MD_C_trafo_Tracts_DTI.mat'];
	img_filename = [];  % for now	
	% 3D volume of the atlas (e.g. AAL) reference template. 
	% Note this volume must be read into MATLAB using the command:A_T=E_DTI_load_nii	
	A_T = template_filename;
	% 3D volume of the atlas (e.g. AAL) labels. As above, the volume need to be loaded into MATLABA_L=E_DTI_load_nii
	A_L = labels_filename;
	% cell array of label names for each region in A_L.
	L = {};	
	% /data/Phillips2/projects/dtistudy/ExploreDTI/Templates/AAL_Label_Names.txt
	% Directory to save the connectivity matrices	
	mat_dir = [ exploreDTI_dir '/CMs'];
	img_suffix =  [] ; % for now.
	% a n_labelsx2 matrix of pairs of label indices to create an individual tract .mat file for that pair of labels.
	selected_labels = [];
	% Use [] to disable. 
	% Use nchoosek(1:n_atlas_labels,2) to output for every pair of labels # for now. 
	ACh = 3;

	% call script
	E_DTI_Network_analysis_exe(dti_filename,tracts_filename,img_filename,A_T,A_L,L,VDims_A,mat_dir,img_suffix,selected_labels,ACh);







