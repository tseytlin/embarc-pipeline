#!/usr/bin/env python
# EMBARC 2.0 Analysis Pipeline (resting)
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import sys
import os                                  
import re
import gold

## Predefined constants ##
conf = gold.Config()
conf.time_repetition  = 2.0
conf.CPU_CORES = 20
conf.fugue_dwell_time = 0.00025	
conf.filter_image_bptf = ' -bptf 50 1.4'	


conf.bet_mask = True
conf.bet_frac = 0.5 #0.6
conf.bet_robust = True
conf.bet_vertical_gradient = 0
		
conf.flirt_cost = 'mutualinfo'
conf.flirt_bins = 256
conf.flirt_dof = 12
conf.flirt_interp = 'trilinear'
conf.flirt_searchr_x = [-180, 180]
conf.flirt_searchr_y = [-180, 180]
conf.flirt_searchr_z = [-180, 180]

conf.coregister_cost_function = "nmi"
conf.coregister_separation = [4, 2]
conf.coregister_tolerance = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001]
conf.coregister_fwhm = [7, 7]
conf.prepare_fieldmap_scanner = "SIEMENS"
conf.prepare_fieldmap_delta_TE = 2.46  

conf.fugue_poly_order = 3

conf.dartel_fwhm = 6 #TODO Check value
conf.dartel_voxel_size =  (2, 2, 2)
conf.susan_brightness_threshold = 750 #200.0
conf.susan_fwhm = 6
		
conf.modelspec_concatenate_runs   = False
conf.modelspec_high_pass_filter_cutoff = 60 # reward only
conf.modelspec_input_units = 'secs'
		
conf.level1design_bases = {'hrf':{'derivs': [0,0]}}
conf.level1design_timing_units = 'secs'
conf.level1estimate_estimation_method = {'Classical' : 1}
conf.contrastestimate_use_derivs = True
conf.level1design_microtime_onset = 1
conf.level1design_microtime_resolution = 16
conf.level1design_model_serial_correlations = 'AR(1)'




useFieldmap=True
noPrint = False

"""
EMBARC 1.0 Input Data Source
"""
def datasource(directory, sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o
	import glob as gl
	
	# define some variables beforehand
	subject=get_subject(directory)
	
	# define templates for datasource
	field_template = dict(func=sequence+'/*.nii',struct='*mprage/*_crop.nii')	
	template_args  = dict(func=[[]],struct=[[]])                
	outfields=['func', 'struct']       


	#DMH: changed fieldmap to be resliced ones (to match dimensions of EPI), but this didn't work - commented out
	if useFieldmap:
			#field_template['fieldmap_mag']   = "field_map/magnitude/mag_resliced.nii"
			#field_template['fieldmap_phase'] ="field_map/phase/phase_resliced.nii"
			field_template['fieldmap_mag']   = "field_map/magnitude/*.nii"
			field_template['fieldmap_phase'] ="field_map/phase/*.nii"
			template_args['fieldmap_mag']  = [[]]
			template_args['fieldmap_phase']  = [[]]
			outfields.append('fieldmap_mag')
			outfields.append('fieldmap_phase')


	# specify input dataset just pass through parameters
	datasource = pe.Node(interface=nio.DataGrabber(
						 infields=['subject_id','sequence'], 
						 outfields=outfields),
	                     name = 'datasource')
	datasource.inputs.base_directory = os.path.abspath(directory)
	datasource.inputs.template = '*'
	datasource.inputs.field_template = field_template
	datasource.inputs.template_args  = template_args
	datasource.inputs.subject_id = subject
	datasource.inputs.sequence = sequence
	datasource.inputs.sort_filelist = True

	return datasource

"""
EMBARC get subject name from a directory
"""
def get_subject(directory):
	m = re.search('Ss([0-9]+)',directory)
	if m:
		return m.group(1)
	return "subject"

def subset(x,i):
	return x[i]	


"""
EMBARC/DIAMOND 2.0 PreProcessing Pipeline
input: 
	func  - functional image
	struct - structural image
	template - template image
	fieldmap_mag - fieldmap magnitute image
	fieldmap_phase	- phase image
output:
	func  - functional processed images
	ufunc - unsmoothed functional image
	mask  - mask image from the functional
	movement - realign movement parameters
	struct - structural processed image
"""
def preprocess_mni(config,useFieldmap=True,name='preprocess2'):
	import nipype.interfaces.spm as spm          # spm
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.fsl.maths as math  
	import nipype.interfaces.utility as util     # utility
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.afni as afni	     # afni
	import nipype.workflows.fmri.spm.preprocess as dartel # preprocess
	from nipype.interfaces.nipy.preprocess import Trim

	fsl.FSLCommand.set_default_output_type('NIFTI')
	
	
	preproc = pe.Workflow(name=name)
	inputFields = ['func','struct']
	if useFieldmap:
		inputFields = ['func','struct','fieldmap_mag','fieldmap_phase']
		
	inputnode = pe.Node(interface=util.IdentityInterface(fields=inputFields),name='input')

	# realign 4D functional
	realign = pe.Node(interface=spm.Realign(), name="realign")
	realign.inputs.register_to_mean = True
	preproc.connect(inputnode,"func",realign,"in_files")

	# skull strip mean functional image
	bet_mean = pe.Node(interface=fsl.BET(), name="bet_mean")
	bet_mean.inputs.mask = config.bet_mask
	bet_mean.inputs.frac = config.bet_frac
	bet_mean.inputs.robust = config.bet_robust
	bet_mean.inputs.vertical_gradient = config.bet_vertical_gradient
	preproc.connect(realign,'mean_image',bet_mean,'in_file') 

	# skull strip mean structural image
	bet_struct = pe.Node(interface=fsl.BET(), name="bet_struct")
	bet_struct.inputs.mask = config.bet_mask
	bet_struct.inputs.frac = config.bet_frac
	bet_struct.inputs.robust = config.bet_robust
	bet_struct.inputs.vertical_gradient = config.bet_vertical_gradient
	preproc.connect(inputnode,'struct',bet_struct,'in_file')	


	# coregister images
	coreg_func2struct = pe.Node(interface=spm.Coregister(),name="coreg_func2struct")
	coreg_func2struct.inputs.jobtype = "estimate"
 	coreg_func2struct.inputs.cost_function = config.coregister_cost_function
	coreg_func2struct.inputs.separation = config.coregister_separation
	coreg_func2struct.inputs.tolerance = config.coregister_tolerance
	coreg_func2struct.inputs.fwhm = config.coregister_fwhm
	preproc.connect(bet_struct,'out_file',coreg_func2struct,'target')
	preproc.connect(realign,'realigned_files',coreg_func2struct,'apply_to_files')
	preproc.connect(bet_mean,'out_file',coreg_func2struct,'source')


	# select first image of mag
	if useFieldmap:
		firstmag = pe.Node(interface=Trim(),name="firstmag")
		firstmag.inputs.begin_index = 0	
		firstmag.inputs.end_index = 1
		preproc.connect(inputnode,'fieldmap_mag',firstmag,'in_file')

		# select first image of phase
		#firstphase = pe.Node(interface=Trim(),name="firstphase")
		#firstphase.inputs.begin_index = 0	
		#firstphase.inputs.end_index = 1
		#preproc.connect(inputnode,'fieldmap_phase',firstphase,'in_file')

		# coregister phase image to structural
		coreg_fphase2struct = pe.Node(interface=spm.Coregister(),name="coreg_fieldmapphase2struct")
		coreg_fphase2struct.inputs.jobtype = "estimate"
		coreg_fphase2struct.inputs.cost_function = config.coregister_cost_function
		coreg_fphase2struct.inputs.separation = config.coregister_separation
		coreg_fphase2struct.inputs.tolerance = config.coregister_tolerance
		coreg_fphase2struct.inputs.fwhm = config.coregister_fwhm
		preproc.connect(bet_struct,'out_file',coreg_fphase2struct,'target')
		preproc.connect(firstmag,'out_file',coreg_fphase2struct,'source')
		preproc.connect(inputnode,'fieldmap_phase',coreg_fphase2struct,'apply_to_files')
		#preproc.connect(firstphase,'out_file',coreg_fphase2struct,'apply_to_files')
	
		# after coreg fphase, invert phase image and pipe into prepare_field
	
		# get maximum from image
		image_max = pe.Node(interface=fsl.ImageStats(),name='image_max')	
		image_max.inputs.op_string = "-R"
		preproc.connect(coreg_fphase2struct,'coregistered_files',image_max,'in_file')

		# invert image using fslmats
		invert_image = pe.Node(interface=math.MathsCommand(),name='invert_image')	
		preproc.connect(coreg_fphase2struct,'coregistered_files',invert_image,'in_file')
		preproc.connect(image_max,('out_stat',create_inversion_args, "-mul -1 -add "),invert_image,'args')	

		# coregister magnitude image to structural
		coreg_fmag2struct = pe.Node(interface=spm.Coregister(),name="coreg_fieldmapmag2struct")
		coreg_fmag2struct.inputs.jobtype = "estimate"
		coreg_fmag2struct.inputs.cost_function = config.coregister_cost_function
		coreg_fmag2struct.inputs.separation = config.coregister_separation
		coreg_fmag2struct.inputs.tolerance = config.coregister_tolerance
		coreg_fmag2struct.inputs.fwhm = config.coregister_fwhm
		preproc.connect(bet_struct,'out_file',coreg_fmag2struct,'target')
		preproc.connect(firstmag,'out_file',coreg_fmag2struct,'source')
		preproc.connect(inputnode,'fieldmap_mag',coreg_fmag2struct, 'apply_to_files')
	
		# skull strip magnitude image
		bet_mag = pe.Node(interface=fsl.BET(), name="bet_mag")
		bet_mag.inputs.mask = config.bet_mask
		bet_mag.inputs.frac = config.bet_frac
		bet_mag.inputs.robust = config.bet_robust
		bet_mag.inputs.vertical_gradient = config.bet_vertical_gradient
		preproc.connect(coreg_fmag2struct,'coregistered_source',bet_mag,'in_file')
	
		# prepare fieldmap
		prepare_field = pe.Node(interface=fsl.PrepareFieldmap(),name='prepare_fieldmap')
		prepare_field.inputs.output_type = "NIFTI"
		prepare_field.inputs.scanner = config.prepare_fieldmap_scanner
		prepare_field.inputs.delta_TE = config.prepare_fieldmap_delta_TE
		preproc.connect(bet_mag,'out_file',prepare_field,'in_magnitude')
		preproc.connect(invert_image,'out_file',prepare_field,'in_phase')

		
		#DMH: ADD BELOW TO SKIP INVERSION STEP
    		#preproc.connect(coreg_fphase2struct,'coregistered_files',prepare_field,'in_phase')	
		#reslice fieldmap: added by DMH to run with BIOS data
		#reslice_fieldmap = pe.Node(interface=fsl.ApplyXfm(), name='reslicenode')
		#reslice_fieldmap.inputs.uses_qform = True
		#reslice_fieldmap.inputs.apply_xfm = False
		#preproc.connect(prepare_field,'out_fieldmap',reslice_fieldmap,'in_file')
		#preproc.connect(inputnode,'func',reslice_fieldmap,'reference')


		# FSL FUGUE
		fugue = pe.Node(interface=fsl.FUGUE(),name='fieldmap_FUGUE')
		fugue.inputs.dwell_time = config.fugue_dwell_time
		fugue.inputs.poly_order = config.fugue_poly_order
		#DMH: Change unwarp direction (from default y) based on Anna's suggestion		
		#fugue.inputs.unwarp_direction = 'z'
		preproc.connect(coreg_func2struct,'coregistered_files',fugue,'in_file')
		preproc.connect(prepare_field,'out_fieldmap',fugue,'fmap_in_file')
		
		#DMH: Link FUGUE to resliced fieldmap		
		#preproc.connect(reslice_fieldmap,'out_file',fugue,'fmap_in_file')		

		#TODO mask from BET struct
		#TODO asym_se_time
		#TODO dwell_time = 0.79 ??? different for encore
		#TODO dwell_to_asym_ratio
	else:
		# convert image using fslmats
		convert_image = pe.Node(interface=math.MathsCommand(),name='convert_image')
		convert_image.inputs.args = "-mul 1"
		preproc.connect(coreg_func2struct,'coregistered_files',convert_image,'in_file')
		
	
	# create dartel template
	dartel_template = dartel.create_DARTEL_template()
	dartel_template.inputs.inputspec.template_prefix = 'Template'
	preproc.connect(inputnode, 'struct',dartel_template,'inputspec.structural_files')


	# now lets do normalization with DARTEL
	norm_func =  pe.Node(interface=spm.DARTELNorm2MNI(modulate=True),name='norm_func')	
	norm_func.inputs.fwhm = config.dartel_fwhm
	norm_func.inputs.voxel_size = config.dartel_voxel_size
	preproc.connect(dartel_template,'outputspec.template_file',norm_func,'template_file')
	preproc.connect(dartel_template, 'outputspec.flow_fields', norm_func, 'flowfield_files')
	if useFieldmap:
		preproc.connect(fugue,'unwarped_file',norm_func,'apply_to_files')
	else:
		preproc.connect(convert_image,'out_file',norm_func,'apply_to_files')
		#preproc.connect(coreg_func2struct,'coregistered_files',norm_func,'apply_to_files')

	# now lets do normalization with DARTEL
	norm_struct =  pe.Node(interface=spm.DARTELNorm2MNI(modulate=True),name='norm_struct')
	norm_struct.inputs.fwhm = config.dartel_fwhm  #TODO Check value
	preproc.connect(dartel_template,'outputspec.template_file',norm_struct,'template_file')
	preproc.connect(dartel_template, 'outputspec.flow_fields', norm_struct, 'flowfield_files')
	preproc.connect(bet_struct,'out_file',norm_struct,'apply_to_files')

	
	# remove spikes
	despike = pe.Node(interface=afni.Despike(), name='despike')
	despike.inputs.outputtype = 'NIFTI'
	#TODO: need params 3 and 5?  Default 2 and 4
	preproc.connect(norm_func,'normalized_files',despike,'in_file')
	
	# calculated brighness threshold for susan (mean image intensity * 0.75)
	image_mean = pe.Node(interface=fsl.ImageStats(),name='image_mean')	
	image_mean.inputs.op_string = "-M"
	#preproc.connect(bet_mean,'out_file',image_mean,'in_file')
	preproc.connect(despike,'out_file',image_mean,'in_file')


	# scale image so that mean 1000/original
	scale_image = pe.Node(interface=math.MathsCommand(),name='scale_image')
	preproc.connect(despike,'out_file',scale_image,'in_file')
	preproc.connect(image_mean,('out_stat',create_scale_args, "-mul 1000 -div "),scale_image,'args')


	# smooth image using SUSAN
	susan = pe.Node(interface=fsl.SUSAN(), name="smooth")
	susan.inputs.brightness_threshold = config.susan_brightness_threshold 
	susan.inputs.fwhm = config.susan_fwhm
	preproc.connect(scale_image,'out_file',susan,'in_file') 
	#preproc.connect(despike,'out_file',susan,'in_file') 
	#preproc.connect(image_mean,('out_stat',create_brightness_threshold),susan,'brightness_threshold') 

	# create a nice mask to output	
	bet_func = pe.Node(interface=fsl.BET(), name="bet_func")
	bet_func.inputs.mask = config.bet_mask
	bet_func.inputs.frac = config.bet_frac
	bet_func.inputs.robust = config.bet_robust
	bet_func.inputs.vertical_gradient = config.bet_vertical_gradient
	preproc.connect(susan,'smoothed_file',bet_func,'in_file')	
	            
	# gather output
	outputnode = pe.Node(interface=util.IdentityInterface(fields=['func','ufunc','mask','movement','struct']),name='output')
	preproc.connect(despike,'out_file',outputnode, 'ufunc')
	preproc.connect(susan,'smoothed_file',outputnode,'func')
	preproc.connect(realign,'realignment_parameters',outputnode,'movement')
	preproc.connect(norm_struct,'normalized_files',outputnode,'struct')
	preproc.connect(bet_func,'mask_file',outputnode,'mask')
	
	return preproc



"""
EMBARC 1.0 Resting Sequence Ex: resting1/resting2
directory - dataset directory
sequence  - name of the sequence
subject   - optional subject name if None, embarc subject will be derived
ds		  - DataSource node for this dataset, if None embarc will be used

"""
def resting(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import CPAC					# import CPAC nuisance
	import nipype.interfaces.fsl as fsl          # fsl
	import wrappers as wrap	
	import nipype.interfaces.afni as afni		 # afni
	from nipype.interfaces.utility import Function
	import nipype.interfaces.io as nio           # Data i/o
	import nipype.algorithms.misc as misc
	import gold	
	
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	subject = get_subject(directory)
	
	conf.modelspec_high_pass_filter_cutoff = 256
	
	# setup some constants
	#resting_roi_names = ['LeftInsula','RightInsula','LeftAmygdala',
	#		     'RightAmygdala','LeftVS','RightVS','LeftVLPFC','RightVLPFC',
	#	       'BilateralAmygdala']
	#resting_roi_images = [conf.ROI_L_insula,conf.ROI_R_insula,conf.ROI_L_amyg,conf.ROI_R_amyg,
	#		conf.ROI_VS_L,conf.ROI_VS_R,conf.ROI_L_VLPFC,conf.ROI_R_VLPFC,
	#		conf.ROI_amygdala_LR]	
 
	# setup some constants (narrowed down, since not all are currently available in 91x109x91 dimensions)
	resting_roi_names = ['LeftAmygdala','RightAmygdala','LeftVLPFC','RightVLPFC','LeftAI','RightAI','LeftDLPFC','RightDLPFC','vmPFC']
	resting_roi_images = [conf.ROI_L_amyg,conf.ROI_R_amyg,conf.ROI_L_VLPFC,conf.ROI_R_VLPFC,conf.ROI_L_ant_insula, 		
			      conf.ROI_R_ant_insula, conf.ROI_L_DLPFC, conf.ROI_R_DLPFC, conf.ROI_BA10]	
 
	ds = datasource(directory,sequence)
	pp = preprocess_mni(conf,useFieldmap)
	
	nu = pe.Node(interface=wrap.Nuisance(), name="nuisance")
	nu.inputs.white_mask = conf.ROI_white
	nu.inputs.time_repetition = conf.time_repetition

	column_select = pe.Node(interface=wrap.ColumnSelect(),name="column_select")
	column_select.inputs.selection = "18,24"
	column_select.inputs.complement = True


	glm = pe.Node(interface=fsl.GLM(), name="glm")
	glm.inputs.out_res_name = "residual.4d.nii"
	
	glm_NGS = pe.Node(interface=fsl.GLM(), name="glm_NGS")
	glm_NGS.inputs.out_res_name = "residual.4d.nii"

	filt = pe.Node(interface=fsl.ImageMaths(), name="filter")
	filt.inputs.op_string = conf.filter_image_bptf
	filt.inputs.terminal_output = 'none'
	
	alff = dict()
	alff_nm = []
	for hl in [[0.01,0.1],[0.01,0.027],[0.027,0.073]]:
		nm = "ALFF"+str(hl[1]).replace("0.","_")
		alff_nm.append(nm)		
		alff[nm] = CPAC.alff.create_alff(wf_name=nm.lower())
		alff[nm].inputs.hp_input.hp = hl[0]
		alff[nm].inputs.lp_input.lp = hl[1]
	
	reho = CPAC.reho.create_reho()
	reho.inputs.inputspec.cluster_size = 27
	
	nc = CPAC.network_centrality.create_resting_state_graphs(wf_name='network_centrality')
	nc.inputs.inputspec.method_option=0
	nc.inputs.inputspec.weight_options=[True, True]	
	nc.inputs.inputspec.threshold_option = 1
	nc.inputs.inputspec.threshold = 0.0744 
	nc.inputs.inputspec.template = conf.OASIS_labels
	zscore =  CPAC.network_centrality.get_zscore(wf_name='z_score')

	sca = dict()
	maskave = dict()
	gunzip = dict()
	
	for mask in ["LeftAmygdala","RightAmygdala","LeftVLPFC","RightVLPFC","LeftAI","RightAI","LeftDLPFC","RightDLPFC","vmPFC"]:
		sca[mask] = CPAC.sca.create_sca(name_sca="sca_"+mask);
		maskave[mask] = pe.Node(interface=afni.Maskave(),name="roi_ave_"+mask)
		maskave[mask].inputs.outputtype = "NIFTI"
		maskave[mask].inputs.quiet= True
		maskave[mask].inputs.mask = resting_roi_images[resting_roi_names.index(mask)]
		gunzip[mask] = pe.Node(interface=misc.Gunzip(),name="gunzip_"+mask)
	
	roiave = pe.MapNode(interface=afni.Maskave(),name="roi_ave",iterfield="mask")
	roiave.inputs.outputtype = "NIFTI"
	roiave.inputs.mask = resting_roi_images
	roiave.inputs.quiet= True

	corroi = pe.Node(interface=wrap.CorrelateROIs(),name="corr_roi")
	corroi.inputs.roi_names = resting_roi_names
	corroi.inputs.task_name = "Resting_State"
	corroi.inputs.out_file = subject+"_"+sequence+"_outcomes_CORR.csv"

	
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	

	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	
	if useFieldmap:	
		task.connect([(ds,pp,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
	else:
		task.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])	
	task.connect([(pp,nu,[('output.ufunc','source'),
				 ('output.mask','brain_mask'),
				 ('output.movement','movement')])])
	


	task.connect(nu,"regressors",glm,"design")
	task.connect(pp,"output.func",glm,"in_file")
	task.connect(glm,"out_res",filt,"in_file")
	
	task.connect(filt,"out_file",roiave,"in_file")
	task.connect(roiave,"out_file",corroi,"in_files")
		
	
	for nm in alff_nm:
		task.connect(glm,'out_res',alff[nm],'inputspec.rest_res')
		task.connect(pp,'output.mask',alff[nm],'inputspec.rest_mask')	

	task.connect(filt,"out_file",reho,"inputspec.rest_res_filt")
	task.connect(pp,"output.mask",reho,"inputspec.rest_mask")
	
	task.connect(glm,'out_res',nc,'inputspec.subject')
	task.connect(nc,'outputspec.centrality_outputs',zscore,'inputspec.input_file')
	task.connect(pp,'output.mask',zscore,'inputspec.mask_file')

	for mask in ["LeftAmygdala","RightAmygdala","LeftVLPFC","RightVLPFC","LeftAI","RightAI","LeftDLPFC","RightDLPFC","vmPFC"]:
		task.connect(filt,"out_file",maskave[mask],"in_file")
		task.connect(filt,"out_file",sca[mask],"inputspec.functional_file")
		task.connect(maskave[mask],"out_file",sca[mask],"inputspec.timeseries_one_d")
		task.connect(sca[mask],("outputspec.Z_score",subset,0),gunzip[mask],'in_file')
		task.connect(sca[mask],"outputspec.Z_score",datasink,"data.sca."+mask)
	
	
	task.connect(reho,"outputspec.z_score",datasink,"data.reho")
	task.connect(corroi,"out_file",datasink,"csv.@par5")
	# alff_Z_img in 0.3.5 now in 0.3.6 falff_img	
	for nm in alff_nm:	
		task.connect(alff[nm],"outputspec.alff_Z_img",datasink,"data."+nm.lower())
	
	task.connect(zscore,"outputspec.z_score_img",datasink,"data.nc")
	

	# print and save the output of the preprocess pipeline
	gold.save_files(task,pp.get_node('output'),datasink,["func","movement","struct","mask"], not noPrint)	
	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task


# run pipeline if used as standalone script
if __name__ == "__main__":	
	opts = "[-fieldmap|-noprint]"
	opt_list = []	

	# get arguments
	if len(sys.argv) < 2:
		print "Usage: resting "+opts+" <subject directory>"
		sys.exit(1)
	
	# logging verbosity
	import time
	import nipype
	import logging
	from nipype import config
	import nipype.interfaces.matlab as mlab 
	
	# pick dataset that we'll be wroking on
	for arg in sys.argv:
		 if arg in opts:
		 	opt_list.append(arg)
	directory = sys.argv[len(sys.argv)-1]
	
	# check directory
	if not os.path.exists(directory):
		print "Error: data directory "+directory+" does not exist"
		sys.exit(1)
	
	# setup logging, display and other config
	disp = os.environ['DISPLAY']
	log_dir = os.path.abspath(directory)+"/logs"
	if not os.path.exists(log_dir):
		os.mkdir(log_dir)
	cfg = dict(logging={'interface_level':'INFO',
						'workflow_level':'INFO',
						'log_to_file': True,
						'log_directory': log_dir},
    		   execution={'stop_on_first_crash': True,
    		   			  'display_variable':disp,
                      	  'hash_method': 'timestamp',
                      	  'keep_inputs': True,
                      	  'remove_unnecessary_outputs': False})
	config.update_config(cfg)
	nipype.logging.update_logging(config)
	log = nipype.logging.getLogger('workflow')
	l = nipype.logging.getLogger('interface').parent.handlers[0]
	l.setFormatter(logging.Formatter('%(name)-2s %(levelname)-2s:\t %(message)s'))
	###########
	
	# setup matlab env
	mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
	mlab.MatlabCommand.set_default_terminal_output('stream')
	#mlab.MatlabCommand.set_default_paths(bin_dir)
	
	
	if "-fieldmap" in opt_list:	
		useFieldmap = True
		opt_list.remove("-fieldmap")
	if "-noprint" in opt_list:	
		noPrint = True
		opt_list.remove("-noprint")

	log.info("\n\nRESTING pipeline ...\n\n")
	t = time.time()		
	resting1 = resting(directory,'lams_resting_state')
	resting1.run()
	log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
			
	log.info("\n\npipeline complete\n\n")
