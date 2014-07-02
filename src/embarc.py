#!/usr/bin/env python
# EMBARC 2.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import sys
import os                                  
import re


# some predefined constants
bin_dir = os.path.dirname(os.path.realpath(__file__))
OASIS_template = bin_dir+"/../validation/templates/OASIS-30_Atropos_template_in_MNI152_2mm.nii.gz"
OASIS_labels = bin_dir+"/../validation/templates/OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_2mm.nii.gz"
ROI_dir = bin_dir+"/../validation/ROI_EMBARC_scaled2mm/"

"""
EMBARC 1.5 PreProcessing Pipeline
input: 
	func  - functional images
	struct - structural image
output:
	func  - functional processed images
	mask  - mask image from the functional
	movement - realign movement parameters
"""
def preprocess():
	import nipype.interfaces.spm as spm          # spm
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.utility as util     # utility
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.afni as afni		 # afni
	
	fsl.FSLCommand.set_default_output_type('NIFTI')
	
	#d = os.path.dirname(os.path.realpath(__file__))+"/../validation/templates/"
	#template = d+"OASIS-30_Atropos_template_in_MNI152_2mm.nii.gz"
	template = OASIS_template
	
	preproc = pe.Workflow(name='preprocess')

	inputnode = pe.Node(interface=util.IdentityInterface(
				fields=['func','struct']),name='input')

	realign = pe.Node(interface=spm.Realign(), name="realign")
	realign.inputs.register_to_mean = True

	func_bet = pe.Node(interface=fsl.BET(), name="bet_mean")
	func_bet.inputs.mask = True
	func_bet.inputs.frac = 0.5
	func_bet.inputs.robust = True
	func_bet.inputs.vertical_gradient = 0
	
	sfunc_bet = pe.Node(interface=fsl.BET(), name="bet_func")
	sfunc_bet.inputs.mask = True
	sfunc_bet.inputs.frac = 0.5
	sfunc_bet.inputs.robust = True
	sfunc_bet.inputs.vertical_gradient = 0
	
	#struct_bet = pe.Node(interface=fsl.BET(), name="struct_bet")
	#struct_bet.inputs.mask = True
	#struct_bet.inputs.frac = 0.5
	#struct_bet.inputs.robust = True
	#struct_bet.inputs.vertical_gradient = 0
	#struct_bet.inputs.center = [133,115,88]	
	
	struct_bet = pe.Node(interface=afni.SkullStrip(), name="skull_strip")
	struct_bet.inputs.args = " -ld 30 -push_to_edge -no_avoid_eyes "
	struct_bet.inputs.outputtype = 'NIFTI'
	
	#coregister = pe.Node(interface=spm.Coregister(), name="coregister")
	#coregister.inputs.jobtype = 'estimate'
	
	# coregister structural to mean functional 
	coregister = pe.Node(interface=fsl.FLIRT(), name='flirt_mean2struct')
	coregister.inputs.cost = 'corratio'
	coregister.inputs.bins = 256
	coregister.inputs.dof = 12
	coregister.inputs.interp = 'trilinear'
	coregister.inputs.searchr_x = [-180, 180]
	coregister.inputs.searchr_y = [-180, 180]
	coregister.inputs.searchr_z = [-180, 180]
	
	# apply transform to realigned 4d functional referenced to structural
	coreg_xfm = pe.Node(interface=fsl.ApplyXfm(),name='apply_func2struct')
	coreg_xfm.inputs.interp = 'trilinear'
	
	# coregister structural to template	            
	flirt = pe.Node(interface=fsl.FLIRT(), name='flirt_struct2template')
	flirt.inputs.cost = 'corratio'
	flirt.inputs.bins = 256
	flirt.inputs.dof = 12
	flirt.inputs.interp = 'trilinear'
	flirt.inputs.searchr_x = [-180, 180]
	flirt.inputs.searchr_y = [-180, 180]
	flirt.inputs.searchr_z = [-180, 180]
	flirt.inputs.reference = template
	
	# apply transform to coregistered 4D functional referenced template
	xfm = pe.Node(interface=fsl.ApplyXfm(),name='apply_func2template')
	xfm.inputs.interp = 'trilinear'
	xfm.inputs.reference = template
		
	# apply transform to coregistered 4D functional referenced template
	sxfm = pe.Node(interface=fsl.ApplyXfm(),name='apply_struct2template')
	sxfm.inputs.interp = 'trilinear'
	sxfm.inputs.reference = template


	# do fnirt
	fnirt_s2t = pe.Node(interface=fsl.FNIRT(), name='fnirt_struct2template')	
	fnirt_s2t.inputs.ref_file = template
	
	warp_f2t  = pe.Node(interface=fsl.ApplyWarp(), name='warp_func2template')	
	warp_f2t.inputs.ref_file = template
	
	warp_s2t  = pe.Node(interface=fsl.ApplyWarp(), name='warp_structtemplate')	
	warp_s2t.inputs.ref_file = template
	
	susan = pe.Node(interface=fsl.SUSAN(), name="smooth")
	susan.inputs.brightness_threshold = 2000.0
	susan.inputs.fwhm = 8

	despike = pe.Node(interface=afni.Despike(), name='despike')
	despike.inputs.outputtype = 'NIFTI'

	outputnode = pe.Node(interface=util.IdentityInterface(
				 fields=['func','ufunc','mask','movement','struct']),name='output')


	preproc.connect([(inputnode,realign,[('func','in_files')]),
                     (realign,func_bet,[('mean_image','in_file')]), 
					 (inputnode,struct_bet,[('struct','in_file')]),
					 
					 (struct_bet,coregister,[('out_file','reference')]),
					 (func_bet,coregister,[('out_file','in_file')]),
				
					 (coregister,coreg_xfm,[('out_matrix_file','in_matrix_file')]),
					 (realign,coreg_xfm,[('realigned_files','in_file')]),
					 (struct_bet,coreg_xfm,[('out_file','reference')]),
				
					 (struct_bet,flirt,[('out_file','in_file')]),
		             (coreg_xfm,xfm,[('out_file','in_file')]),
		             (flirt,xfm,[('out_matrix_file','in_matrix_file')]),
		           	 (struct_bet,sxfm,[('out_file','in_file')]),
		             (flirt,sxfm,[('out_matrix_file','in_matrix_file')]),
		           
					 (sxfm,fnirt_s2t,[('out_file','in_file')]),
		             (xfm,warp_f2t,[('out_file','in_file')]),
		             (fnirt_s2t,warp_f2t,[('field_file','field_file')]),
		             
		             (sxfm,warp_s2t,[('out_file','in_file')]),
		             (fnirt_s2t,warp_s2t,[('field_file','field_file')]),
		             
		             #(xfm,despike,[('out_file','in_file')]),
		           	 (warp_f2t,despike,[('out_file','in_file')]),
		             (despike,susan, [('out_file', 'in_file')]), 
		             (despike,outputnode, [('out_file', 'ufunc')]), 
		             (susan,outputnode,[('smoothed_file','func')]),
		             (susan,sfunc_bet,[('smoothed_file','in_file')]),
		             (sfunc_bet,outputnode,[('mask_file','mask')]),
		             (realign,outputnode,[('realignment_parameters','movement')]),
		             (warp_s2t,outputnode,[('out_file','struct')])
		             ])
	return preproc

# convert eprime to design matrix
def eprime2dm(eprime, pppi):
   	import os
   	import re
   	import scipy.io as sp
	import glob as gl
	import numpy
	import nipype.interfaces.matlab as mlab 
	from nipype.interfaces.base import Bunch
    
	m = re.search("eprime_([a-z]+)\.txt",eprime)
	sequence = m.group(1)
    
	mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
	
	# execute matlab script
	m = mlab.MatlabCommand()
	m.inputs.mfile = False
	m.inputs.script = sequence+"_eprime2dm(\'"+eprime+"\');"
	m.run();

	# extract variables from .mat files
	mat = os.path.join(os.path.dirname(eprime),'nDM*.mat')
	mat = gl.glob(mat)
	if len(mat) > 0:
		mat = mat[0]
	dm = sp.loadmat(mat,squeeze_me=True)
		
	names  = []
	for a in dm.get('names').tolist():
		names.append(str(a))
	onsets = []  
	for a in dm.get('onsets').tolist():
		if type(a) == numpy.ndarray:
			onsets.append(a.tolist())
		else:
			onsets.append([a])
	durations = []
	for a in dm.get('durations').tolist():
		if type(a) == numpy.ndarray:
			durations.append(a.tolist())
		else:
			durations.append([a])
	# load up values and convert them
	# for PPPI remove last column for reward PPI; for ert last 3 columns
	# error, posterror, misc
	if pppi:
		trim = 3
		if sequence == "reward":
			trim = 1
		names = names[0:(len(names)-trim)]
		durations = durations[0:(len(durations)-trim)]
		onsets = onsets[0:(len(onsets)-trim)]
	
	# create bunch to return
	bunch = Bunch(conditions=names,onsets=onsets,durations=durations)
	if 'pmod' in dm:
		pmod = []
		for i in range(0,len(dm.get('pmod'))):
			if len(dm['pmod']['name'][i]) == 0:
				pmod.append(None)
			else:
				name = str(dm['pmod']['name'][i])
				param = dm['pmod']['param'][i].tolist()
				poly = dm['pmod']['poly'][i]
				pmod.append(Bunch(name=[name],param=[param],poly=[poly]))
		bunch.pmod = pmod
	
	return bunch

# predefined ert contrast estimates
def get_contrasts(eprime):	
	import re
	m = re.search("eprime_([a-z]+)\.txt",eprime)
	sequence = m.group(1)
	
	# extract sequence
	if sequence == "ert":
		# ERT contrasts
		cont1 = ('iI_cI','T', ['cI', 'iI'],[-1, 1])
		cont2 = ('Conflict','T',['cC','cI','iC','iI'],[-1, 1,-1,1])
		return [cont1,cont2]
	elif sequence == "reward":
		# Reward contrasts
		cont1 = ('RewardExpectancy','T', ['anticipationxanti^1'],[1])
		cont2 = ('PredictionError','T', ['outcomexsignedPE^1'],[1])
		return [cont1,cont2]
	return []


"""
EMBARC 1.0 Level1 Analysis Pipeline
"""
def level1analysis(pppi):
	import nipype.interfaces.spm as spm          # spm
	import nipype.interfaces.matlab as mlab      # how to run matlab
	import nipype.interfaces.utility as util     # utility
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.algorithms.modelgen as model   # model specification
	
	mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
	name = 'level1'
	if pppi:
		name = "level1_pppi"
	
	
	l1analysis = pe.Workflow(name=name)
	inputnode = pe.Node(interface=util.IdentityInterface(
						fields=['movement','func','struct','mask','behav']),
						name='input')


	modelspec = pe.Node(interface=model.SpecifySPMModel(), name= "modelspec")
	modelspec.inputs.concatenate_runs   = True
	modelspec.inputs.time_repetition = 2
	modelspec.inputs.high_pass_filter_cutoff = 60
	modelspec.inputs.input_units = 'secs'
	
	level1design = pe.Node(interface=spm.Level1Design(), name= "level1design")
	level1design.inputs.bases = {'hrf':{'derivs': [0,0]}}
	level1design.inputs.timing_units = 'secs'
	level1design.inputs.interscan_interval = 2

	level1estimate = pe.Node(interface=spm.EstimateModel(), name="level1estimate")
	level1estimate.inputs.estimation_method = {'Classical' : 1}
	
	# no need for contrast for pppi model
	contrastestimate = pe.Node(interface = spm.EstimateContrast(), name="contrastestimate")
	
	outputnode = pe.Node(interface=util.IdentityInterface(
			fields=['contrasts','con_images','spmT_images']),name='output')
	
	
	l1analysis.connect([(inputnode,modelspec,[('movement','realignment_parameters'),
					  						('func','functional_runs'),
					  						(('behav', eprime2dm, pppi),'subject_info')]),	
					  (modelspec,level1design,[('session_info','session_info')])])
	if pppi:
		l1analysis.connect(level1design,'spm_mat_file',outputnode,'contrasts')
	else:	            
		l1analysis.connect([
			(level1design,level1estimate,[('spm_mat_file','spm_mat_file')]),
			(inputnode,contrastestimate,[(('behav',get_contrasts),'contrasts')]),
		    (level1estimate,contrastestimate,[('spm_mat_file','spm_mat_file'),
		                                      ('beta_images','beta_images'),
		                                      ('residual_image','residual_image')]),
		    (contrastestimate,outputnode,[('spm_mat_file','contrasts'),
		              		              			('con_images','con_images'),
		              		              			('spmT_images','spmT_images')])])
	return l1analysis		             

"""
EMBARC get subject name from a directory
"""
def get_subject(directory):
	m = re.search('embarc_CU_([A-Z0-9]+_\dR\d)_mri',directory)
	if m:
		return m.group(1)
	return "subject"

"""
EMBARC 1.0 Input Data Source
"""
def datasource(directory, sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o

	# define some variables beforehand
	subject=get_subject(directory)
	
	# func template
	func_template = 'dicom_bold_%s/%s_bold_%s.nii'
	if sequence == 'asl':
		func_template = 'dicom_%s/%s_%s_*.nii'

	# define templates for datasource
	field_template = dict(func=func_template,
                          struct='dicom_anatomical/%s_anatomical.nii')
	template_args  = dict(func=[['sequence','subject_id','sequence']],
    					  struct=[['subject_id']])                

	# add behavior file to task oriented design
	if sequence == 'ert' or sequence == 'reward':
		field_template['behav'] = 'dicom_bold_%s/eprime_%s.txt'
		template_args['behav']  = [['sequence','sequence']]

	# specify input dataset just pass through parameters
	datasource = pe.Node(interface=nio.DataGrabber(
						 infields=['subject_id','sequence'], 
						 outfields=['func', 'struct','behav']),
	                     name = 'datasource')
	datasource.inputs.base_directory = os.path.abspath(directory)
	datasource.inputs.template = '*'
	datasource.inputs.field_template = field_template
	datasource.inputs.template_args  = template_args
	datasource.inputs.subject_id = subject
	datasource.inputs.sequence = sequence
	datasource.inputs.sort_filelist = True

	return datasource

def subset(x,i):
	return x[i]	

"""
EMBARC 1.0 Task Sequence Ex: ert/reward
"""
def task(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	subject = get_subject(directory)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# get components
	ds = datasource(directory,sequence)
	pp = preprocess()
	l1 = level1analysis(False);
	l2 = level1analysis(True);
		
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])
	task.connect([(pp,l1,[('output.func','input.func'),
						 ('output.mask','input.mask'),
						 ('output.movement','input.movement')])])
	task.connect([(pp,l2,[('output.func','input.func'),
						 ('output.mask','input.mask'),
						 ('output.movement','input.movement')])])
	task.connect([(ds,l1,[('behav','input.behav')])])
	task.connect([(ds,l2,[('behav','input.behav')])])
	
	# define ROIs from  level1 analysis
	rois = []
	l1_rois = []
	l1_names = []
	task_name = ""
	task_units = ""
	
	# ERT ROIs
	if sequence == 'ert':
		task_name = "EmotionRecognitionTask"
		task_units = "ParameterEstimate"
		l1_names = ["iI_cI","conflict"]
		# 1st level ROIs
		l1_rois = [[("LeftAmygdala_iI_minus_cI",ROI_dir+"L_amyg_2mm.img"),
					("RightAmygdala_iI_minus_cI",ROI_dir+"R_amyg_2mm.img"),
					("LeftInsula_iI_minus_cI",ROI_dir+"L_insula_2mm.img"),
					("RightInsula_iI_minus_cI",ROI_dir+"R_insula_2mm.img"),
					("pgCing_iI_minus_cI",ROI_dir+"pgACC_2mm.img"),
					("dCing_iI_minus_cI",ROI_dir+"dACC_2mm.img"),
					("sgCing_iI_minus_cI",ROI_dir+"sgACC_2mm.img"),
					("LeftMFG_iI_minus_cI",ROI_dir+"L_MFG_compensate_2mm.img"),
					("RightMFG_iI_minus_cI",ROI_dir+"R_MFG_compensate_2mm.img")],
					
				   [("LeftAmyg_I_minus_C",ROI_dir+"L_amyg_2mm.img"),
					("RightAmyg_I_minus_C",ROI_dir+"R_amyg_2mm.img"),
					("LeftInsula_I_minus_C",ROI_dir+"L_insula_2mm.img"),
					("RightInsula_I_minus_C",ROI_dir+"R_insula_2mm.img"),
					("pgCing_I_minus_C",ROI_dir+"pgACC_2mm.img"),
					("dCing_I_minus_C",ROI_dir+"dACC_2mm.img"),
					("sgCing_I_minus_C",ROI_dir+"sgACC_2mm.img")]]
		# ROIs for PPI
		rois = [("ER_pgACC",ROI_dir+"pgACC_2mm.img",
					[("LeftAmyg_iI_minus_cI_pgAcc_PPI",ROI_dir+"L_amyg_2mm.img"),
					 ("RightAmyg_iI_minus_cI_pgAcc_PPI",ROI_dir+"R_amyg_2mm.img")]),
				("ER_L_amyg",ROI_dir+"L_amyg_2mm.img",
					[("pgCing_iI_minus_cI_LAmyg_PPI",ROI_dir+"pgACC_2mm.img"),
					 ("dCing_iI_minus_cI_LAmyg_PPI",ROI_dir+"dACC_2mm.img"),
					 ("sgCing_iI_minus_cI_LAmyg_PPI",ROI_dir+"sgACC_2mm.img")]),
				("ER_R_amyg",ROI_dir+"R_amyg_2mm.img",
					[("pgCing_iI_minus_cI_RAmyg_PPI",ROI_dir+"pgACC_2mm.img"),
					 ("dCing_iI_minus_cI_RAmyg_PPI",ROI_dir+"dACC_2mm.img"),
					 ("sgCing_iI_minus_cI_RAmyg_PPI",ROI_dir+"sgACC_2mm.img")])]				

	elif sequence == 'reward':
		task_name = "RewardTask"
		task_units = "ParameterEstimate"
		l1_names = ["anticipation","outcome"]
		l1_rois = [[("BA10_anticipation",ROI_dir+"BA10_5mm.img"),
					("LeftBA9_anticipation",ROI_dir+"BA09_left_2mm.img"),
					("RightBA9_anticipation",ROI_dir+"BA09_right_2mm.img"),
					("LeftVS_anticipation",ROI_dir+"VS_left_2mm.img"),
					("RightVS_anticipation",ROI_dir+"VS_right_2mm.img"),
					("LeftBA47_anticipation",ROI_dir+"Left_VLPFC_5mm.img"),
					("RightBA47_anticipation",ROI_dir+"Right_VLPFC_5mm.img"),
					("BR1_anticipation",ROI_dir+"BR1_2mm.img"),
					("BR2_anticipation",ROI_dir+"BR2_2mm.img"),
					("BR3_anticipation",ROI_dir+"BR3_2mm.img"),
					("BR4_anticipation",ROI_dir+"BR4_2mm.img")],
					
				   [("LeftBA9_outcome",ROI_dir+"BA09_left_2mm.img"),
					("RightBA9_outcome",ROI_dir+"BA09_right_2mm.img"),
					("LeftVS_outcome",ROI_dir+"VS_left_2mm.img"),
					("RightVS_outcome",ROI_dir+"VS_right_2mm.img"),
					("BR1_outcome",ROI_dir+"BR1_2mm.img"),
					("BR2_outcome",ROI_dir+"BR2_2mm.img"),
					("BR3_outcome",ROI_dir+"BR3_2mm.img"),
					("BR4_outcome",ROI_dir+"BR4_2mm.img")]]
		# ROIs for PPI
		rois = [("Reward_VS",ROI_dir+"VS_L&R.nii",
					[("BR1_anticipation_PPI",ROI_dir+"BR1_2mm.img"),
					 ("BR2_anticipation_PPI",ROI_dir+"BR2_2mm.img"),
					 ("BR3_anticipation_PPI",ROI_dir+"BR3_2mm.img"),
					 ("BR4_anticipation_PPI",ROI_dir+"BR4_2mm.img")])]
			
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	task.connect(l1,"output.con_images",datasink,"data.l1")
	
	for i in range(0,len(l1_rois)):
		# extract ROIs from  level1 analysis
		extract = pe.Node(interface=wrap.ROIExtractor(), name="extract_level1_rois_"+l1_names[i])
		extract.inputs.roi_images = list(zip(*l1_rois[i])[1])
		extract.inputs.average = 'none'
		extract.inputs.interpelation = 0
		task.connect(l1,("output.con_images",subset,i),extract,'source')
		
		# save CSV file
		save_csv = pe.Node(name="save_csv_"+l1_names[i],
			interface=Function(input_names=["task","units","names","ext","output"],
			output_names=["csv_file"],function=wrap.save_csv))
		save_csv.inputs.task = task_name
		save_csv.inputs.units = task_units
		save_csv.inputs.names = list(zip(*l1_rois[i])[0])
		save_csv.inputs.output = subject+"_"+sequence+"_outcomes_"+l1_names[i]+".csv"
		task.connect(extract,"mat_file",save_csv,"ext")
		task.connect(save_csv,"csv_file",datasink,"csv.@par"+l1_names[i])
		
		
	# do PPPI on a set of ROIs
	
	# now do gPPI analysis
	for roi in rois:
		pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
		pppi.inputs.voi_name = roi[0]
		pppi.inputs.voi_file = roi[1]
		pppi.inputs.subject = subject
		
		estimate = pe.Node(interface=spm.EstimateModel(), name="estimate_"+roi[0])
		estimate.inputs.estimation_method = {'Classical' : 1}

		contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
		
		extract = pe.Node(interface=wrap.ROIExtractor(), name="extract_"+roi[0])
		extract.inputs.roi_images = list(zip(*roi[2])[1])
		extract.inputs.average = 'none'
		extract.inputs.interpelation = 0	
	
		# save CSV file
		save_csv = pe.Node(name="save_csv_"+roi[0],
			interface=Function(input_names=["task","units","names","ext","output"],
			output_names=["csv_file"],function=wrap.save_csv))
		save_csv.inputs.task = task_name
		save_csv.inputs.units = task_units
		save_csv.inputs.names = list(zip(*roi[2])[0])
		save_csv.inputs.output = subject+"_"+sequence+"_outcomes_"+roi[0]+".csv"
		
		task.connect(l2,'output.contrasts',pppi,'spm_mat_file')
		task.connect(pppi,'spm_mat_file',estimate,'spm_mat_file')
		task.connect(ds,('behav',get_contrasts),contrast,'contrasts') # diff contrasts
		task.connect(estimate,'spm_mat_file',contrast,'spm_mat_file')
		task.connect(estimate,'beta_images',contrast,'beta_images')
		task.connect(estimate,'residual_image',contrast,'residual_image')
		task.connect(contrast,('con_images',subset,0),extract,'source')
		task.connect(extract,"mat_file",save_csv,"ext")
		task.connect(save_csv,"csv_file",datasink,"csv.@par"+roi[0])
		task.connect(contrast,('con_images',subset,0),datasink,"data.l1_pppi_"+roi[0])
	
	task.connect(pp,"output.func",datasink,"data.functional")
	task.connect(pp,"output.struct",datasink,"data.structural")
	task.connect(pp,"output.mask",datasink,"data.mask")
	
	# print
	p_ru = pe.Node(interface=wrap.Print(), name="print_realign_params")
	p_ru.inputs.out_file = "realignment_parameters.ps"
	
	p_dm = pe.Node(interface=wrap.Print(), name="print_design_matrix")
	p_dm.inputs.out_file = "level1_design_matrix.ps"
	
	p_pdm = pe.Node(interface=wrap.Print(), name="print_pppi_design_matrix")
	p_pdm.inputs.out_file = "level1_pppi_design_matrix.ps"
	
	p_struct = pe.Node(interface=wrap.Print(), name="print_anatomical")
	p_struct.inputs.out_file = "brain_anatomical.ps"
	
	p_func = pe.Node(interface=wrap.Print(), name="print_func")
	p_func.inputs.out_file = "brain_func.ps"
	
	p_mask = pe.Node(interface=wrap.Print(), name="print_mask")
	p_mask.inputs.out_file = "brain_mask.ps"
	
	
	task.connect(pp,"output.movement",p_ru,"in_file")
	task.connect(l1,"output.contrasts",p_dm,"in_file")
	task.connect(l2,"output.contrasts",p_pdm,"in_file")
	task.connect(pp,"output.struct",p_struct,"in_file")
	task.connect(pp,"output.func",p_func,"in_file")
	task.connect(pp,"output.mask",p_mask,"in_file")
	
	task.connect(p_ru,"out_file",datasink,"ps.@par1")
	task.connect(p_dm,"out_file",datasink,"ps.@par2")
	task.connect(p_pdm,"out_file",datasink,"ps.@par3")
	task.connect(p_func,"out_file",datasink,"ps.@par4")
	task.connect(p_struct,"out_file",datasink,"ps.@par5")
	task.connect(p_mask,"out_file",datasink,"ps.@par6")
	
		
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task


"""
EMBARC 1.0 Resting Sequence Ex: resting1/resting2
"""
def resting(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import CPAC									 # import CPAC nuisance
	import nipype.interfaces.fsl as fsl          # fsl
	import wrappers as wrap	
	import nipype.interfaces.afni as afni		 # afni
	from nipype.interfaces.utility import Function
	import nipype.interfaces.io as nio           # Data i/o
	import nipype.algorithms.misc as misc
		
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	subject = get_subject(directory)
	
	# get components
	ds = datasource(directory,sequence)
	pp = preprocess()
	nu = pe.Node(interface=wrap.Nuisance(), name="nuisance")
	nu.inputs.white_mask = ROI_dir+"white_test3.nii"
	glm = pe.Node(interface=fsl.GLM(), name="glm")
	glm.inputs.out_res_name = "residual.4d.nii"
	filt = pe.Node(interface=fsl.ImageMaths(), name="filter")
	filt.inputs.op_string = ' -bptf 128 12.5 '
	filt.inputs.terminal_output = 'none'
	
	alff = CPAC.alff.create_alff(wf_name='alff')
	alff.inputs.hp_input.hp = 0.01
	alff.inputs.lp_input.lp = 0.08
	
	reho = CPAC.reho.create_reho()
	reho.inputs.inputspec.cluster_size = 27
	
	nc = CPAC.network_centrality.create_resting_state_graphs(wf_name='network_centrality')
	nc.inputs.centrality_options.method_options=[True, True]
	nc.inputs.centrality_options.weight_options=[True, True]
	nc.inputs.inputspec.threshold_option = 1
	nc.inputs.inputspec.threshold = 0.0744 
	nc.inputs.inputspec.template = OASIS_labels
	zscore =  CPAC.network_centrality.get_zscore(wf_name='z_score')
	
	
	sca = CPAC.sca.create_sca(name_sca='sca');
	#tonifti = pe.Node(interface=misc.CreateNifti(),name="to_nifti")
	#tonifti.inputs.data_file = ROI_dir+'L_insula_2mm.img'
	#tonifti.inputs.header_file = ROI_dir+'L_insula_2mm.hdr'
	
	maskave = pe.Node(interface=afni.Maskave(),name="mask_ave")
	maskave.inputs.outputtype = "NIFTI"
	maskave.inputs.quiet= True
	gunzip = pe.Node(interface=misc.Gunzip(),name="gunzip")
	
	#maskave.inputs.mask = ROI_dir+'L_insula_2mm.img'
	
	task_name = "Resting_First_Block"
	if sequence == 'resting2':
		task_name = "Resting_Second_Block"
		
	task_units = "Z_Score"
	resting_roi_names = ['LeftInsula','RightInsula','LeftAmygdala',
						'RightAmygdala','LeftVS','RightVS','LeftBA9','RightBA9',
						'BR1','BR2','BR3','BR4','BR9']
	resting_roi_images = [ROI_dir+'L_insula_2mm.img',ROI_dir+'R_insula_2mm.img',
					ROI_dir+'L_amyg_2mm.img',ROI_dir+'R_amyg_2mm.img',
					ROI_dir+'VS_left_2mm.img',ROI_dir+'VS_right_2mm.img',
					ROI_dir+'BA09_left_2mm.img',ROI_dir+'BA09_right_2mm.img',
					ROI_dir+'BR1_2mm.img',ROI_dir+'BR2_2mm.img',
					ROI_dir+'BR3_2mm.img',ROI_dir+'BR4_2mm.img',
					ROI_dir+'BR9_2mm.img' ]
	

	extract = dict()
	save_csv = dict()
	
	for nm in ["ALFF","ReHO","NC","SCA"]: 
		ext = pe.Node(interface=wrap.ROIExtractor(),name="extract_"+nm)
		ext.inputs.roi_images = resting_roi_images
		ext.inputs.average = 'none'
		ext.inputs.interpelation = 0
		extract[nm] = ext
		
		csv = pe.Node(name="save_"+nm,interface=Function(
			input_names=["task","units","names","ext","output"],
			output_names=["csv_file"],function=wrap.save_csv))
		csv.inputs.task = task_name
		csv.inputs.units = task_units
		csv.inputs.names = [nm+"_"+ s for s in resting_roi_names]
		csv.inputs.output = subject+"_"+sequence+"_outcomes_"+nm+".csv"
		save_csv[nm] = csv
	
	
	
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	

	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	
	task.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])
	task.connect([(pp,nu,[('output.ufunc','source'),
						 ('output.mask','brain_mask'),
						 ('output.movement','movement')])])
	task.connect(nu,"regressors",glm,"design")
	task.connect(pp,"output.func",glm,"in_file")
	task.connect(glm,"out_res",filt,"in_file")
	
	task.connect(glm,'out_res',alff,'inputspec.rest_res')
	task.connect(pp,'output.mask',alff,'inputspec.rest_mask')	

	task.connect(filt,"out_file",reho,"inputspec.rest_res_filt")
	task.connect(pp,"output.mask",reho,"inputspec.rest_mask")
	
	task.connect(glm,'out_res',nc,'inputspec.subject')
	task.connect(nc,'outputspec.centrality_outputs',zscore,'inputspec.input_file')
	task.connect(pp,'output.mask',zscore,'inputspec.mask_file')

	task.connect(filt,"out_file",maskave,"in_file")
	task.connect(pp,"output.mask",maskave,"mask")
	#task.connect(tonifti,"nifti_file",maskave,"mask")
	
	task.connect(filt,"out_file",sca,"inputspec.functional_file")
	task.connect(maskave,"out_file",sca,"inputspec.timeseries_one_d")
	
	task.connect(alff,"outputspec.alff_Z_img",extract["ALFF"],'source')
	task.connect(extract["ALFF"],"mat_file",save_csv["ALFF"],"ext")

	task.connect(reho,"outputspec.z_score",extract["ReHO"],'source')
	task.connect(extract["ReHO"],"mat_file",save_csv["ReHO"],"ext")
	
	task.connect(zscore,("outputspec.z_score_img",subset,0),extract["NC"],'source')
	task.connect(extract["NC"],"mat_file",save_csv["NC"],"ext")
	
	task.connect(sca,("outputspec.Z_score",subset,0),gunzip,'in_file')
	task.connect(gunzip,"out_file",extract["SCA"],'source')
	
	task.connect(extract["SCA"],"mat_file",save_csv["SCA"],"ext")
	
	task.connect(save_csv["ALFF"],"csv_file",datasink,"csv.@par1")
	task.connect(save_csv["ReHO"],"csv_file",datasink,"csv.@par2")
	task.connect(save_csv["NC"],"csv_file",datasink,"csv.@par3")
	task.connect(save_csv["SCA"],"csv_file",datasink,"csv.@par4")

	task.connect(alff,"outputspec.alff_Z_img",datasink,"data.alff")
	task.connect(reho,"outputspec.z_score",datasink,"data.reho")
	task.connect(sca,"outputspec.Z_score",datasink,"data.sca")
	task.connect(zscore,"outputspec.z_score_img",datasink,"data.nc")
	task.connect(pp,"output.func",datasink,"data.functional")
	task.connect(pp,"output.struct",datasink,"data.structural")
	task.connect(pp,"output.mask",datasink,"data.mask")

	# print
	p_ru = pe.Node(interface=wrap.Print(), name="print_realign_params")
	p_ru.inputs.out_file = "realignment_parameters.ps"
	
	p_struct = pe.Node(interface=wrap.Print(), name="print_anatomical")
	p_struct.inputs.out_file = "brain_anatomical.ps"
	
	p_func = pe.Node(interface=wrap.Print(), name="print_func")
	p_func.inputs.out_file = "brain_func.ps"
	
	p_mask = pe.Node(interface=wrap.Print(), name="print_mask")
	p_mask.inputs.out_file = "brain_mask.ps"
	
	
	task.connect(pp,"output.movement",p_ru,"in_file")
	task.connect(pp,"output.struct",p_struct,"in_file")
	task.connect(pp,"output.func",p_func,"in_file")
	task.connect(pp,"output.mask",p_mask,"in_file")
	
	task.connect(p_ru,"out_file",datasink,"ps.@par1")
	task.connect(p_func,"out_file",datasink,"ps.@par4")
	task.connect(p_struct,"out_file",datasink,"ps.@par5")
	task.connect(p_mask,"out_file",datasink,"ps.@par6")


	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task

def imcalc_expression(x):
	return "i1/"+str(x)	


def flt(directory,sequence):
	import nipype.interfaces.base as base
	import wrappers as wrap
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	subject = get_subject(directory)
		
	work = pe.Workflow(name=sequence)
	work.base_dir = base_dir
	
	flt = pe.Node(wrap.FLT(),name="FLT")
	flt.inputs.in_file = os.path.abspath(directory+"/bhv_flt/eprime_flt.txt")
	flt.inputs.out_file = os.path.abspath(directory+"/bhv_flt/"+subject+"_flt.csv")
	
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	work.connect(flt,"out_file",datasink,"csv")
	work.write_graph(dotfilename=sequence+"-workflow")
	
	return work

"""
EMBARC 2.0 ASL Sequence
"""
def asl(directory,sequence):
	import wrappers as wrap
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.spm as spm 
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.afni as afni		 # afni
	from nipype.interfaces.utility import Function
	import nipype.interfaces.io as nio           # Data i/o

	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	subject = get_subject(directory)
	
	# get components
	ds = datasource(directory,sequence)
	#pp = preprocess()

	realign = pe.Node(interface=spm.Realign(), name="realign")
	realign.inputs.register_to_mean = True
	realign.inputs.write_which = [2,1]
	
	asl_node = pe.Node(interface=wrap.ASL(), name="asl_script")
	
	# now lets do some standard preprocessing
	func_bet = pe.Node(interface=fsl.BET(), name="func_bet")
	func_bet.inputs.mask = True
	func_bet.inputs.frac = 0.5
	func_bet.inputs.robust = True
	func_bet.inputs.vertical_gradient = 0
	
	#struct_bet = pe.Node(interface=fsl.BET(), name="struct_bet")
	#struct_bet.inputs.mask = True
	#struct_bet.inputs.frac = 0.5
	#struct_bet.inputs.robust = True
	#struct_bet.inputs.vertical_gradient = 0
	#struct_bet.inputs.center = [133,115,88]	
	
	struct_bet = pe.Node(interface=afni.SkullStrip(), name="skull_strip")
	struct_bet.inputs.args = " -ld 30 -push_to_edge -no_avoid_eyes "
	struct_bet.inputs.outputtype = 'NIFTI'
	
	
	coregister = pe.Node(interface=fsl.FLIRT(), name='flirt_mean2struct')
	coregister.inputs.cost = 'corratio'
	coregister.inputs.bins = 256
	coregister.inputs.dof = 12
	coregister.inputs.interp = 'trilinear'
	coregister.inputs.searchr_x = [-180, 180]
	coregister.inputs.searchr_y = [-180, 180]
	coregister.inputs.searchr_z = [-180, 180]
	
	coreg_xfm = pe.Node(interface=fsl.ApplyXfm(),name='apply_func2struct')
	coreg_xfm.inputs.interp = 'trilinear'
	
	flirt = pe.Node(interface=fsl.FLIRT(), name='flirt_struct2template')
	flirt.inputs.cost = 'corratio'
	flirt.inputs.bins = 256
	flirt.inputs.dof = 12
	flirt.inputs.interp = 'trilinear'
	flirt.inputs.searchr_x = [-180, 180]
	flirt.inputs.searchr_y = [-180, 180]
	flirt.inputs.searchr_z = [-180, 180]
	flirt.inputs.reference = OASIS_template
	
	xfm = pe.Node(interface=fsl.ApplyXfm(),name='apply_func2template')
	xfm.inputs.interp = 'trilinear'
	xfm.inputs.reference = OASIS_template
	
		
	# apply transform to coregistered 4D functional referenced template
	sxfm = pe.Node(interface=fsl.ApplyXfm(),name='apply_struct2template')
	sxfm.inputs.interp = 'trilinear'
	sxfm.inputs.reference = OASIS_template


	# do fnirt
	fnirt_s2t = pe.Node(interface=fsl.FNIRT(), name='fnirt_struct2template')	
	fnirt_s2t.inputs.ref_file = OASIS_template
	
	warp_f2t  = pe.Node(interface=fsl.ApplyWarp(), name='warp_func2template')	
	warp_f2t.inputs.ref_file = OASIS_template
	
	warp_s2t  = pe.Node(interface=fsl.ApplyWarp(), name='warp_structtemplate')	
	warp_s2t.inputs.ref_file = OASIS_template
		
	
	susan = pe.Node(interface=fsl.SUSAN(), name="smooth")
	susan.inputs.brightness_threshold = 2000.0
	susan.inputs.fwhm = 8
	susan.inputs.output_type = 'NIFTI'

	imcalc = pe.Node(interface=wrap.ImageCalc(), name="imcalc")
	
	rois = [ROI_dir+"BA10_5mm.img",
			ROI_dir+"BR1.nii",
			ROI_dir+"BR2.nii",
			ROI_dir+"BR3.nii",
			ROI_dir+"BR4.nii",
			ROI_dir+"BR9.nii",
			ROI_dir+"pgACC_0_42_4_8mm.img",
			ROI_dir+"sgACC_0_24_-8_8mm.img",
			ROI_dir+"dACC_0_24_38_8mm.img",
			ROI_dir+"FIRST_L_amyg_small.nii",
			ROI_dir+"FIRST_R_amyg_small.nii",
			ROI_dir+"BA09_left.img",
			ROI_dir+"BA09_right.img",
			ROI_dir+"Left_VLPFC_5mm.img",
			ROI_dir+"Right_VLPFC_5mm.img",
			ROI_dir+"L_ant_insula.img",
			ROI_dir+"R_ant_insula.img",
			ROI_dir+"L_ant_MFG_emoconflict_MDDcompensate_p001.img",
			ROI_dir+"R_ant_MFG_emoconflict_MDDcompensate_p001.img",
			ROI_dir+"VS_left.img",
			ROI_dir+"VS_right.img",
			ROI_dir+"putamen_left.img",
			ROI_dir+"putamen_right.img",
			ROI_dir+"caudate_body_left.img",
			ROI_dir+"caudate_body_right.img",
			ROI_dir+"caudate_head_left.img",
			ROI_dir+"caudate_head_right.img"]
	
	roi_names =['BA10_CBF','BR01_CBF','BR02_CBF','BR03_CBF','BR04_CBF',
				'BR09_CBF','pgCing_CBF','sgCing_CBF','dCing_CBF','LeftAmygdala_CBF',
				'RightAmygdala_CBF','LeftBA09_CBF','RightBA09_CBF','LeftBA47_CBF',
				'RightBA47_CBF','LeftInsula_CBF','RightInsula_CBF','LeftMFG_CBF',
				'RightMFG_CBF','LeftVS_CBF','RightVS_CBF','LeftPutamen_CBF',
				'RightPutamen_CBF','LeftCaudateBody_CBF','RightCaudateBody_CBF',
				'LeftCaudateHead_CBF','RightCaudateHead_CBF']
				
	roi = pe.Node(interface=wrap.ROIExtractor(), name="roi_extract")
	roi.inputs.roi_images = rois
	roi.inputs.average = 'vox'
	roi.inputs.interpelation = 0
	
	save_csv = pe.Node(name="save_csv",
			interface=Function(input_names=["task","units","names","ext","output"],
			output_names=["csv_file"],function=wrap.save_csv))
	save_csv.inputs.task = "Cerebral_Blood_Flow"
	save_csv.inputs.units = "mL/min/100g"
	save_csv.inputs.names = roi_names
	save_csv.inputs.output = subject+"_"+sequence+"_outcomes_asl.csv"
	
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	outputnode = pe.Node(interface=util.IdentityInterface(fields=['csv_file']),name='output')
	
	
	asl = pe.Workflow(name=sequence)
	asl.base_dir = base_dir
	
	#pp.disconnect()
	#asl.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])
	#asl.connect(pp,"realign.realigned_files",asl_node,"in_files")
	#asl.connect(asl_node,'cbf_image',pp,"apply_func2struct.in_file")
	#asl.connect(pp,"susan.smoothed_file",imcalc,"in_file")
	
	asl.connect(ds,'func',realign,'in_files')
	asl.connect(realign,'realigned_files',asl_node,'in_files')
	
	asl.connect(realign,'mean_image',func_bet,'in_file') 
	asl.connect(ds,'struct',struct_bet,'in_file')
	
	asl.connect(struct_bet,'out_file',coregister,'reference')	
	asl.connect(func_bet,'out_file',coregister,'in_file')
		
	asl.connect(coregister,'out_matrix_file',coreg_xfm,'in_matrix_file')
	asl.connect(asl_node,'cbf_image',coreg_xfm,'in_file')
	asl.connect(struct_bet,'out_file',coreg_xfm,'reference')
	asl.connect(struct_bet,'out_file',flirt,'in_file')
	asl.connect(coreg_xfm,'out_file',xfm,'in_file')
	asl.connect(flirt,'out_matrix_file',xfm,'in_matrix_file')

	asl.connect(struct_bet,"out_file",sxfm,"in_file")
	asl.connect(flirt,'out_matrix_file',sxfm,'in_matrix_file')         
	
	asl.connect(sxfm,'out_file',fnirt_s2t,'in_file')
	asl.connect(xfm,'out_file',warp_f2t,'in_file')
	asl.connect(fnirt_s2t,'field_file',warp_f2t,'field_file')
	asl.connect(sxfm,'out_file',warp_s2t,'in_file')
	asl.connect(fnirt_s2t,'field_file',warp_s2t,'field_file')
	
	asl.connect(warp_f2t,'out_file',susan, 'in_file') 
	#asl.connect(xfm,'out_file',susan, 'in_file') 
	asl.connect(susan,'smoothed_file',imcalc,'in_file')
	asl.connect(asl_node,('cbf_value',imcalc_expression),imcalc,'expression')

	asl.connect(imcalc,"out_file",roi,'source')
	asl.connect(roi,"mat_file",save_csv,'ext')
	asl.connect(save_csv,"csv_file",outputnode,"csv_file")
	asl.connect(save_csv,"csv_file",datasink,"csv")
	
	
	# print
	p_ru = pe.Node(interface=wrap.Print(), name="print_realign_params")
	p_ru.inputs.out_file = "realignment_parameters.ps"
	
	p_asl = pe.Node(interface=wrap.Print(), name="print_cbf")
	p_asl.inputs.out_file = "brain_cbf.ps"
			
	asl.connect(realign,"realign_parameters",p_ru,"in_file")
	asl.connect(asl_node,'cbf_image',p_asl,"in_file")
	
	asl.connect(p_ru,"out_file",datasink,"ps.@par1")
	asl.connect(p_asl,"out_file",datasink,"ps.@par2")

	asl.write_graph(dotfilename=sequence+"-workflow")
	return asl;
	

# check sequence
def check_sequence(opt_list,directory,seq):
	seq_dir = "/dicom_bold_"+seq
	if seq == "asl":
		seq_dir = "/dicom_"+seq
	elif seq == "flt":
		seq_dir = "/bhv_flt"

	# check if directory exists
	if not os.path.exists(directory+seq_dir):
		print "Error: data directory for "+seq+" does not exists, skipping .."
		print "Missing directory: "+directory+seq_dir
		return False

	# if no sequence specified, check QC failed condition
	if len(opt_list) == 0:
		files = []
		files.append(directory+seq_dir+"/FAIL.txt")
		files.append(directory+seq_dir+"/FAIL_checked.txt")
		if seq != "flt":
			files.append(directory+"/dicom_anatomical/FAIL.txt")
			files.append(directory+"/dicom_anatomical/FAIL_checked.txt")
		
		for f in files:
			if os.path.exists(f):
				print "Error: looks like "+seq+" has failed QA, skipping .."
				print "QA file: "+f
				return False
		return True
	# else if sequence specified, do it and ignore failed condition	
	elif "-"+seq in opt_list:
		return True
	# else 	
	return False




# run pipeline if used as standalone script
if __name__ == "__main__":	
	opts = "[-asl|-ert|-resting1|-resting2|-reward|-flt]"
	opt_list = []
	
	# get arguments
	if len(sys.argv) < 2:
		print "Usage: embarc.py "+opts+" <embarc subject directory>"
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
	
	if "subject" == get_subject(directory):
		print "Error: "+directory+" is not a valid EMBARC data directory"
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
                      	  'hash_method': 'timestamp'})
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
	
	if check_sequence(opt_list,directory,"ert"):
		log.info("\n\nERT pipeline ...\n\n")
		t = time.time()		
		ert = task(directory,'ert')
		ert.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
		#ert.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
	
	if check_sequence(opt_list,directory,"reward"):
		log.info("\n\nREWARD pipeline ...\n\n")
		t = time.time()		
		reward = task(directory,'reward')
		reward.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"resting1"):
		log.info("\n\nRESTING1 pipeline ...\n\n")
		t = time.time()		
		resting1 = resting(directory,'resting1')
		resting1.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"resting2"):
		log.info("\n\nRESTING2 pipeline ...\n\n")
		t = time.time()		
		resting2 = resting(directory,'resting2')
		resting2.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
		
	if check_sequence(opt_list,directory,"asl"):
		log.info("\n\nASL pipeline ...\n\n")
		t = time.time()		
		asl = asl(directory,'asl')
		asl.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
	
	if check_sequence(opt_list,directory,"flt"):
		log.info("\n\nFLT pipeline ...\n\n")
		t = time.time()		
		flt = flt(directory,"flt")
		flt.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
		
	log.info("\n\npipeline complete\n\n")
