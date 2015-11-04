#!/usr/bin/env python
# EMBARC 2.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import os   


"""
params is the structure with parameters that is study specific
"""
class Config:
	def __init__(self):
		import sys
		import os                                  
		import re
		
		self.CPU_CORES = 16
		self.time_repetition  = 1.5
		
		self.bet_mask = True
		self.bet_frac = 0.5
		self.bet_robust = True
		self.bet_vertical_gradient = 0
		
		self.flirt_cost = 'corratio'
		self.flirt_bins = 256
		self.flirt_dof = 12
		self.flirt_interp = 'trilinear'
		self.flirt_searchr_x = [-180, 180]
		self.flirt_searchr_y = [-180, 180]
		self.flirt_searchr_z = [-180, 180]
		
		self.susan_brightness_threshold = 200.0
		self.susan_fwhm = 6
		
		self.modelspec_concatenate_runs   = False
		self.modelspec_high_pass_filter_cutoff = 256
		self.modelspec_input_units = 'secs'
		
		self.level1design_bases = {'hrf':{'derivs': [1,0]}}
		self.level1design_timing_units = 'secs'
		self.level1estimate_estimation_method = {'Classical' : 1}
		self.contrastestimate_use_derivs = True
			
		bin_dir = os.path.dirname(os.path.realpath(__file__))
		data_dir = os.environ.get("SUPPORT_DATA_DIR",bin_dir+"/../data")

		self.OASIS_template = data_dir+"/templates/OASIS-30_Atropos_template_in_MNI152_2mm.nii.gz"
		self.OASIS_labels = data_dir+"/templates/OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_2mm.nii.gz"

		ROI_dir = data_dir+"/OASIS_rois/"
		ROI_suffix = "_oasis.nii"
		self.FNIRT_config = os.getenv("FSLDIR")+"/etc/flirtsch/T1_2_MNI152_2mm.cnf"
		self.subject_site=""
		########################
		## predefine all ROIs ##
		########################
		self.ROI_white = ROI_dir+"white_matter_mask"+ROI_suffix
		self.ROI_L_insula = ROI_dir+"left_insula"+ROI_suffix
		self.ROI_R_insula = ROI_dir+"right_insula"+ROI_suffix
		self.ROI_L_amyg = ROI_dir+"left_amygdala"+ROI_suffix
		self.ROI_R_amyg = ROI_dir+"left_amygdala"+ROI_suffix
		self.ROI_VS_L = ROI_dir+"left_accumbens_area"+ROI_suffix
		self.ROI_VS_R = ROI_dir+"right_accumbens_area"+ROI_suffix
		self.ROI_VS_LR = ROI_dir+"bilateral_accumbens_area"+ROI_suffix
		self.ROI_BA9_L = ROI_dir+"left_rostral_middle_frontal"+ROI_suffix
		self.ROI_BA9_R = ROI_dir+"right_rostral_middle_frontal"+ROI_suffix
		self.ROI_BA10 = ROI_dir+"medial_brodmann_area_10"+ROI_suffix
		self.ROI_BR1 = ROI_dir+"beckmann_region_1"+ROI_suffix
		self.ROI_BR2 = ROI_dir+"beckmann_region_2"+ROI_suffix
		self.ROI_BR3 = ROI_dir+"beckmann_region_3"+ROI_suffix
		self.ROI_BR4 = ROI_dir+"beckmann_region_4"+ROI_suffix
		self.ROI_BR9 = ROI_dir+"beckmann_region_9"+ROI_suffix 
		self.ROI_PG_ACC = ROI_dir+"pregenual_ACC"+ROI_suffix
		self.ROI_D_ACC = ROI_dir+"dorsal_ACC"+ROI_suffix
		self.ROI_SG_ACC = ROI_dir+"subgenual_ACC"+ROI_suffix
		self.ROI_L_MFG_compensate = ROI_dir+"left_anterior_MFG_compensate"+ROI_suffix
		self.ROI_R_MFG_compensate = ROI_dir+"right_anterior_MFG_compensate"+ROI_suffix
		self.ROI_L_VLPFC = ROI_dir+"left_brodmann_area_47"+ROI_suffix
		self.ROI_R_VLPFC = ROI_dir+"right_brodmann_area_47"+ROI_suffix
		self.ROI_L_ant_insula =	ROI_dir+"left_anterior_insula"+ROI_suffix
		self.ROI_R_ant_insula = ROI_dir+"right_anterior_insula"+ROI_suffix
		self.ROI_putamen_L = ROI_dir+"left_putamen"+ROI_suffix
		self.ROI_putamen_R = ROI_dir+"right_putamen"+ROI_suffix
		self.ROI_caudate_body_L = ROI_dir+"left_caudate"+ROI_suffix
		self.ROI_caudate_body_R = ROI_dir+"right_caudate"+ROI_suffix
		self.ROI_caudate_head_L = ROI_dir+"left_caudate"+ROI_suffix
		self.ROI_caudate_head_R = ROI_dir+"right_caudate"+ROI_suffix
		



"""
EMBARC/DIAMOND 2.0 PreProcessing Pipeline
input: 
	func  - functional image
	struct - structural image
	template - template image
output:
	func  - functional processed images
	ufunc - unsmoothed functional image
	mask  - mask image from the functional
	movement - realign movement parameters
	struct - structural processed image
"""
def preprocess(config):
	import nipype.interfaces.spm as spm          # spm
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.utility as util     # utility
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.afni as afni	     # afni
	
	fsl.FSLCommand.set_default_output_type('NIFTI')
	
	
	preproc = pe.Workflow(name='preprocess')

	inputnode = pe.Node(interface=util.IdentityInterface(fields=['func','struct']),name='input')

	realign = pe.Node(interface=spm.Realign(), name="realign")
	realign.inputs.register_to_mean = True
	preproc.connect(inputnode,"func",realign,"in_files")


	bet_mean = pe.Node(interface=fsl.BET(), name="bet_mean")
	bet_mean.inputs.mask = config.bet_mask
	bet_mean.inputs.frac = config.bet_frac
	bet_mean.inputs.robust = config.bet_robust
	bet_mean.inputs.vertical_gradient = config.bet_vertical_gradient
	preproc.connect(realign,'mean_image',bet_mean,'in_file') 

	bet_struct = pe.Node(interface=fsl.BET(), name="bet_struct")
	bet_struct.inputs.mask = config.bet_mask
	bet_struct.inputs.frac = config.bet_frac
	bet_struct.inputs.robust = config.bet_robust
	bet_struct.inputs.vertical_gradient = config.bet_vertical_gradient
	preproc.connect(inputnode,'struct',bet_struct,'in_file')	

	# flirt_meam2struct structural to mean functional 
	flirt_m2s = pe.Node(interface=fsl.FLIRT(), name='flirt_mean2struct')
	flirt_m2s.inputs.cost = config.flirt_cost
	flirt_m2s.inputs.bins = config.flirt_bins
	flirt_m2s.inputs.dof = config.flirt_dof
	flirt_m2s.inputs.interp = config.flirt_interp
	flirt_m2s.inputs.searchr_x = config.flirt_searchr_x
	flirt_m2s.inputs.searchr_y = config.flirt_searchr_y
	flirt_m2s.inputs.searchr_z = config.flirt_searchr_z
	preproc.connect(bet_struct,'out_file',flirt_m2s,'reference')
	preproc.connect(bet_mean,'out_file',flirt_m2s,'in_file')
 
	# flirt_struct2template structural to template	            
	flirt_s2t = pe.Node(interface=fsl.FLIRT(), name='flirt_struct2template')
	flirt_s2t.inputs.cost = config.flirt_cost
	flirt_s2t.inputs.bins = config.flirt_bins
	flirt_s2t.inputs.dof = config.flirt_dof
	flirt_s2t.inputs.interp = config.flirt_interp
	flirt_s2t.inputs.searchr_x = config.flirt_searchr_x
	flirt_s2t.inputs.searchr_y = config.flirt_searchr_y
	flirt_s2t.inputs.searchr_z = config.flirt_searchr_z
	flirt_s2t.inputs.reference = config.OASIS_template
	preproc.connect(bet_struct,'out_file',flirt_s2t,'in_file')
	

	# it is expensive to apply func2 struct right away
	# so better to concatanate matricies
	cat_f2s = pe.Node(interface=fsl.ConvertXFM(),name='concat_func2struct')
	cat_f2s.inputs.concat_xfm = True
	preproc.connect(flirt_m2s,'out_matrix_file',cat_f2s,'in_file')
	preproc.connect(flirt_s2t,'out_matrix_file',cat_f2s,'in_file2')

	# apply transform to realigned 4d functional referenced to structural
	apply_f2s = pe.Node(interface=fsl.ApplyXfm(),name='apply_func2struct')
	apply_f2s.inputs.interp =  config.flirt_interp
	apply_f2s.inputs.reference = config.OASIS_template	
	preproc.connect(cat_f2s,'out_file',apply_f2s,'in_matrix_file')
	preproc.connect(realign,'realigned_files',apply_f2s,'in_file')
		
	# apply transform to flirt_mean2structed 4D functional referenced template
	apply_s2t = pe.Node(interface=fsl.ApplyXfm(),name='apply_struct2template')
	apply_s2t.inputs.interp =  config.flirt_interp
	apply_s2t.inputs.reference = config.OASIS_template
	preproc.connect(bet_struct,'out_file',apply_s2t,'in_file')
	preproc.connect(flirt_s2t,'out_matrix_file',apply_s2t,'in_matrix_file')
	

	# apply transform to flirt_mean2structed 4D functional referenced template
	apply_us2t = pe.Node(interface=fsl.ApplyXfm(),name='apply_ustruct2template')
	apply_us2t.inputs.interp =  config.flirt_interp
	apply_us2t.inputs.reference = config.OASIS_template
	preproc.connect(flirt_s2t,'out_matrix_file',apply_us2t,'in_matrix_file')
	preproc.connect(inputnode,'struct',apply_us2t,'in_file')
	
	# do fnirt
	fnirt_s2t = pe.Node(interface=fsl.FNIRT(), name='fnirt_struct2template')	
	fnirt_s2t.inputs.ref_file = config.OASIS_template
	fnirt_s2t.inputs.config_file = config.FNIRT_config
	preproc.connect(apply_us2t,'out_file',fnirt_s2t,'in_file')

	# apply field
	warp_f2t  = pe.Node(interface=fsl.ApplyWarp(), name='warp_func2template')	
	warp_f2t.inputs.ref_file = config.OASIS_template
	preproc.connect(apply_f2s,'out_file',warp_f2t,'in_file')
	preproc.connect(fnirt_s2t,'field_file',warp_f2t,'field_file')
	
	# apply field
	warp_s2t  = pe.Node(interface=fsl.ApplyWarp(), name='warp_struct2template')	
	warp_s2t.inputs.ref_file = config.OASIS_template
	preproc.connect(apply_s2t,'out_file',warp_s2t,'in_file')
	preproc.connect(fnirt_s2t,'field_file',warp_s2t,'field_file')
	
	# remove spikes
	despike = pe.Node(interface=afni.Despike(), name='despike')
	despike.inputs.outputtype = 'NIFTI'
	preproc.connect(warp_f2t,'out_file',despike,'in_file')
	
	# smooth image using SUSAN
	susan = pe.Node(interface=fsl.SUSAN(), name="smooth")
	susan.inputs.brightness_threshold = config.susan_brightness_threshold 
	susan.inputs.fwhm = config.susan_fwhm
	preproc.connect(despike,'out_file',susan,'in_file') 

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
	preproc.connect(warp_s2t,'out_file',outputnode,'struct')
	preproc.connect(bet_func,'mask_file',outputnode,'mask')
	
	return preproc

"""
load Matlab Design Matrix (nDM file)
mat_file -> nDM SPM design matrix file
trim     -> trim design matrix by number of columns

"""
def load_design_matrix(mat_file,trim=0):
	import os
	import re
	import scipy.io as sp
	import glob as gl
	import numpy
	import nipype.interfaces.matlab as mlab 
	from nipype.interfaces.base import Bunch

	# convert numpy data array
	def convert_numpy(ar,toString = False):
		lst = []
		for a in ar.tolist():
			if toString:
				lst.append(str(a))
			elif type(a) == numpy.ndarray:
				if a.size > 1:			
					lst.append(a.tolist())
				else:
					lst.append([a.tolist()])
			else:
				lst.append([a])
		return lst

	# if list of mat_files, then do a list
	mat_files = []	
	if isinstance(mat_file,list):
		mat_files = mat_file
	else:
		mat_files = [mat_file]
		
	bunches = []

	# go over mat files
	for mat_file in mat_files:
		# load design matrix 
		dm = sp.loadmat(mat_file,squeeze_me=True)

		names  = convert_numpy(dm.get('names'),True)
		onsets = convert_numpy(dm.get('onsets'))
		durations = convert_numpy(dm.get('durations'))

		# load up values and convert them
		# for PPPI remove last column for reward PPI; for ert last 3 columns
		# error, posterror, misc
		if trim > 0:
			names = names[0:(len(names)-trim)]
			durations = durations[0:(len(durations)-trim)]
			onsets = onsets[0:(len(onsets)-trim)]

		# create bunch to return
		bunch = Bunch(conditions=names,onsets=onsets,durations=durations)
		if 'pmod' in dm:
			pmod = []
			for i in range(0,len(dm.get('pmod'))):
				if dm['pmod']['name'][i].size == 0:
					pmod.append(None)
				else:
					name = str(dm['pmod']['name'][i])
					param = dm['pmod']['param'][i].tolist()
					poly = dm['pmod']['poly'][i]
					pmod.append(Bunch(name=[name],param=[param],poly=[poly]))
			bunch.pmod = pmod
		bunches.append(bunch)

	return bunches


"""
 Generic method to convert an E-Prime file or other tasks to nDM 
 design matrix
"""
def create_design_matrix(matlab_function, eprime_file):
	import os
   	import re
	import glob as gl
	import nipype.interfaces.matlab as mlab 
	
	# execute matlab script to generate nDM file
	m = mlab.MatlabCommand()
	m.inputs.mfile = False
	m.inputs.script = matlab_function+"(\'"+eprime_file+"\');"
	m.run();

	# get nDM file (again)
	mat = os.path.join(os.path.dirname(eprime_file),'nDM*.mat')
	mat = gl.glob(mat)
	if len(mat) > 0:
		mat = mat[0]
	return mat

"""
EMBARC 1.5 Level1 Analysis Pipeline
"""
def level1analysis(config,trim=0,name='level1'):
	import nipype.interfaces.spm as spm          # spm
	import nipype.interfaces.matlab as mlab      # how to run matlab
	import nipype.interfaces.utility as util     # utility
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.algorithms.modelgen as model   # model specification
	
	#mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
	l1analysis = pe.Workflow(name=name)
	inputnode = pe.Node(interface=util.IdentityInterface(fields=['movement','func','design_matrix','contrasts']),name='input')

	# specify design matrix model
	modelspec = pe.Node(interface=model.SpecifySPMModel(), name= "modelspec")
	modelspec.inputs.concatenate_runs   = config.modelspec_concatenate_runs
	modelspec.inputs.time_repetition = config.time_repetition
	modelspec.inputs.high_pass_filter_cutoff = config.modelspec_high_pass_filter_cutoff
	modelspec.inputs.input_units = config.modelspec_input_units
	l1analysis.connect(inputnode,'movement',modelspec,'realignment_parameters')
	l1analysis.connect(inputnode,'func',modelspec,'functional_runs')
	l1analysis.connect(inputnode,('design_matrix',load_design_matrix,trim),modelspec,'subject_info')

	# create design matrix
	level1design = pe.Node(interface=spm.Level1Design(), name= "level1design")
	level1design.inputs.bases = config.level1design_bases
	level1design.inputs.timing_units = config.level1design_timing_units
	level1design.inputs.interscan_interval = config.time_repetition
	l1analysis.connect(modelspec,'session_info',level1design,'session_info')

	# level 1 estimate
	level1estimate = pe.Node(interface=spm.EstimateModel(), name="level1estimate")
	level1estimate.inputs.estimation_method = config.level1estimate_estimation_method
	l1analysis.connect(level1design,'spm_mat_file',level1estimate,'spm_mat_file')	

	# no need for contrast for pppi model
	contrastestimate = pe.Node(interface = spm.EstimateContrast(), name="contrastestimate")
	contrastestimate.inputs.use_derivs = config.contrastestimate_use_derivs
	
	l1analysis.connect(inputnode,'contrasts',contrastestimate,'contrasts')
	l1analysis.connect(level1estimate,'spm_mat_file',contrastestimate,'spm_mat_file')
	l1analysis.connect(level1estimate,'beta_images',contrastestimate,'beta_images'),
	l1analysis.connect(level1estimate,'residual_image',contrastestimate,'residual_image')

	# output
	outputnode = pe.Node(interface=util.IdentityInterface(fields=['spm_mat_file','con_images','spmT_images']),name='output')
	l1analysis.connect(contrastestimate,'spm_mat_file',outputnode,'spm_mat_file')
	l1analysis.connect(contrastestimate,'con_images',outputnode,'con_images')
	l1analysis.connect(contrastestimate,'spmT_images',outputnode,'spmT_images')
	              		              			

	return l1analysis		             



"""
EMBARC 1.0 Task Sequence Ex: ert/reward
create generic task analysis
pppi_trim_dm - trim the last N columns of design matrix for PPPI analysis
pppi_rois    - list of tuples (roi_name,roi_file) to create a workflow for 

"""
def task(pppi_trim_dm, pppi_rois):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	
		
	# connect components into a pipeline
	task = pe.Workflow(name="task")
	l1 = level1analysis();
	l2 = level1analysis(pppi_trim_dm,"level1_pppi");
	#l2.name = "level1_pppi"

	fields=['subject','func','movement','design_matrix','contrasts','pppi_contrasts']
	inputnode = pe.Node(interface=util.IdentityInterface(fields=fields),name='input')
	task.connect([(inputnode,l1,[('func','input.func'),('design_matrix','input.design_matrix'),
				     ('movement','input.movement'),('contrasts','input.contrasts')])])
	task.connect([(inputnode,l2,[('func','input.func'),('design_matrix','input.design_matrix'),
				     ('movement','input.movement'),('contrasts','input.contrasts')])])

	fields = ['spm_mat_file','con_images']	
	for roi in pppi_rois:
		fields.append("pppi_"+roi[0]+"_con_images")

	outputnode = pe.Node(interface=util.IdentityInterface(fields=fields),name='output')	
	task.connect(l1,"output.con_images",outputnode,"con_images")
	task.connect(l1,"output.spm_mat_file",outputnode,"spm_mat_file")


	# now do gPPI analysis
	for roi in pppi_rois:
		pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
		pppi.inputs.voi_name = roi[0]
		pppi.inputs.voi_file = roi[1]
		task.connect(l2,'output.spm_mat_file',pppi,'spm_mat_file')
		task.connect(inputnode,'subject',pppi,'subject')
		
		#contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
		#task.connect(inputnode,'pppi_contrasts',contrast,'contrasts')
		#task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
		#task.connect(pppi,'beta_images',contrast,'beta_images')
		#task.connect(pppi,'residual_image',contrast,'residual_image')
		#task.connect(contrast,'con_images',outputnode,"pppi_"+roi[0]+"_con_images")
	
	return task

	
"""
add printing and saving of files to a workflow
workflow - where everything will be added
node     - node from which to get output files
files	 - files that need to be printed and saved
"""
def print_save_files(workflow,node,datasink,files):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	# go over parameters
	for param in files:
		prnt= pe.Node(interface=wrap.Print(), name="print_"+param)
		prnt.inputs.out_file = param+".ps"
		workflow.connect(node,param,datasink,"data."+param)
		workflow.connect(node,param,prnt,"in_file")	
		workflow.connect(prnt,"out_file",datasink,"ps.@par"+param)

"""
extract and save a set of ROIs
"""
def extract_save_rois(name,csv_name,task_name,task_units,roi_list,workflow,datasink):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	
	extract = pe.Node(interface=wrap.ROIExtractor(), name="extract_"+name)
	extract.inputs.roi_images = list(zip(*roi_list)[1])
	extract.inputs.average = 'none'
	extract.inputs.interpelation = 0
	#extract.inputs.source = source_image
	
	# save CSV file
	save_csv = pe.Node(name="save_csv_"+name,
		interface=Function(input_names=["task","units","names","ext","output"],
		output_names=["csv_file"],function=wrap.save_csv))
	save_csv.inputs.task = task_name
	save_csv.inputs.units = task_units
	save_csv.inputs.names = list(zip(*roi_list)[0])
	save_csv.inputs.output = csv_name
	

	workflow.connect(extract,"mat_file",save_csv,"ext")
	workflow.connect(save_csv,"csv_file",datasink,"csv.@par"+name)

	return extract

