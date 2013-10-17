#!/usr/bin/env python
# EMBARC 2.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import sys
import os                                  
import re


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
	
	d = os.path.dirname(os.path.realpath(__file__))+"/../validation/templates/"
	template = d+"OASIS-30_Atropos_template_in_MNI152_2mm.nii.gz"
	
	preproc = pe.Workflow(name='preprocess')

	inputnode = pe.Node(interface=util.IdentityInterface(
				fields=['func','struct']),name='input')

	realign = pe.Node(interface=spm.Realign(), name="realign")
	realign.inputs.register_to_mean = True

	func_bet = pe.Node(interface=fsl.BET(), name="func_bet")
	func_bet.inputs.mask = True
	func_bet.inputs.frac = 0.5
	func_bet.inputs.robust = True
	func_bet.inputs.vertical_gradient = 0
	
	struct_bet = pe.Node(interface=fsl.BET(), name="struct_bet")
	struct_bet.inputs.mask = True
	struct_bet.inputs.frac = 0.5
	struct_bet.inputs.robust = True
	struct_bet.inputs.vertical_gradient = 0
	struct_bet.inputs.center = [133,115,88]	
	
	#coregister = pe.Node(interface=spm.Coregister(), name="coregister")
	#coregister.inputs.jobtype = 'estimate'
	coregister = pe.Node(interface=fsl.FLIRT(), name='coregister')
	coregister.inputs.cost = 'corratio'
	coregister.inputs.bins = 256
	coregister.inputs.dof = 12
	coregister.inputs.interp = 'trilinear'
	coregister.inputs.searchr_x = [-180, 180]
	coregister.inputs.searchr_y = [-180, 180]
	coregister.inputs.searchr_z = [-180, 180]
	
	coreg_xfm = pe.Node(interface=fsl.ApplyXfm(),name='coreg_xfm')
	coreg_xfm.inputs.interp = 'trilinear'
	#coreg_xfm.inputs.reference = template
	
	flirt = pe.Node(interface=fsl.FLIRT(), name='flirt')
	flirt.inputs.cost = 'corratio'
	flirt.inputs.bins = 256
	flirt.inputs.dof = 12
	flirt.inputs.interp = 'trilinear'
	flirt.inputs.searchr_x = [-180, 180]
	flirt.inputs.searchr_y = [-180, 180]
	flirt.inputs.searchr_z = [-180, 180]
	flirt.inputs.reference = template
	
	xfm = pe.Node(interface=fsl.ApplyXfm(),name='apply_xfm')
	#xfm.inputs.cost = 'corratio'
	#xfm.inputs.bins = 256
	#xfm.inputs.dof = 12
	xfm.inputs.interp = 'trilinear'
	#xfm.inputs.searchr_x = [-180, 180]
	#xfm.inputs.searchr_y = [-180, 180]
	#xfm.inputs.searchr_z = [-180, 180]
	xfm.inputs.reference = template
	
	susan = pe.Node(interface=fsl.SUSAN(), name="smooth")
	susan.inputs.brightness_threshold = 2000.0
	susan.inputs.fwhm = 8

	despike = pe.Node(interface=afni.Despike(), name='despike')
	despike.inputs.outputtype = 'NIFTI'

	outputnode = pe.Node(interface=util.IdentityInterface(
				 fields=['func','mask','movement']),name='output')


	preproc.connect([(inputnode,realign,[('func','in_files')]),
                     (realign,func_bet,[('mean_image','in_file')]), 
					 (inputnode,struct_bet,[('struct','in_file')]),
					 
					 #(struct_bet,coregister,[('out_file','source')]),
					 #(func_bet,coregister,[('out_file','target')]),
					 
					 (struct_bet,coregister,[('out_file','reference')]),
					 (func_bet,coregister,[('out_file','in_file')]),
					 
					 (coregister,coreg_xfm,[('out_matrix_file','in_matrix_file')]),
					 (realign,coreg_xfm,[('realigned_files','in_file')]),
					 (struct_bet,coreg_xfm,[('out_file','reference')]),
					 
					 #(realign,coregister,[('realigned_files','apply_to_files')]),
					 
					 #(coregister,flirt,[('coregistered_source','in_file')]),
					 (struct_bet,flirt,[('out_file','in_file')]),
		             #(struct_bet,flirt,[('out_file','in_file')]),
		             #(realign,xfm,[('realigned_files','in_file')]),
		             (coreg_xfm,xfm,[('out_file','in_file')]),
		             (flirt,xfm,[('out_matrix_file','in_matrix_file')]),
		             
		             (xfm,despike,[('out_file','in_file')]),
		             (despike,susan, [('out_file', 'in_file')]), 
		             #(xfm, susan, [('out_file', 'in_file')]), 
		             (susan,outputnode,[('smoothed_file','func')]),
		             (func_bet,outputnode,[('mask_file','mask')]),
		             (realign,outputnode,[('realignment_parameters','movement')])
		             ])
	return preproc



# convert eprime to design matrix
def eprime2dm(eprime):
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
	m.inputs.script = sequence+"_eprime2dm(\'"+eprime+"\')"
	m.run();

	# extract variables from .mat files
	mat = os.path.join(os.path.dirname(eprime),'nDM*.mat')
	mat = gl.glob(mat)
	if len(mat) > 0:
		mat = mat[0]
	dm = sp.loadmat(mat,squeeze_me=True)
		
	# load up values and convert them
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
		#cont1 = ('RewardExpectancy','T', ['anticipationxanti^1'],[0, 0, 1])
		#cont2 = ('PredictionError','T', ['outcomexsignedPE^1'],[0, 0, 0, 0, 1])
		cont1 = ('RewardExpectancy','T', ['anticipationxanti^1'],[1])
		cont2 = ('PredictionError','T', ['outcomexsignedPE^1'],[1])
		return [cont1,cont2]
	return []


"""
EMBARC 1.0 Level1 Analysis Pipeline
"""
def level1analysis():
	import nipype.interfaces.spm as spm          # spm
	import nipype.interfaces.matlab as mlab      # how to run matlab
	import nipype.interfaces.utility as util     # utility
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.algorithms.modelgen as model   # model specification
	
	mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
	
	l1analysis = pe.Workflow(name='level1analysis')
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

	contrastestimate = pe.Node(interface = spm.EstimateContrast(), name="contrastestimate")
	
	outputnode = pe.Node(interface=util.IdentityInterface(
			fields=['contrasts','con_images','spmT_images']),name='output')
	
	
	l1analysis.connect([(inputnode,modelspec,[('movement','realignment_parameters'),
					  						('func','functional_runs'),
					  						(('behav', eprime2dm),'subject_info')]),	
					  (modelspec,level1design,[('session_info','session_info')]),
					  #(inputnode,level1design,[('mask','mask_image')]),
		              (level1design,level1estimate,[('spm_mat_file','spm_mat_file')]),
		              (inputnode,contrastestimate,[(('behav',get_contrasts),'contrasts')]),
		              (level1estimate,contrastestimate,[('spm_mat_file','spm_mat_file'),
		                                              ('beta_images','beta_images'),
		                                              ('residual_image','residual_image')]),
		              (contrastestimate,outputnode,[('spm_mat_file','contrasts'),
		              		              			('con_images','con_images'),
		              		              			('spmT_images','spmT_images')])		   
		              ])
	return l1analysis		             

"""
EMBARC 1.0 Input Data Source
"""
def datasource(directory, sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o

	# define some variables beforehand
	m = re.search('embarc_CU_([A-Z0-9]+_\dR\d)_mri',directory)
	subject=m.group(1)

	# specify input dataset just pass through parameters
	datasource = pe.Node(interface=nio.DataGrabber(
						 infields=['subject_id','sequence'], 
						 outfields=['func', 'struct','behav']),
	                     name = 'datasource')
	datasource.inputs.base_directory = os.path.abspath(directory)
	datasource.inputs.template = '*'
	datasource.inputs.field_template = dict(func='dicom_bold_%s/%s_bold_%s.nii',
                                    		struct='dicom_anatomical/%s_anatomical.nii',
                                    		behav='dicom_bold_%s/eprime_%s.txt')
	datasource.inputs.template_args  = dict(func=[['sequence','subject_id','sequence']],
										    struct=[['subject_id']],
										    behav=[['sequence','sequence']])
	datasource.inputs.subject_id = subject
	datasource.inputs.sequence = sequence

	return datasource


"""
EMBARC 1.0 Task Sequence Ex: ert/reward
"""
def task(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
		
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	
	# get components
	ds = datasource(directory,sequence)
	pp = preprocess()
	l1 = level1analysis();
	
	
	# connect components into a pipeline
	ert = pe.Workflow(name=sequence)
	ert.base_dir = base_dir
	ert.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])
	ert.connect([(pp,l1,[('output.func','input.func'),
						 ('output.mask','input.mask'),
						 ('output.movement','input.movement')])])
	ert.connect([(ds,l1,[('behav','input.behav')])])
	ert.write_graph(dotfilename=sequence+"-workflow") #,graph2use='flat')
	return ert



# run pipeline if used as standalone script
if __name__ == "__main__":	
	# get arguments
	if len(sys.argv) < 2:
		print "Usage: embarc.py <embarc subject directory>"
		sys.exit(1)
	
	# logging verbosity
	import nipype
	import logging
	from nipype import config
	
	# pick dataset that we'll be wroking on
	directory = sys.argv[1]
	
	
	
	# setup logging, display and other config
	disp = os.environ['DISPLAY']
	log_dir = directory+"/logs"
	if not os.path.exists(log_dir):
		os.mkdir(log_dir)
	cfg = dict(logging={'workflow_level':'INFO',
						'interface_level':'INFO',
						'log_to_file': True,
						'log_directory': log_dir},
    		   execution={'stop_on_first_crash': True,
    		   			  'display_variable':disp,
                      	  'hash_method': 'timestamp'})
	config.update_config(cfg)
	log = nipype.logging.getLogger('workflow')
	l = nipype.logging.getLogger('interface').parent.handlers[0]
	l.setFormatter(logging.Formatter('%(name)-2s %(levelname)-2s:\t %(message)s'))
	###########
	
	
	log.info("\n\nERT pipeline ...\n\n")
	ert = task(directory,'ert')
	ert.run()
	#ert.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
	
	#do reward
	log.info("\n\nREWARD pipeline ...\n\n")
	reward = task(directory,'reward')
	reward.run()
