#!/usr/bin/env python
# DIAMOND 1.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import sys
import os                                  
import re
#import embarc
import gold

## Predefined constants ##
conf = gold.Config()
conf.CPU_CORES = 16
conf.time_repetition  = 1.5
		
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
	
conf.fugue_dwell_time = 0.00075 # prizma scanner is default 
conf.fugue_poly_order = 3

conf.dartel_fwhm = 6 #TODO Check value
conf.dartel_voxel_size =  (2, 2, 2)
conf.susan_brightness_threshold = 750 #200.0
conf.susan_fwhm = 6

conf.filter_image_bptf = ' -bptf 83.33 6.67' #hard coded filters, but dependent on TR

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




# default value to use fieldmap in the pipeline
useFieldmap=False

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
	outfields=['func', 'struct']
	field_template = dict(func=sequence+"/*"+sequence+".[ni]*",struct="anat/*_anat_crop.nii")
	template_args  = dict(func=[[]],struct=[[]])                

	if useFieldmap:
		field_template['fieldmap_mag']   = "field_map_hr/*_mag.nii"
		field_template['fieldmap_phase'] ="field_map_hr/*_phase.nii"
		template_args['fieldmap_mag']  = [[]]
		template_args['fieldmap_phase']  = [[]]
		outfields.append('fieldmap_mag')
		outfields.append('fieldmap_phase')
			


	# add behavior file to task oriented design
	#if sequence.startswith('reward') or sequence.startswith('efnback') or sequence.startswith('dynamic_faces'):
	#	field_template['behav'] = sequence+"/*task.txt"
	#	template_args['behav']  = [[]]

	# specify input dataset just pass through parameters
	datasource = pe.Node(interface=nio.DataGrabber(
						 infields=['subject_id','sequence'], 
						 outfields=outfields),
	                     name = 'datasource_'+sequence)
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
	m = re.search('encore_[0-9]+.([0-9]+)',directory)
	if m:
		return m.group(1)
	return "subject"

def subset(x,i):
	return x[i]	


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
	
	
	
	# setup some constants
	resting_roi_names = ['LeftInsula','RightInsula','LeftAmygdala',
			     'RightAmygdala','LeftVS','RightVS','LeftBA9','RightBA9',
			     'BR1','BR2','BR3','BR4','BR9','LeftBA47','RightBA47','LeftPutamen','RightPutamen',
			     'LeftCaudateHead','RightCaudateHead'] #, 'leftVLPFC'

	resting_roi_images = [conf.ROI_L_insula,conf.ROI_R_insula,conf.ROI_L_amyg,conf.ROI_R_amyg,
			conf.ROI_VS_L,conf.ROI_VS_R,conf.ROI_BA9_L,conf.ROI_BA9_R,
			conf.ROI_BR1,conf.ROI_BR2,conf.ROI_BR3,conf.ROI_BR4,conf.ROI_BR9,
			conf.ROI_L_VLPFC,conf.ROI_R_VLPFC, conf.ROI_putamen_L, conf.ROI_putamen_R,
			conf.ROI_caudate_head_L,conf.ROI_caudate_head_R] #, 				conf.ROI_leftVLPFC


	ds = datasource(directory,sequence)
	pp = gold.preprocess_mni(conf,useFieldmap)
	
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
	for hl in [[0.01,0.1],[0.0,0.04],[0.04,0.08],[0.08,0.12],[0.12,0.16],[0.16,0.20],[0.20,0.24]]: #[[0.01,0.1],[0.01,0.027],[0.027,0.073]]: # add more frequencies
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
	
	task.connect(nu,"regressors",column_select,"in_file")
	task.connect(column_select,"out_file",glm_NGS,"design")
	task.connect(pp,"output.func",glm_NGS,"in_file")

	task.connect(filt,"out_file",roiave,"in_file")
	task.connect(roiave,"out_file",corroi,"in_files")
		
	
	for nm in alff_nm:
		#task.connect(glm,'out_res',alff[nm],'inputspec.rest_res')
		task.connect(glm_NGS,'out_res',alff[nm],'inputspec.rest_res')
		task.connect(pp,'output.mask',alff[nm],'inputspec.rest_mask')	

	task.connect(filt,"out_file",reho,"inputspec.rest_res_filt")
	task.connect(pp,"output.mask",reho,"inputspec.rest_mask")
	
	task.connect(glm,'out_res',nc,'inputspec.subject')
	task.connect(nc,'outputspec.centrality_outputs',zscore,'inputspec.input_file')
	task.connect(pp,'output.mask',zscore,'inputspec.mask_file')

	for mask in resting_roi_names:
		#["BR9","LeftVS","RightVS","BR2","BR3"]: # add left vlpfc
		sca[mask] = CPAC.sca.create_sca(name_sca="sca_"+mask);
		maskave[mask] = pe.Node(interface=afni.Maskave(),name="roi_ave_"+mask)
		maskave[mask].inputs.outputtype = "NIFTI"
		maskave[mask].inputs.quiet= True
		maskave[mask].inputs.mask = resting_roi_images[resting_roi_names.index(mask)]
		gunzip[mask] = pe.Node(interface=misc.Gunzip(),name="gunzip_"+mask)

		task.connect(filt,"out_file",maskave[mask],"in_file")
		task.connect(filt,"out_file",sca[mask],"inputspec.functional_file")
		task.connect(maskave[mask],"out_file",sca[mask],"inputspec.timeseries_one_d")
		task.connect(sca[mask],("outputspec.Z_score",subset,0),gunzip[mask],'in_file')
		task.connect(sca[mask],"outputspec.Z_score",datasink,"data.sca."+mask)
	
	
	task.connect(reho,"outputspec.z_score",datasink,"data.reho")
	task.connect(corroi,"out_file",datasink,"csv.@par5")
	for nm in alff_nm:	
		task.connect(alff[nm],"outputspec.alff_Z_img",datasink,"data."+nm.lower())
	
	task.connect(zscore,"outputspec.z_score_img",datasink,"data.nc")
	

	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task

# run simply preprocessing
def preprocess(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o
	import gold

	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# get components
	ds = datasource(directory,sequence)
	pp = gold.preprocess_mni(conf)
	#pp.get_node("input").inputs.template = embarc.OASIS_template	
	
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])
				
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
		
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task


# check sequence
def check_sequence(opt_list,directory,seq):
	seq_dir = seq
	# check if directory exists
	if not (os.path.exists(directory+"/"+seq_dir) or os.path.exists(directory+"/"+seq_dir+"_1")):
		print "Error: data directory for "+seq+" does not exists, skipping .."
		print "Missing directory: "+directory+seq_dir
		return False

	# if no sequence specified, check QC failed condition
	#if len(opt_list) == 0:
		#files = []
		#files.append(directory+seq_dir+"/FAIL.txt")
		#files.append(directory+seq_dir+"/FAIL_checked.txt")
		#if seq != "flt":
		#	files.append(directory+"/dicom_anatomical/FAIL.txt")
		#	files.append(directory+"/dicom_anatomical/FAIL_checked.txt")
		
		#for f in files:
		#	if os.path.exists(f):
		#		print "Error: looks like "+seq+" has failed QA, skipping .."
		#		print "QA file: "+f
		#		return False
	#	return True
	# else if sequence specified, do it and ignore failed condition	
	#elif "-"+seq in opt_list:
	#	return True
	# else 	
	return True


# run pipeline if used as standalone script
if __name__ == "__main__":	
	opts = "[-resting_state_hr|-fieldmap|-trio]"
	opt_list = []
	
	# get arguments
	if len(sys.argv) < 2:
		print "Usage: encore.py "+opts+" <encore subject directory>"
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
		print "Error: "+directory+" is not a valid ENCORE data directory"
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
	mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodisplay -nosplash")
	mlab.MatlabCommand.set_default_terminal_output('stream')
	#mlab.MatlabCommand.set_default_paths(bin_dir)
	
	if "-fieldmap" in opt_list:	
		useFieldmap = True

	# change dwell time based on scanner
	if "-trio" in opt_list:
		opt_list.remove("-trio")
		conf.fugue_dwell_time = 0.00088   # For encore Trio (currently in gold.py)
	else:
		conf.fugue_dwell_time = 0.00075   # For encore Prisma + Impress Prisma



	if check_sequence(opt_list,directory,"resting_state_hr"):
		log.info("\n\nRESTING pipeline ...\n\n")
		t = time.time()		
		resting1 = resting(directory,"resting_state_hr")
		resting1.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
		
	log.info("\n\npipeline complete\n\n")
