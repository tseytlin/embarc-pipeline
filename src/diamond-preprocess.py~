#!/usr/bin/env python
# DIAMOND 1.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import sys
import os                                  
import re



"""
EMBARC 1.0 Input Data Source
"""
def datasource(directory, sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o

	# define some variables beforehand
	subject=get_subject(directory)
	
	# define templates for datasource
	field_template = dict(func=sequence+"/*.img",struct="anat/T1MPRAGE*[0-9].nii",behav=sequence+"/*Task*.txt")
	template_args  = dict(func=[[]],struct=[[]],behav=[[]])                

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

"""
EMBARC get subject name from a directory
"""
def get_subject(directory):
	m = re.search('diamond_[0-9]+.([0-9]+)',directory)
	if m:
		return m.group(1)
	return "subject"

"""
EMBARC 2.0 Reward sequence
"""
def reward(directory,sequence):
	import embarc
	params = dict()
	params["Task Name"] =  "RewardTask"
	params["Task Units"] = "ParameterEstimate"
	params["Level1 Names"] = ["anticipation","outcome"]
	# 1st level ROIs
	params["Level1 ROIs"] =[[("BA10_anticipation",embarc.ROI_BA10),
					("LeftBA9_anticipation",embarc.ROI_BA9_L),
					("RightBA9_anticipation",embarc.ROI_BA9_R),
					("LeftVS_anticipation",embarc.ROI_VS_L),
					("RightVS_anticipation",embarc.ROI_VS_R),
					("LeftBA47_anticipation",embarc.ROI_L_VLPFC),
					("RightBA47_anticipation",embarc.ROI_R_VLPFC),
					("BR1_anticipation",embarc.ROI_BR1),
					("BR2_anticipation",embarc.ROI_BR2),
					("BR3_anticipation",embarc.ROI_BR3),
					("BR4_anticipation",embarc.ROI_BR4)],
					
				   [("LeftBA9_outcome",embarc.ROI_BA9_L),
					("RightBA9_outcome",embarc.ROI_BA9_R),
					("LeftVS_outcome",embarc.ROI_VS_L),
					("RightVS_outcome",embarc.ROI_VS_R),
					("BR1_outcome",embarc.ROI_BR1),
					("BR2_outcome",embarc.ROI_BR2),
					("BR3_outcome",embarc.ROI_BR3),
					("BR4_outcome",embarc.ROI_BR4)]]
	# ROIs for PPI
	params["PPI ROIs"] = [("Reward_VS",embarc.ROI_VS_LR,
					[("BR1_anticipation_PPI",embarc.ROI_BR1),
					 ("BR2_anticipation_PPI",embarc.ROI_BR2),
					 ("BR3_anticipation_PPI",embarc.ROI_BR3),
					 ("BR4_anticipation_PPI",embarc.ROI_BR4)])]
	
	# Reward Level 1 models
	params["Level1 Model"] = embarc.design_matrix();
	params["PPPI Model"]   = embarc.design_matrix();
	

	# Reward contrasts
	cont1 = ('RewardExpectancy','T', ['anticipationxanti^1'],[1])
	cont2 = ('PredictionError','T', ['outcomexsignedPE^1'],[1])
	params["Contrasts"] = [cont1,cont2]

		
	ds = datasource(directory,sequence)
	return embarc.task(directory,sequence,get_subject(directory),ds,params)

# run resting sequence
def resting(directory):
	import embarc
	ds = datasource(directory,"resting_state")
	return embarc.resting(directory,"resting_state",get_subject(directory),ds)


# run simply preprocessing
def preprocess(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	import embarc	

	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# get components
	ds = datasource(directory,sequence)
	pp = embarc.preprocess()
		
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])
	
			
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
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








# check sequence
def check_sequence(opt_list,directory,seq):
	seq_dir = "/"+seq	
	# check if directory exists
	if not os.path.exists(directory+seq_dir):
		print "Error: data directory for "+seq+" does not exists, skipping .."
		print "Missing directory: "+directory+seq_dir
		return False

	# if no sequence specified, check QC failed condition
	if len(opt_list) == 0:
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
		return True
	# else if sequence specified, do it and ignore failed condition	
	elif "-"+seq in opt_list:
		return True
	# else 	
	return False


# run pipeline if used as standalone script
if __name__ == "__main__":	
	opts = "[-dynamic_faces|-efnback_1|-efnback_2|-reward_1|-reward_2|-resting_state]"
	opt_list = []
	
	# get arguments
	if len(sys.argv) < 2:
		print "Usage: diamond.py "+opts+" <diamond subject directory>"
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
	
	
	if check_sequence(opt_list,directory,"reward_1"):
		log.info("\n\nREWARD 1 pipeline ...\n\n")
		t = time.time()		
		reward = preprocess(directory,"reward_1")
		reward.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"reward_2"):
		log.info("\n\nREWARD 2 pipeline ...\n\n")
		t = time.time()		
		reward = preprocess(directory,"reward_2")
		reward.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
	if check_sequence(opt_list,directory,"efnback_1"):
		log.info("\n\nEfnback 1 pipeline ...\n\n")
		t = time.time()		
		efnback = preprocess(directory,"efnback_1")
		efnback.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"efnback_2"):
		log.info("\n\nEfnback 2 pipeline ...\n\n")
		t = time.time()		
		efnback = preprocess(directory,"efnback_2")
		efnback.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
	
	if check_sequence(opt_list,directory,"dynamic_faces"):
		log.info("\n\nDynamic_Faces 2 pipeline ...\n\n")
		t = time.time()		
		df = preprocess(directory,"dynamic_faces")
		df.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"resting_state"):
		log.info("\n\nRESTING1 pipeline ...\n\n")
		t = time.time()		
		resting1 = resting(directory)
		resting1.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
		
	log.info("\n\npipeline complete\n\n")
