#!/usr/bin/env python
# DIAMOND 1.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import sys
import os                                  
import re
import embarc


"""
EMBARC 1.0 Input Data Source
"""
def datasource(directory, sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o

	# define some variables beforehand
	subject=get_subject(directory)
	
	# define templates for datasource
	field_template = dict(func=sequence+"/*.img",struct="anat/T1MPRAGE*[0-9].nii")
	template_args  = dict(func=[[]],struct=[[]])                


	# add behavior file to task oriented design
	if sequence.startswith('reward'):
		field_template['behav'] = sequence+"/*Reward_Task*.txt"
		template_args['behav']  = [[]]
	elif sequence.startswith('efnback'):
		field_template['behav'] = sequence+"/EFNBACK_NewEye*.txt"
		template_args['behav']  = [[]]
	elif sequence.startswith('dynamic_faces'):
		field_template['behav'] = sequence+"/subject*.txt"
		template_args['behav']  = [[]]

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

def subset(x,i):
	return x[i]	

"""
EMBARC 2.0 Reward sequence
"""
def reward(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	import gold	

	subject = get_subject(directory)
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# some hard-coded sequence specific components
	# Reward contrasts
	cont1 = ('RewardExpectancy','T', ['anticipationxanti^1'],[1])
	cont2 = ('PredictionError','T', ['outcomexsignedPE^1'],[1])
	contrasts = [cont1,cont2]

	cont1 = ('RewardExpectancy','T', ['PPI_anticipationxanti^1'],[1])
	cont2 = ('PredictionError','T', ['PPI_outcomexsignedPE^1'],[1])
	pppi_contrasts = [cont1,cont2]

	# get components
	pp = gold.preprocess(embarc.OASIS_template)
	ts = gold.task(1,[("Reward_VS",embarc.ROI_VS_LR)])
	ts.inputs.input.subject = subject
	ts.inputs.input.contrasts = contrasts
	ts.inputs.input.pppi_contrasts = pppi_contrasts

	#dm = pe.
		
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])
	task.connect([(pp,ts,[('output.func','input.func'),('output.movement','input.movement')])])
	task.connect(dm,"design_matrix",ts,"design_matrix")
	
	# define datasink
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
	gold.print_save_files(task,ts.get_node('output'),datasink,("spm_mat_file","con_images","pppi_Reward_VS_con_images"))	
	
	# extract specific ROIs	
	anticipation_rois =    [("BA10_anticipation",embarc.ROI_BA10),
				("LeftBA9_anticipation",embarc.ROI_BA9_L),
				("RightBA9_anticipation",embarc.ROI_BA9_R),
				("LeftVS_anticipation",embarc.ROI_VS_L),
				("RightVS_anticipation",embarc.ROI_VS_R),
				("LeftBA47_anticipation",embarc.ROI_L_VLPFC),
				("RightBA47_anticipation",embarc.ROI_R_VLPFC),
				("BR1_anticipation",embarc.ROI_BR1),
				("BR2_anticipation",embarc.ROI_BR2),
				("BR3_anticipation",embarc.ROI_BR3),
				("BR4_anticipation",embarc.ROI_BR4)]
	outcome_rois =   [("LeftBA9_outcome",embarc.ROI_BA9_L),
			  ("RightBA9_outcome",embarc.ROI_BA9_R),
			  ("LeftVS_outcome",embarc.ROI_VS_L),
			  ("RightVS_outcome",embarc.ROI_VS_R),
			  ("BR1_outcome",embarc.ROI_BR1),
			  ("BR2_outcome",embarc.ROI_BR2),
			  ("BR3_outcome",embarc.ROI_BR3),
			  ("BR4_outcome",embarc.ROI_BR4)]
	pppi_rois  = 	[("BR1_anticipation_PPI",embarc.ROI_BR1),
			 ("BR2_anticipation_PPI",embarc.ROI_BR2),
			 ("BR3_anticipation_PPI",embarc.ROI_BR3),
			 ("BR4_anticipation_PPI",embarc.ROI_BR4)]


	csv_name = subject+"_"+sequence+"_outcomes_anticipation.csv"
	ext = gold.extract_save_rois("anticipation",csv_name,"RewardTask","ParameterEstimate",anticipation_rois,task,datasink)
	task.connect(ts,("output.con_images",subset,0),ext,"source")

	csv_name = subject+"_"+sequence+"_outcomes_outcome.csv"
	ext = gold.extract_save_rois("outcome",csv_name,"RewardTask","ParameterEstimate",outcome_rois,task,datasink)
	task.connect(ts,("output.con_images",subset,1),ext,"source")

	csv_name = subject+"_"+sequence+"_outcomes_anticipation.csv"
	ext = gold.extract_save_rois("anticipation_pppi",csv_name,"RewardTask","ParameterEstimate",pppi_rois,task,datasink)
	task.connect(ts,("output.pppi_Reward_VS_con_images",subset,0),ext,"source")
	

	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task

"""
EMBARC 2.0 EFNBACK sequence
"""
def efnback(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	import gold	

	subject = get_subject(directory)
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# some hard-coded sequence specific components
	# Reward contrasts
	cont1 = ('RewardExpectancy','T', ['anticipationxanti^1'],[1])
	cont2 = ('PredictionError','T', ['outcomexsignedPE^1'],[1])
	contrasts = [cont1,cont2]

	cont1 = ('RewardExpectancy','T', ['PPI_anticipationxanti^1'],[1])
	cont2 = ('PredictionError','T', ['PPI_outcomexsignedPE^1'],[1])
	pppi_contrasts = [cont1,cont2]

	# get components
	ds1 = datasource(directory,sequence+"_1")
	ds2 = datasource(directory,sequence+"_2")		
	pp1 = gold.preprocess(embarc.OASIS_template)
	pp1.name = "preprocess_1"	
	pp2 = gold.preprocess(embarc.OASIS_template)
	pp2.name = "preprocess_2"

	# merge inputs
	m1 = pe.Node(interface=util.Merge(), name="merge")
	m2 = pe.Node(interface=fsl.Merge(), name="fsl_merge")

	# actual task
	ts = gold.task(1,[("Reward_VS",embarc.ROI_VS_LR)])
	ts.inputs.input.subject = subject
	ts.inputs.input.contrasts = contrasts
	ts.inputs.input.pppi_contrasts = pppi_contrasts

	# create DesignMatrix
	dm = pe.Node(name="create_DM",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm.inputs.matlab_function = "efnback_task2dm"

	
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds1,pp1,[('func','input.func'),('struct','input.struct')])])
	task.connect([(ds2,pp2,[('func','input.func'),('struct','input.struct')])])
	
	task.connect([(pp,ts,[('output.func','input.func'),('output.movement','input.movement')])])
	task.connect(ds,'behav',dm,"eprime_file")	
	task.connect(dm,"design_matrix",ts,"input.design_matrix")
	
	# define datasink
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
	gold.print_save_files(task,ts.get_node('output'),datasink,("spm_mat_file","con_images","pppi_Reward_VS_con_images"))	
	
	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task

"""
GOLD 2.0 dynamic faces sequence
"""
def dynamic_faces(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	import gold
	import embarc	

	subject = get_subject(directory)
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# some hard-coded sequence specific components
	contrasts = []	
	contrasts.append(("Anger","T",["Anger*bf(1)"],[1]))
	contrasts.append(("Anger_td","T",["Anger*bf(2)"],[1]))
	contrasts.append(("Fear","T",["Fear*bf(1)"],[1]))
	contrasts.append(("Fear_td","T",["Fear*bf(2)"],[1]))
	contrasts.append(("Sad","T",["Sad*bf(1)"],[1]))
	contrasts.append(("Sad_td","T",["Sad*bf(2)"],[1]))
	contrasts.append(("Happy","T",["Happy*bf(1)"],[1]))
	contrasts.append(("Happy_td","T",["Happy*bf(2)"],[1]))
	contrasts.append(("IDmorph","T",["IDmorph*bf(1)"],[1]))
	contrasts.append(("IDmorph_td","T",["IDmorph*bf(2)"],[1]))
	contrasts.append(("Shape","T",["Shape*bf(1)"],[1]))
	contrasts.append(("Shape_td","T",["Shape*bf(2)"],[1]))
	contrasts.append(("Anger > Shape","T",["Anger*bf(1)","Shape*bf(1)"],[1,-1]))
	contrasts.append(("Fear > Shape","T",["Fear*bf(1)","Shape*bf(1)"],[1,-1]))
	contrasts.append(("Sad > Shape","T",["Sad*bf(1)","Shape*bf(1)"],[1,-1]))
	contrasts.append(("Happy > Shape","T",["Happy*bf(1)","Shape*bf(1)"],[1,-1]))
	contrasts.append(("Emotion > Shape","T",["Anger*bf(1)","Fear*bf(1)","Sad*bf(1)","Happy*bf(1)","Shape*bf(1)"],[.25,.25,.25,.25,-1]))

	# get components
	ds = datasource(directory,sequence)	
	pp = gold.preprocess(embarc.OASIS_template)
	l1 = gold.level1analysis();
	l1.inputs.input.contrasts = contrasts
	
	# mCompCor
	cc = pe.Node(interface=wrap.mCompCor(), name="mCompCor")
	cc.inputs.white_mask = embarc.ROI_white

	# create DesignMatrix
	dm = pe.Node(name="create_DM",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm.inputs.matlab_function = "dynfaces_task2dm"

	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])
	task.connect([(pp,cc,[('output.ufunc','source'),('output.mask','brain_mask'),('output.movement','movement')])])
	task.connect(cc,"regressors",l1,"input.movement")	
	task.connect(pp,'input.func',l1,'input.func')
	task.connect(ds,'behav',dm,"eprime_file")	
	task.connect(dm,'design_matrix',l1,"input.design_matrix")
	
	# define datasink
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
	gold.print_save_files(task,l1.get_node('output'),datasink,("spm_mat_file","con_images"))	
	
	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task



# run resting sequence
def resting(directory):
	import embarc
	ds = datasource(directory,"resting_state")
	return embarc.resting(directory,"resting_state",get_subject(directory),ds)


# run simply preprocessing
def preprocess(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o
	import embarc
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
	pp = gold.preprocess(embarc.OASIS_template)
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
	opts = "[-dynamic_faces|-efnback|-reward|-resting_state]"
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
		log.info("\n\nREWARD pipeline ...\n\n")
		t = time.time()		
		reward = reward(directory,"reward")
		reward.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
	if check_sequence(opt_list,directory,"resting_state"):
		log.info("\n\nRESTING pipeline ...\n\n")
		t = time.time()		
		resting1 = resting(directory)
		resting1.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"efnback"):
		log.info("\n\nEFNBACK pipeline ...\n\n")
		t = time.time()		
		efnback = efnback(directory,"efnback")
		efnback.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
	if check_sequence(opt_list,directory,"dynamic_faces"):
		log.info("\n\nDynamic_Faces 2 pipeline ...\n\n")
		t = time.time()		
		df = dynamic_faces(directory,"dynamic_faces")
		df.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
		
	log.info("\n\npipeline complete\n\n")
