#!/usr/bin/env python
# DIAMOND 1.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import sys
import os                                  
import re
import embarc
import gold

## Predefined constants ##
conf = gold.Config()



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
	field_template = dict(func=sequence+"/*"+sequence+".img",struct="anat/*_anat.nii")
	template_args  = dict(func=[[]],struct=[[]])                


	# add behavior file to task oriented design
	#if sequence.startswith('reward') or sequence.startswith('efnback') or sequence.startswith('dynamic_faces'):
	#	field_template['behav'] = sequence+"/*task.txt"
	#	template_args['behav']  = [[]]

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
	m = re.search('encore_[0-9]+.([0-9]+)',directory)
	if m:
		return m.group(1)
	return "subject"

def subset(x,i):
	return x[i]	


# run resting sequence
def resting(directory):
	import embarc
	ds = datasource(directory,"resting_state_hr")
	return embarc.resting(directory,"resting_state_hr",get_subject(directory),ds)


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
	pp = gold.preprocess(conf)
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
	opts = "[-resting_state_hr]"
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
	
	
	
	if check_sequence(opt_list,directory,"resting_state_hr"):
		log.info("\n\nRESTING pipeline ...\n\n")
		t = time.time()		
		resting1 = resting(directory)
		resting1.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
		
	log.info("\n\npipeline complete\n\n")
