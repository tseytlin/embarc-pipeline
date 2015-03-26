#!/usr/bin/env python
# EMBARC 2.0 Analysis Pipeline (resting)
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
	field_template = dict(func='resting_state/resting_state+lams_resting_state.nii.gz',
                          struct='mprage/MPRAGE+lams_sagittal_mprage.nii.gz')
	
	template_args  = dict(func=[[]],struct=[[]])                

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
	m = re.search('Ss([0-9]+)',directory)
	if m:
		return m.group(1)
	return "subject"


"""
EMBARC 1.0 Resting Sequence Ex: resting1/resting2
"""
def resting(directory):
	import embarc
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.algorithms.misc as misc
	sequence = "resting"
	subject = get_subject(directory)
	ds = datasource(directory,sequence)
	task = embarc.resting(directory,sequence,subject,ds)
	
	gunzip = pe.Node(interface=misc.Gunzip(),name="gunzip_func")
	
	
	task.disconnect([(ds,task.get_node('preprocess'),[('func','input.func'),('struct','input.struct')])])
	
	task.connect(task.get_node('datasource'),'func',gunzip,'in_file')
	task.connect(gunzip,'out_file',task.get_node('preprocess'),'input.func')
	task.connect(task.get_node('datasource'),'struct',task.get_node('preprocess'),'input.struct')
	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')

	return task

# run pipeline if used as standalone script
if __name__ == "__main__":	
	
	# get arguments
	if len(sys.argv) < 1:
		print "Usage: resting "+opts+" <subject directory>"
		sys.exit(1)
	
	# logging verbosity
	import time
	import nipype
	import logging
	from nipype import config
	import nipype.interfaces.matlab as mlab 
	directory = sys.argv[1]
	
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
	
	log.info("\n\nRESTING pipeline ...\n\n")
	t = time.time()		
	resting1 = resting(directory)
	resting1.run()
	log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
			
	log.info("\n\npipeline complete\n\n")
