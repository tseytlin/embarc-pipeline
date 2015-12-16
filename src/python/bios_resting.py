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

"""
EMBARC 1.0 Input Data Source
"""
def datasource(directory, sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o

	# define some variables beforehand
	subject=get_subject(directory)
	
	# define templates for datasource
	field_template = dict(func=sequence+'/*.nii',
                          struct='*mprage/*_crop.nii')
	
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
			     'BR1','BR2','BR3','BR4','BR9']
	resting_roi_images = [conf.ROI_L_insula,conf.ROI_R_insula,conf.ROI_L_amyg,conf.ROI_R_amyg,
			conf.ROI_VS_L,conf.ROI_VS_R,conf.ROI_BA9_L,conf.ROI_BA9_R,
			conf.ROI_BR1,conf.ROI_BR2,conf.ROI_BR3,conf.ROI_BR4,conf.ROI_BR9]
	
	ds = datasource(directory,sequence)
	pp = gold.preprocess(conf)
	
	nu = pe.Node(interface=wrap.Nuisance(), name="nuisance")
	nu.inputs.white_mask = conf.ROI_white

	glm = pe.Node(interface=fsl.GLM(), name="glm")
	glm.inputs.out_res_name = "residual.4d.nii"
	
	filt = pe.Node(interface=fsl.ImageMaths(), name="filter")
	#filt.inputs.op_string = ' -bptf 128 12.5 '
	filt.inputs.op_string = ' -bptf 37 4.167'
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
	zscore =  CPAC.network_centrality.get_cent_zscore(wf_name='z_score')
	
	sca = dict()
	maskave = dict()
	#gunzip = dict()
	
	for mask in ["BR9","LeftVS","RightVS","BR2","BR3"]:
		sca[mask] = CPAC.sca.create_sca(name_sca="sca_"+mask);
		maskave[mask] = pe.Node(interface=afni.Maskave(),name="roi_ave_"+mask)
		maskave[mask].inputs.outputtype = "NIFTI"
		maskave[mask].inputs.quiet= True
		maskave[mask].inputs.mask = resting_roi_images[resting_roi_names.index(mask)]
		#gunzip[mask] = pe.Node(interface=misc.Gunzip(),name="gunzip_"+mask)
	
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

	for mask in ["BR9","LeftVS","RightVS","BR2","BR3"]:
		task.connect(filt,"out_file",maskave[mask],"in_file")
		task.connect(filt,"out_file",sca[mask],"inputspec.functional_file")
		task.connect(maskave[mask],"out_file",sca[mask],"inputspec.timeseries_one_d")
		#task.connect(sca[mask],("outputspec.Z_score",subset,0),gunzip[mask],'in_file')
		task.connect(sca[mask],"outputspec.Z_score",datasink,"data.sca."+mask)
	
	
	task.connect(reho,"outputspec.raw_reho_map",datasink,"data.reho")
	task.connect(corroi,"out_file",datasink,"csv.@par5")
	for nm in alff_nm:	
		task.connect(alff[nm],"outputspec.alff_img",datasink,"data."+nm.lower())
	
	task.connect(zscore,"outputspec.z_score_img",datasink,"data.nc")
	

	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
	
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
	resting1 = resting(directory,'lams_resting_state')
	resting1.run()
	log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))
			
	log.info("\n\npipeline complete\n\n")
