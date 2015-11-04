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
	
	# figure out if this is new or old naming convention file
	orig = False	
	m = gl.glob(os.path.join(directory,"anat/T1MPRAGE*[0-9].nii"))
	if len(m) > 0:
		orig = True

	if orig:
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
	else:

		# define templates for datasource
		field_template = dict(func=sequence+"/*"+sequence+".img",struct="anat/*_anat.nii")
		template_args  = dict(func=[[]],struct=[[]])                


		# add behavior file to task oriented design
		if sequence.startswith('reward') or sequence.startswith('efnback') or sequence.startswith('dynamic_faces'):
			field_template['behav'] = sequence+"/*task.txt"
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
	contrasts = []
	contrasts.append(('RewardExpectancy','T', ['anticipationxanti^1'],[1]))
	contrasts.append(('PredictionError','T', ['outcomexsignedPE^1'],[1]))

	ppi_contrasts = []
	ppi_contrasts.append(('RewardExpectancy','T', ['PPI_anticipationxanti^1'],[1]))
	ppi_contrasts.append(('PredictionError','T', ['PPI_outcomexsignedPE^1'],[1]))

	# get components
	ds1 = datasource(directory,sequence+"_1")
	ds2 = datasource(directory,sequence+"_2")
	
	pp1 = gold.preprocess(conf)
	pp2 = gold.preprocess(conf)
		
	l1 = gold.level1analysis(conf);
	l1.inputs.input.contrasts = contrasts

	l2 = level1analysis(1,"level1_pppi");
	l2.inputs.input.contrasts = contrasts

	# create DesignMatrix
	dm1 = pe.Node(name="create_DM1",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm1.inputs.matlab_function = "reward_eprime2dm"

	dm2 = pe.Node(name="create_DM2",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm2.inputs.matlab_function = "reward_eprime2dm"

	# setup merge points
	merge_func = pe.Node(name="merge_func",interface=util.Merge())
	merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge())
	merge_move = pe.Node(name="merge_movement",interface=util.Merge())
		
		
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds1,pp1,[('func','input.func'),('struct','input.struct')])])
	task.connect([(ds2,pp2,[('func','input.func'),('struct','input.struct')])])
	task.connect(ds1,'behav',dm1,"eprime_file")	
	task.connect(ds2,'behav',dm2,"eprime_file")	
	task.connect(pp1,'output.func',merge_func,'in1')
	task.connect(pp2,'output.func',merge_func,'in2')
	task.connect(pp1,'output.movement',merge_move,'in1')
	task.connect(pp2,'output.movement',merge_move,'in2')
	task.connect(dm1,'design_matrix',merge_nDM,'in1')
	task.connect(dm2,'design_matrix',merge_nDM,'in2')

	task.connect(merge_move,"out",l1,"input.movement")	
	task.connect(mere_func,'out',l1,'input.func')
	task.connect(merge_nDM,'out',l1,"input.design_matrix")
	
	task.connect(merge_move,"out",l2,"input.movement")	
	task.connect(mere_func,'out',l2,'input.func')
	task.connect(merge_nDM,'out',l2,"input.design_matrix")

	# define datasink
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir

	
	# now define PPI
	pppi_rois  = [("Reward_VS",conf.ROI_VS_LR)]	

	# now do gPPI analysis
	for roi in pppi_rois:
		pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
		pppi.inputs.voi_name = roi[0]
		pppi.inputs.voi_file = roi[1]
		pppi.inputs.subject = subject
		task.connect(l,'output.spm_mat_file',pppi,'spm_mat_file')
		
		contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
		contrast.inputs.contrasts = ppi_contrasts
		task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
		task.connect(pppi,'beta_images',contrast,'beta_images')
		task.connect(pppi,'residual_image',contrast,'residual_image')
		task.connect(contrast,'con_images',datasink,"data.pppi_"+roi[0]+"_con_images")
		task.connect(pppi,"spm_mat_file",datasink,"data.ppi_"+roi[0]+"_spm_file")
	

	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp1.get_node('output'),datasink,("func1","movement1","struct1","mask1"))	
	gold.print_save_files(task,pp2.get_node('output'),datasink,("func2","movement2","struct2","mask2"))	
	gold.print_save_files(task,l1.get_node('output'),datasink,("spm_mat_file","con_images","pppi_Reward_VS_con_images"))	


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
	contrasts = []	
	
	# get components
	ds1 = datasource(directory,sequence+"_1")
	ds2 = datasource(directory,sequence+"_2")
	
	pp1 = gold.preprocess(conf)
	pp2 = gold.preprocess(conf)
		
	l1 = gold.level1analysis(conf);
	l1.inputs.input.contrasts = contrasts
	
	# create DesignMatrix
	dm1 = pe.Node(name="create_DM1",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm1.inputs.matlab_function = "efnback_eprime2dm"

	dm2 = pe.Node(name="create_DM2",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm2.inputs.matlab_function = "efnback_eprime2dm"

	# setup merge points
	merge_func = pe.Node(name="merge_func",interface=util.Merge())
	merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge())
	merge_move = pe.Node(name="merge_movement",interface=util.Merge())
		
		
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds1,pp1,[('func','input.func'),('struct','input.struct')])])
	task.connect([(ds2,pp2,[('func','input.func'),('struct','input.struct')])])
	task.connect(ds1,'behav',dm1,"eprime_file")	
	task.connect(ds2,'behav',dm2,"eprime_file")	
	task.connect(pp1,'output.func',merge_func,'in1')
	task.connect(pp2,'output.func',merge_func,'in2')
	task.connect(pp1,'output.movement',merge_move,'in1')
	task.connect(pp2,'output.movement',merge_move,'in2')
	task.connect(dm1,'design_matrix',merge_nDM,'in1')
	task.connect(dm2,'design_matrix',merge_nDM,'in2')

	task.connect(merge_move,"out",l1,"input.movement")	
	task.connect(mere_func,'out',l1,'input.func')
	task.connect(merge_nDM,'out',l1,"input.design_matrix")
	
	# define datasink
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir

	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp1.get_node('output'),datasink,("func1","movement1","struct1","mask1"))	
	gold.print_save_files(task,pp2.get_node('output'),datasink,("func2","movement2","struct2","mask2"))	
	gold.print_save_files(task,l1.get_node('output'),datasink,("spm_mat_file","con_images"))	

	# now define PPI
	ppi_contrasts = []	
	#ppi_contrasts.append(("Anger","T",["PPI_Anger"],[1]))
	#ppi_contrasts.append(("Anger_td","T",["PPI_Anger"],[1]))
	#pppi_rois  = 	[("left_amygdala",embarc.ROI_L_amyg),
	#		 ("right_amygdala",embarc.ROI_R_amyg),
	#		 ("left_VLPFC",embarc.ROI_L_VLPFC),
	#		 ("right_VLPFC",embarc.ROI_R_VLPFC)]	

	# now do gPPI analysis
	#for roi in pppi_rois:
	#	pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
	#	pppi.inputs.voi_name = roi[0]
	#	pppi.inputs.voi_file = roi[1]
	#	pppi.inputs.subject = subject
	#	task.connect(l1,'output.spm_mat_file',pppi,'spm_mat_file')
	#	
	#	contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
	#	contrast.inputs.contrasts = ppi_contrasts
	#	task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
	#	task.connect(pppi,'beta_images',contrast,'beta_images')
	#	task.connect(pppi,'residual_image',contrast,'residual_image')
	#	task.connect(contrast,'con_images',datasink,"data.pppi_"+roi[0]+"_con_images")
	#	task.connect(pppi,"spm_mat_file",datasink,"data.ppi_"+roi[0]+"_spm_file")
	
	
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
	pp = gold.preprocess(conf)
	l1 = gold.level1analysis(conf);
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
	task.connect(pp,'output.func',l1,'input.func')
	task.connect(ds,'behav',dm,"eprime_file")	
	task.connect(dm,'design_matrix',l1,"input.design_matrix")
	
	# define datasink
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir

	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
	gold.print_save_files(task,l1.get_node('output'),datasink,("spm_mat_file","con_images"))	

	# now define PPI
	ppi_contrasts = []	
	ppi_contrasts.append(("Anger","T",["PPI_Anger"],[1]))
	ppi_contrasts.append(("Fear","T",["PPI_Fear"],[1]))
	ppi_contrasts.append(("Sad","T",["PPI_Sad"],[1]))
	ppi_contrasts.append(("Happy","T",["PPI_Happy"],[1]))

	#ppi_contrasts.append(("Anger_td","T",["PPI_Anger"],[1]))
	
	pppi_rois  = 	[("left_amygdala",embarc.ROI_L_amyg),
			 ("right_amygdala",embarc.ROI_R_amyg),
			 ("left_VLPFC",embarc.ROI_L_VLPFC),
			 ("right_VLPFC",embarc.ROI_R_VLPFC)]	

	# now do gPPI analysis
	for roi in pppi_rois:
		pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
		pppi.inputs.voi_name = roi[0]
		pppi.inputs.voi_file = roi[1]
		pppi.inputs.subject = subject
		task.connect(l1,'output.spm_mat_file',pppi,'spm_mat_file')
		
		contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
		contrast.inputs.contrasts = ppi_contrasts
		task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
		task.connect(pppi,'beta_images',contrast,'beta_images')
		task.connect(pppi,'residual_image',contrast,'residual_image')
		task.connect(contrast,'con_images',datasink,"data.pppi_"+roi[0]+"_con_images")
		task.connect(pppi,"spm_mat_file",datasink,"data.ppi_"+roi[0]+"_spm_file")
	
	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task




# run resting sequence
def resting(directory):
	import embarc
	ds = datasource(directory,"resting_state")
	return embarc.resting(directory,"resting_state",get_subject(directory),ds)




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
	nu.inputs.white_mask = ROI_white

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
	#nc.inputs.centrality_options.method_options=[True, True]
	#nc.inputs.centrality_options.weight_options=[True, True]
	nc.inputs.inputspec.method_option=0
	nc.inputs.inputspec.weight_options=[True, True]	
	nc.inputs.inputspec.threshold_option = 1
	nc.inputs.inputspec.threshold = 0.0744 
	nc.inputs.inputspec.template = OASIS_labels
	zscore =  CPAC.network_centrality.get_zscore(wf_name='z_score')
	
	sca = dict()
	maskave = dict()
	gunzip = dict()
	
	for mask in ["BR9","LeftVS","RightVS","BR2","BR3"]:
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
	corroi.inputs.task_name = task_name
	corroi.inputs.out_file = subject+"_"+sequence+"_outcomes_CORR.csv"

	extract = dict()
	save_csv = dict()
	
	#for nm in [alff_nm[0],alff_nm[1],alff_nm[2],"ReHO","NC","SCA_BR9","SCA_LeftVS","SCA_RightVS"]: 
	#	ext = pe.Node(interface=wrap.ROIExtractor(),name="extract_"+nm)
	#	ext.inputs.roi_images = resting_roi_images
	#	ext.inputs.average = 'none'
	#	ext.inputs.interpelation = 0
	#	extract[nm] = ext
	#	
	#	csv = pe.Node(name="save_"+nm,interface=Function(
	#		input_names=["task","units","names","ext","output"],
	#		output_names=["csv_file"],function=wrap.save_csv))
	#	csv.inputs.task = task_name
	#	csv.inputs.units = task_units
	#	csv.inputs.names = [nm+"_"+ s for s in resting_roi_names]
	#	csv.inputs.output = subject+"_"+sequence+"_outcomes_"+nm+".csv"
	#	save_csv[nm] = csv
	
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

	for mask in ["BR9","LeftVS","RightVS"]:
		task.connect(filt,"out_file",maskave[mask],"in_file")
		#task.connect(pp,"output.mask",maskave,"mask")
		task.connect(filt,"out_file",sca[mask],"inputspec.functional_file")
		task.connect(maskave[mask],"out_file",sca[mask],"inputspec.timeseries_one_d")
		task.connect(sca[mask],("outputspec.Z_score",subset,0),gunzip[mask],'in_file')
		#task.connect(gunzip[mask],"out_file",extract["SCA_"+mask],'source')
		task.connect(sca[mask],"outputspec.Z_score",datasink,"data.sca."+mask)
		#task.connect(extract["SCA_"+mask],"mat_file",save_csv["SCA_"+mask],"ext")
		#task.connect(save_csv["SCA_"+mask],"csv_file",datasink,"csv.@par"+mask)

	#for nm in alff_nm:
	#	#task.connect(alff[nm],"outputspec.alff_Z_img",extract[nm],'source')
	#	task.connect(alff[nm],"outputspec.falff_Z_img",extract[nm],'source')
	#	task.connect(extract[nm],"mat_file",save_csv[nm],"ext")

	#task.connect(reho,"outputspec.z_score",extract["ReHO"],'source')
	#task.connect(extract["ReHO"],"mat_file",save_csv["ReHO"],"ext")
	
	#task.connect(zscore,("outputspec.z_score_img",subset,0),extract["NC"],'source')
	#task.connect(extract["NC"],"mat_file",save_csv["NC"],"ext")
	
	#for nm in alff_nm:
	#	task.connect(save_csv[nm],"csv_file",datasink,"csv.@par"+nm)
	#task.connect(save_csv["ReHO"],"csv_file",datasink,"csv.@par2")
	#task.connect(save_csv["NC"],"csv_file",datasink,"csv.@par3")
	#task.connect(save_csv["SCA_RightVS"],"csv_file",datasink,"csv.@par4")
	
	task.connect(reho,"outputspec.z_score",datasink,"data.reho")
	task.connect(corroi,"out_file",datasink,"csv.@par5")
	for nm in alff_nm:	
		task.connect(alff[nm],"outputspec.alff_Z_img",datasink,"data."+nm.lower())
	
	task.connect(zscore,"outputspec.z_score_img",datasink,"data.nc")
	

	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
	#gold.print_save_files(task,l1.get_node('output'),datasink,("spm_mat_file","con_images"))	



	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task



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
	if not (os.path.exists(directory+seq_dir) or os.path.exists(directory+seq_dir+"_1")):
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
	
	

	
	if check_sequence(opt_list,directory,"reward"):
		log.info("\n\nREWARD pipeline ...\n\n")
		t = time.time()		
		reward = reward(directory,"reward")
		reward.run()
		#reward = preprocess(directory,"reward_2")
		#reward.run()
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
		efnback = preprocess(directory,"efnback_1")
		efnback.run()
		efnback = preprocess(directory,"efnback_2")
		efnback.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
	if check_sequence(opt_list,directory,"dynamic_faces"):
		log.info("\n\nDynamic_Faces 2 pipeline ...\n\n")
		t = time.time()		
		df = dynamic_faces(directory,"dynamic_faces")
		df.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
		
	log.info("\n\npipeline complete\n\n")
