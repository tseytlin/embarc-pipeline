#!/usr/bin/env python
# DIAMOND 1.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import sys
import os                                  
import re
import gold

## Predefined constants ##
conf = gold.Config()

# default value to use fieldmap in the pipeline
useFieldmap=False
noPrint = False


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


	outfields=['func', 'struct']

	if orig:
		# define templates for datasource
		field_template = dict(func=sequence+"/*.img",struct="anat/T1MPRAGE*[0-9].nii")
		template_args  = dict(func=[[]],struct=[[]])                


		# add behavior file to task oriented design
		if sequence.startswith('reward'):
			field_template['behav'] = sequence+"/*reward_*_task*.txt"
			template_args['behav']  = [[]]
		elif sequence.startswith('efnback'):
			field_template['behav'] = sequence+"/EFNBACK_NewEye*.txt"
			template_args['behav']  = [[]]
		elif sequence.startswith('dynamic_faces'):
			field_template['behav'] = sequence+"/subject*.txt"
			template_args['behav']  = [[]]
	else:

		# define templates for datasource
		field_template = dict(func=sequence+"/*"+sequence+".img",struct="anat/*_anat_crop.nii")
		template_args  = dict(func=[[]],struct=[[]])                

		

		# add behavior file to task oriented design
		if sequence.startswith('reward') or sequence.startswith('efnback') or sequence.startswith('dynamic_faces'):
			field_template['behav'] = sequence+"/*task.txt"
			template_args['behav']  = [[]]
			outfields.append('behav')

		if sequence.startswith('asl'):
			field_template['func_ref'] = sequence+"/*"+sequence+"_ref.img"
			template_args['func_ref']  = [[]]
			outfields.append('func_ref')


		if useFieldmap:
			field_template['fieldmap_mag']   = "field_map/*_mag.nii"
			field_template['fieldmap_phase'] ="field_map/*_phase.nii"
			template_args['fieldmap_mag']  = [[]]
			template_args['fieldmap_phase']  = [[]]
			outfields.append('fieldmap_mag')
			outfields.append('fieldmap_phase')
			
	
	# specify input dataset just pass through parameters
	datasource = pe.Node(interface=nio.DataGrabber(
						 infields=['subject_id','sequence'], 
						 outfields=outfields),
	                     name = "datasource_"+sequence)
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

	conf.modelspec_high_pass_filter_cutoff = 60	
	conf.level1design_bases = {'hrf':{'derivs': [0,0]}}

	# some hard-coded sequence specific components
	contrasts = []
	contrasts.append(('RewardExpectancy','T', ['anticipationxreward_expectancy^1*bf(1)'],[1]))
	contrasts.append(('UncertainExpectancy','T', ['anticipationxuncertainty^1*bf(1)'],[1]))
	contrasts.append(('PredictionError','T', ['outcomexprediction_error^1*bf(1)'],[1]))

	ppi_contrasts = []
	ppi_contrasts.append(('RewardExpectancy','T', ['PPI_anticipationxreward_expectancy^1'],[1]))
	ppi_contrasts.append(('PredictionError','T', ['PPI_outcomexprediction_error^1'],[1]))

	# get components
	ds1 = datasource(directory,sequence+"_1")
	ds2 = datasource(directory,sequence+"_2")
	
	pp1 = gold.preprocess2(conf,useFieldmap,"preprocess_1")
	pp2 = gold.preprocess2(conf,useFieldmap,"preprocess_2")
		
	l1 = gold.level1analysis(conf);
	l1.inputs.input.contrasts = contrasts

	l2 = gold.level1analysis(conf,0,"level1_pppi");
	l2.inputs.input.contrasts = contrasts

	# create DesignMatrix
	dm1 = pe.Node(name="create_DM1",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm1.inputs.matlab_function = "reward_eprime2dm"

	dm2 = pe.Node(name="create_DM2",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm2.inputs.matlab_function = "reward_eprime2dm"

	# setup merge points
	merge_func = pe.Node(name="merge_func",interface=util.Merge(2))
	merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge(2))
	merge_move = pe.Node(name="merge_movement",interface=util.Merge(2))
	
	# mCompCor
	cc1 = pe.Node(interface=wrap.mCompCor(), name="mCompCor1")
	cc1.inputs.white_mask = conf.ROI_white
	# mCompCor
	cc2 = pe.Node(interface=wrap.mCompCor(), name="mCompCor2")
	cc2.inputs.white_mask = conf.ROI_white
	
		
	
	
		
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	if useFieldmap:	
		task.connect([(ds1,pp1,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
		task.connect([(ds2,pp2,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
	else:
		task.connect([(ds1,pp1,[('func','input.func'),('struct','input.struct')])])
		task.connect([(ds2,pp2,[('func','input.func'),('struct','input.struct')])])
	task.connect(ds1,'behav',dm1,"eprime_file")	
	task.connect(ds2,'behav',dm2,"eprime_file")	

	task.connect([(pp1,cc1,[('output.ufunc','source'),('output.mask','brain_mask'),('output.movement','movement')])])
	task.connect([(pp1,cc2,[('output.ufunc','source'),('output.mask','brain_mask'),('output.movement','movement')])])
	task.connect(cc1,"regressors",merge_move,"in1")	
	task.connect(cc2,"regressors",merge_move,"in2")


	task.connect(pp1,'output.func',merge_func,'in1')
	task.connect(pp2,'output.func',merge_func,'in2')
	#task.connect(pp1,'output.movement',merge_move,'in1')
	#task.connect(pp2,'output.movement',merge_move,'in2')
	task.connect(dm1,'design_matrix',merge_nDM,'in1')
	task.connect(dm2,'design_matrix',merge_nDM,'in2')


	task.connect(merge_move,"out",l1,"input.movement")	
	task.connect(merge_func,'out',l1,'input.func')
	task.connect(merge_nDM,'out',l1,"input.design_matrix")
	
	task.connect(merge_move,"out",l2,"input.movement")	
	task.connect(merge_func,'out',l2,'input.func')
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
		task.connect(l2,'output.spm_mat_file',pppi,'spm_mat_file')
		
		contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
		contrast.inputs.contrasts = ppi_contrasts
		task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
		task.connect(pppi,'beta_images',contrast,'beta_images')
		task.connect(pppi,'residual_image',contrast,'residual_image')
		task.connect(contrast,'con_images',datasink,"data.pppi_"+roi[0]+"_con_images")
		task.connect(pppi,"spm_mat_file",datasink,"data.ppi_"+roi[0]+"_spm_file")
	

	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp1.get_node('output'),datasink,("struct","mask"))	
	gold.print_save_files(task,l1.get_node('input'),datasink,("func","movement"))	
	gold.print_save_files(task,l1.get_node('output'),datasink,("spm_mat_file","con_images"))


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

	conf.modelspec_high_pass_filter_cutoff = 128
	conf.level1design_bases = {'hrf':{'derivs': [0,0]}}

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
	contrasts.append(("0back fear-0back noface","T",["zerofear*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("0back fear-0back neutral","T",["zerofear*bf(1)","zeroneutral*bf(1)"],[1,-1]))
	contrasts.append(("0back happy-0back noface","T",["zerohappy*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("0back happy-0back neutral","T",["zerohappy*bf(1)","zeroneutral*bf(1)"],[1,-1]))
	contrasts.append(("0back neutral-0back noface","T",["zeroneutral*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("2back fear-2back noface","T",["twofear*bf(1)","twoblank*bf(1)"],[1,-1]))
	contrasts.append(("2back fear-2back neutral","T",["twofear*bf(1)","twoneutral*bf(1)"],[1,-1]))
	contrasts.append(("2back happy-2back noface","T",["twohappy*bf(1)","twoblank*bf(1)"],[1,-1]))
	contrasts.append(("2back happy-2back neutral","T",["twohappy*bf(1)","twoneutral*bf(1)"],[1,-1]))
	contrasts.append(("2back neutral-2back noface","T",["twoneutral*bf(1)","twoblank*bf(1)"],[1,-1]))
	contrasts.append(("2back noface-0back noface","T",["twoblank*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("2back neutral-0back neutral","T",["twoneutral*bf(1)","zeroneutral*bf(1)"],[1,-1]))
	contrasts.append(("2back fear-0back fear","T",["twofear*bf(1)","zerofear*bf(1)"],[1,-1]))
	contrasts.append(("2back happy-0back happy","T",["twohappy*bf(1)","zerohappy*bf(1)"],[1,-1]))
	contrasts.append(("2back emotion-2back noface","T",["twofear*bf(1)","twohappy*bf(1)","twoblank*bf(1)"],[.5,.5,-1]))
	contrasts.append(("0back emotion-0back noface","T",["zerofear*bf(1)","zerohappy*bf(1)","zeroblank*bf(1)"],[.5,.5,-1]))
	contrasts.append(("2back emotion-2back neutral","T",["twofear*bf(1)","twohappy*bf(1)","twoneutral*bf(1)"],[.5,.5,-1]))
	contrasts.append(("0back emotion-0back neutral","T",["zerofear*bf(1)","zerohappy*bf(1)","zeroneutral*bf(1)"],[.5,.5,-1]))	
	

		


	# get components
	ds1 = datasource(directory,sequence+"_1")
	ds2 = datasource(directory,sequence+"_2")
	
	pp1 = gold.preprocess2(conf,useFieldmap,"preprocess_1")
	pp2 = gold.preprocess2(conf,useFieldmap,"preprocess_2")
	
	
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
	merge_func = pe.Node(name="merge_func",interface=util.Merge(2))
	merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge(2))
	merge_move = pe.Node(name="merge_movement",interface=util.Merge(2))
	
	# mCompCor
	cc1 = pe.Node(interface=wrap.mCompCor(), name="mCompCor1")
	cc1.inputs.white_mask = conf.ROI_white
	# mCompCor
	cc2 = pe.Node(interface=wrap.mCompCor(), name="mCompCor2")
	cc2.inputs.white_mask = conf.ROI_white
	
		
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	if useFieldmap:	
		task.connect([(ds1,pp1,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
		task.connect([(ds2,pp2,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
	else:
		task.connect([(ds1,pp1,[('func','input.func'),('struct','input.struct')])])
		task.connect([(ds2,pp2,[('func','input.func'),('struct','input.struct')])])

	task.connect([(pp1,cc1,[('output.ufunc','source'),('output.mask','brain_mask'),('output.movement','movement')])])
	task.connect([(pp1,cc2,[('output.ufunc','source'),('output.mask','brain_mask'),('output.movement','movement')])])
	task.connect(cc1,"regressors",merge_move,"in1")	
	task.connect(cc2,"regressors",merge_move,"in2")	

	task.connect(ds1,'behav',dm1,"eprime_file")	
	task.connect(ds2,'behav',dm2,"eprime_file")	
	task.connect(pp1,'output.func',merge_func,'in1')
	task.connect(pp2,'output.func',merge_func,'in2')
	#task.connect(pp1,'output.movement',merge_move,'in1')
	#task.connect(pp2,'output.movement',merge_move,'in2')
	task.connect(dm1,'design_matrix',merge_nDM,'in1')
	task.connect(dm2,'design_matrix',merge_nDM,'in2')

	task.connect(merge_move,"out",l1,"input.movement")	
	task.connect(merge_func,'out',l1,'input.func')
	task.connect(merge_nDM,'out',l1,"input.design_matrix")
	
	# define datasink
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir

	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp1.get_node('output'),datasink,("struct","mask"))	
	gold.print_save_files(task,l1.get_node('input'),datasink,("func","movement"))	
	gold.print_save_files(task,l1.get_node('output'),datasink,("spm_mat_file","con_images"))	

	
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
	
	conf.modelspec_high_pass_filter_cutoff = 256
	conf.level1design_bases = {'hrf':{'derivs': [1,0]}} #TODO check

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
	pp = gold.preprocess2(conf,useFieldmap)
	l1 = gold.level1analysis(conf);
	l1.inputs.input.contrasts = contrasts
	
	# mCompCor
	cc = pe.Node(interface=wrap.mCompCor(), name="mCompCor")
	cc.inputs.white_mask = conf.ROI_white

	# create DesignMatrix
	dm = pe.Node(name="create_DM",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm.inputs.matlab_function = "dynfaces_task2dm"

	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	if useFieldmap:	
		task.connect([(ds,pp,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
	else:
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
	#ppi_contrasts.append(("Shape","T",["PPI_Shape"],[1]))
 	ppi_contrasts.append(("Emotion > Shape","T",["PPI_Anger","PPI_Fear","PPI_Sad","PPI_Happy","PPI_Shape"],[.25,.25,.25,.25,-1]))
	ppi_contrasts.append(("EmotionNeg > Shape","T",["PPI_Anger","PPI_Fear","PPI_Sad","PPI_Shape"],[.33,.33,.33,-1]))
	ppi_contrasts.append(("EmotionNeg > Happy","T",["PPI_Anger","PPI_Fear","PPI_Sad","PPI_Happy"],[.33,.33,.33,-1]))
	#ppi_contrasts.append(("Emotion > Shape","T",["Anger(1)","Fear(1)","Sad(1)","Happy(1)","Shape(1)"],[.25,.25,.25,.25,-1]))
	#ppi_contrasts.append(("Anger_td","T",["PPI_Anger"],[1]))
	
	pppi_rois  = 	[("left_amygdala",conf.ROI_L_amyg),
			 ("right_amygdala",conf.ROI_R_amyg),
			 ("left_VLPFC",conf.ROI_L_VLPFC),
			 ("right_VLPFC",conf.ROI_R_VLPFC),	
 			 ("beckmann_region_1",conf.ROI_BR1),
			 ("beckmann_region_2",conf.ROI_BR2),
			 ("beckmann_region_3",conf.ROI_BR3),
			 ("beckmann_region_4",conf.ROI_BR4),
			 ("bilateral_amygdala",conf.ROI_amygdala_LR)	]

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
	
	conf.modelspec_high_pass_filter_cutoff = 256
	
	# setup some constants
	resting_roi_names = ['LeftInsula','RightInsula','LeftAmygdala',
			     'RightAmygdala','LeftVS','RightVS','LeftBA9','RightBA9',
			     'BR1','BR2','BR3','BR4','BR9', 'leftVLPFC']
	resting_roi_images = [conf.ROI_L_insula,conf.ROI_R_insula,conf.ROI_L_amyg,conf.ROI_R_amyg,
			conf.ROI_VS_L,conf.ROI_VS_R,conf.ROI_BA9_L,conf.ROI_BA9_R,
			conf.ROI_BR1,conf.ROI_BR2,conf.ROI_BR3,conf.ROI_BR4,conf.ROI_BR9, conf.ROI_leftVLPFC]
	
	ds = datasource(directory,sequence)
	pp = gold.preprocess2(conf,useFieldmap)
	
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
	#filt.inputs.op_string = ' -bptf 128 12.5 '
	filt.inputs.op_string = ' -bptf 37 4.167' #TODO hard coded filters, but dependent on TR
	filt.inputs.terminal_output = 'none'
	
	alff = dict()
	alff_nm = []
	for hl in [[0.01,0.1],[0.0,0.04],[0.04,0.08],[0.08,0.12],[0.12,0.16],[0.16,0.20],[0.20,0.24]]: 
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
	
	for mask in ["BR9","LeftVS","RightVS","BR2","BR3", "leftVLPFC"]:
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

	for mask in ["BR9","LeftVS","RightVS","BR2","BR3", "leftVLPFC"]:
		task.connect(filt,"out_file",maskave[mask],"in_file")
		task.connect(filt,"out_file",sca[mask],"inputspec.functional_file")
		task.connect(maskave[mask],"out_file",sca[mask],"inputspec.timeseries_one_d")
		task.connect(sca[mask],("outputspec.Z_score",subset,0),gunzip[mask],'in_file')
		task.connect(sca[mask],"outputspec.Z_score",datasink,"data.sca."+mask)
	
	
	task.connect(reho,"outputspec.z_score",datasink,"data.reho")
	task.connect(corroi,"out_file",datasink,"csv.@par5")
	# alff_Z_img in 0.3.5 now in 0.3.6 falff_img	
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
	pp = gold.preprocess2(conf)
	#pp = gold.preprocess2_no_fieldmap(conf)	
	#pp.get_node("input").inputs.template = embarc.OASIS_template	
	
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds,pp,[('func','input.func'),('struct','input.struct'),
	('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
				
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
		
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task

"""
ASL sequence
"""
def asl(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o
	import gold
	import nipype.interfaces.spm as spm          # spm
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.fsl.maths as math  
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.afni as afni	     # afni
	import nipype.workflows.fmri.spm.preprocess as dartel # preprocess
	import wrappers as wrap		
	from nipype.interfaces.nipy.preprocess import Trim

	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# get components
	inputnode = datasource(directory,sequence)
	
	# connect components into a pipeline
	preproc = pe.Workflow(name=name)
	preproc.base_dir = base_dir
	
	# lets do a workflow	
	fsl.FSLCommand.set_default_output_type('NIFTI')
	
	
	# realign 4D functional
	realign = pe.Node(interface=spm.Realign(), name="realign")
	realign.inputs.register_to_mean = True
	preproc.connect(inputnode,"func",realign,"in_files")

	# skull strip mean functional image
	bet_mean = pe.Node(interface=fsl.BET(), name="bet_mean")
	bet_mean.inputs.mask = config.bet_mask
	bet_mean.inputs.frac = config.bet_frac
	bet_mean.inputs.robust = config.bet_robust
	bet_mean.inputs.vertical_gradient = config.bet_vertical_gradient
	preproc.connect(realign,'mean_image',bet_mean,'in_file') 

	# skull strip mean structural image
	bet_struct = pe.Node(interface=fsl.BET(), name="bet_struct")
	bet_struct.inputs.mask = config.bet_mask
	bet_struct.inputs.frac = config.bet_frac
	bet_struct.inputs.robust = config.bet_robust
	bet_struct.inputs.vertical_gradient = config.bet_vertical_gradient
	preproc.connect(inputnode,'struct',bet_struct,'in_file')	


	# coregister images
	coreg_func2struct = pe.Node(interface=spm.Coregister(),name="coreg_func2struct")
	coreg_func2struct.inputs.jobtype = "estimate"
 	coreg_func2struct.inputs.cost_function = config.coregister_cost_function
	coreg_func2struct.inputs.separation = config.coregister_separation
	coreg_func2struct.inputs.tolerance = config.coregister_tolerance
	coreg_func2struct.inputs.fwhm = config.coregister_fwhm
	preproc.connect(bet_struct,'out_file',coreg_func2struct,'target')
	preproc.connect(realign,'realigned_files',coreg_func2struct,'apply_to_files')
	preproc.connect(bet_mean,'out_file',coreg_func2struct,'source')

	# coregister reference images
	coreg_ref2struct = pe.Node(interface=spm.Coregister(),name="coreg_ref2struct")
	coreg_ref2struct.inputs.jobtype = "estimate"
 	coreg_ref2struct.inputs.cost_function = config.coregister_cost_function
	coreg_ref2struct.inputs.separation = config.coregister_separation
	coreg_ref2struct.inputs.tolerance = config.coregister_tolerance
	coreg_ref2struct.inputs.fwhm = config.coregister_fwhm
	preproc.connect(bet_struct,'out_file',coreg_ref2struct,'target')
	preproc.connect(inputnode,'func_ref',coreg_ref2struct,'apply_to_files')
	preproc.connect(inputnode,'func_ref',coreg_ref2struct,'source')


	# invoke ASL script to create CBF mean image
	pcasl = pe.Node(interface=wrap.pCASL(), name="pCASL")
	preproc.connect(coreg_func2struct,'coregistered_files',pcasl,'in_file')
	preproc.connect(coreg_ref2struct,'coregistered_files',pcasl,'ref_file')

	# convert image using fslmats
	convert_image = pe.Node(interface=math.MathsCommand(),name='convert_image')
	convert_image.inputs.args = "-mul 1"
	preproc.connect(pcasl,'cbf_image',convert_image,'in_file')
		
	
	# create dartel template
	dartel_template = dartel.create_DARTEL_template()
	dartel_template.inputs.inputspec.template_prefix = 'Template'
	preproc.connect(inputnode, 'struct',dartel_template,'inputspec.structural_files')


	# now lets do normalization with DARTEL
	norm_func =  pe.Node(interface=spm.DARTELNorm2MNI(modulate=True),name='norm_func')	
	norm_func.inputs.fwhm = config.dartel_fwhm
	norm_func.inputs.voxel_size = config.dartel_voxel_size
	preproc.connect(dartel_template,'outputspec.template_file',norm_func,'template_file')
	preproc.connect(dartel_template, 'outputspec.flow_fields', norm_func, 'flowfield_files')
	preproc.connect(convert_image,'out_file',norm_func,'apply_to_files')
	
	# now lets do normalization with DARTEL
	norm_struct =  pe.Node(interface=spm.DARTELNorm2MNI(modulate=True),name='norm_struct')
	norm_struct.inputs.fwhm = config.dartel_fwhm  #TODO Check value
	preproc.connect(dartel_template,'outputspec.template_file',norm_struct,'template_file')
	preproc.connect(dartel_template, 'outputspec.flow_fields', norm_struct, 'flowfield_files')
	preproc.connect(bet_struct,'out_file',norm_struct,'apply_to_files')

		
	# calculated brighness threshold for susan (mean image intensity * 0.75)
	image_mean = pe.Node(interface=fsl.ImageStats(),name='image_mean')	
	image_mean.inputs.op_string = "-M"
	preproc.connect(norm_func,'normalized_files',image_mean,'in_file')


	# smooth image using SUSAN
	susan = pe.Node(interface=fsl.SUSAN(), name="smooth")
	susan.inputs.fwhm = config.susan_fwhm
	preproc.connect(norm_func,'normalized_files',susan,'in_file') 
	preproc.connect(image_mean,('out_stat',create_brightness_threshold),susan,'brightness_threshold') 

		            
	# gather output
	outputnode = pe.Node(interface=util.IdentityInterface(fields=['func','ufunc','mask','movement','struct']),name='output')
	preproc.connect(despike,'out_file',outputnode, 'ufunc')
	preproc.connect(susan,'smoothed_file',outputnode,'func')
	preproc.connect(realign,'realignment_parameters',outputnode,'movement')
	preproc.connect(norm_struct,'normalized_files',outputnode,'struct')
	preproc.connect(bet_func,'mask_file',outputnode,'mask')

	
	# datasing			
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	# print and save the output of the preprocess pipeline
	gold.save_files(preproc,pp.get_node('output'),datasink,("func","movement","struct","mask"),!noPrint)	
		
	preproc.write_graph(dotfilename=sequence+"-workflow"))
	return preproc




# check sequence
def check_sequence(opt_list,directory,seq):
	seq_dir = seq
	if seq == "preprocess" and "-"+seq in opt_list:
		return True	

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

########################################################################################
#
# run pipeline if used as standalone script
#
########################################################################################
if __name__ == "__main__":	
	opts = "[-dynamic_faces|-efnback|-reward|-resting_state|-asl|-fieldmap|-noprint]"
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
	if not os.path.exists(directory) and directory != "preprocess":
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
	
	if "-fieldmap" in opt_list:	
		useFieldmap = True
	if "-noprint" in opt_list:	
		noPrint = True
	
	if check_sequence(opt_list,directory,"reward"):
		log.info("\n\nREWARD pipeline ...\n\n")
		t = time.time()		
		reward = reward(directory,"reward")
		reward.run()		
		#reward.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
	if check_sequence(opt_list,directory,"resting_state"):
		log.info("\n\nRESTING pipeline ...\n\n")
		t = time.time()		
		resting = resting(directory,"resting_state")
		resting.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		#resting.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"efnback"):
		log.info("\n\nEFNBACK pipeline ...\n\n")
		t = time.time()		
		efnback = efnback(directory,"efnback")
		#efnback.run()		
		efnback.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
	if check_sequence(opt_list,directory,"dynamic_faces"):
		log.info("\n\nDynamic_Faces pipeline ...\n\n")
		t = time.time()		
		df = dynamic_faces(directory,"dynamic_faces")
		#df.run()		
		df.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"asl"):
		log.info("\n\nASL pipeline ...\n\n")
		t = time.time()		
		workflow = asl(directory,"asl")
		workflow.run()		
		#workflow.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
		
	log.info("\n\npipeline complete\n\n")
