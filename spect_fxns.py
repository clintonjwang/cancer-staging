import helper_fxns as hf
import nibabel as nib
import numpy as np
import os
import pandas as pd
import pyelastix
import SimpleITK as sitk
import shutil


###########################
### ???
###########################

def transform_niis(in_img_path, transform_paths, patient_nii_paths, target_imgs=None, out_img_path=None, overwrite=True):
	"""Transforms image based on previous transform and scaling to target_dims."""
	
	if out_img_path is None:
		out_img_path = add_to_filename(in_img_path, "-reg")
	if (not overwrite) and os.path.exists(out_img_path):
		print(out_img_path, "already exists. Skipping transform.")
		return False
		
	temp_path = add_to_filename(in_img_path, "-temp")
		
	shutil.copyfile(in_img_path, temp_path)
	
	for index, transform_path in enumerate(transform_paths):
		if target_imgs is None:
			target_img_type = transform_path[transform_path.find('_')+1:transform_path.rfind('_')]
			target_img, _ = hf.ni_load(patient_nii_paths[target_img_type])
		else:
			target_img = target_imgs[index]

		transform(temp_path, transform_path)
		
		save_nii(hf.rescale(hf.ni_load(temp_path)[0], target_img.shape)[0], temp_path)
	
	if temp_path != out_img_path:
		if os.path.exists(out_img_path):
			os.remove(out_img_path)
		os.rename(temp_path, out_img_path)
	
	return True

def check_dcm_paths(patient_id, patient_df):
	blmri_path = "Z:\\Isa\\blmri\\"+str(patient_df.loc[patient_id, "BL-MRI"])
	fumri_path = "Z:\\Isa\\fumri\\"+str(patient_df.loc[patient_id, "FU1/2-MRI"])
	if not os.path.exists(blmri_path+"\\T1_AP"):
		print(blmri_path+"\\T1_AP")
	if not os.path.exists(blmri_path+"\\T1_BL"):
		print(blmri_path+"\\T1_AP")
	if not os.path.exists(fumri_path+"\\T1_AP"):
		print(fumri_path+"\\T1_AP")
	if not os.path.exists(fumri_path+"\\T1_BL"):
		print(fumri_path+"\\T1_AP")

def save_nii(img, dest, dims=(1,1,1)):
	affine = np.eye(4)
	for i in range(3):
		affine[i,i] = dims[i]
	if len(img.shape) == 4:
		nii = nib.Nifti1Image(img[:,::-1,:,:], affine)
	else:
		nii = nib.Nifti1Image(img[:,::-1,:], affine)
	nib.save(nii, dest)

def is_segmented(patient_id):
	if os.path.exists(base_dir + "BL-segs"):
		print(patient_id)

def reg_niis(fixed_img_type, moving_img_type, patient_nii_paths, fixed_img=None,
			 moving_img=None, out_img_path=None, out_transform_path=None, overwrite=False):
	"""Registers images. 
	"""
	
	if out_img_path is None:
		out_img_path = add_to_filename(patient_nii_paths[moving_img_type], "-reg")
		
	if out_transform_path is None:
		out_transform_path = patient_nii_paths['base'] + moving_img_type+'_'+fixed_img_type+'_transform.hdf5'
		
	temp_path = add_to_filename(patient_nii_paths[moving_img_type],"-scaled")
		
	if (not overwrite) and os.path.exists(out_img_path):
		print(out_img_path, "already exists. Skipping registration.")
		return None
	
	if fixed_img is None:
		fixed_img, _ = hf.ni_load(patient_nii_paths[fixed_img_type])
	if moving_img is None:
		moving_img, _ = hf.ni_load(patient_nii_paths[moving_img_type])
		
	if moving_img.shape == fixed_img.shape:
		print(fixed_img_type, "and", moving_img_type, "have the same shape. Skipping registration.")
		return None
	
	moving_img_scaled, _ = hf.rescale(moving_img, fixed_img.shape)
	save_nii(moving_img_scaled, temp_path)
	
	reg_img(patient_nii_paths[fixed_img_type], temp_path, out_transform_path, out_img_path)
	
	os.remove(temp_path)
	
	return out_img_path, out_transform_path

def set_paths(patient_id, patient_df, dcm_paths, nii_paths, mask_paths):
	if not os.path.exists("Z:\\Isa\\"+str(patient_id)):
		os.makedirs("Z:\\Isa\\"+str(patient_id))
		
	base_dir = "Z:\\Isa\\"+str(patient_id)+"\\"

	dcm_paths[patient_id] = {"ct": "Z:\\Isa\\spect\\"+str(patient_df.loc[patient_id, "SPECT"])+"\\CT",
			"fused": "Z:\\Isa\\spect\\"+str(patient_df.loc[patient_id, "SPECT"])+"\\Fused",
			"spect": "Z:\\Isa\\spect\\"+str(patient_df.loc[patient_id, "SPECT"])+"\\SPECT",
			"fumri-art": "Z:\\Isa\\fumri\\"+str(patient_df.loc[patient_id, "FU1/2-MRI"])+"\\T1_AP",
			"fumri-pre": "Z:\\Isa\\fumri\\"+str(patient_df.loc[patient_id, "FU1/2-MRI"])+"\\T1_BL",
			"blmri-art": "Z:\\Isa\\blmri\\"+str(patient_df.loc[patient_id, "BL-MRI"])+"\\T1_AP",
			"blmri-pre": "Z:\\Isa\\blmri\\"+str(patient_df.loc[patient_id, "BL-MRI"])+"\\T1_BL"}

	nii_paths[patient_id] = {'base': base_dir}
	img_types = ['blmri-art', 'blmri-pre', 'fumri-art', 'fumri-pre', 'ct', 'spect', 'fused', 'fused-ch1', 'fused-ch2', 'fused-ch3']
	for img_type in img_types:
		nii_paths[patient_id][img_type] = base_dir + img_type + ".nii"

	mask_paths[patient_id] = {"liver": base_dir + "BL-segs\\BL-Liver.ids",
				"tumor": base_dir + "BL-segs\\BL-Tumor.ids",
				"necrosis-bl": base_dir + "BL-segs\\necrosis.ids",
				"viable-tumor-bl": base_dir + "BL-segs\\viable_tumor.ids",
				"necrosis-fu": base_dir + "FU-segs\\necrosis.ids",
				"viable-tumor-fu": base_dir + "FU-segs\\viable_tumor.ids",
				"necrosis-fu-scaled": base_dir + "FU-segs\\necrosis-scaled.ids",
				"viable-tumor-fu-scaled": base_dir + "FU-segs\\viable_tumor-scaled.ids"}

	if not os.path.exists(dcm_paths[patient_id]["spect"]):
		spect_path = "Z:\\Isa\\spect\\"+str(patient_df.loc[patient_id, "SPECT"])
		for fn in os.listdir(spect_path):
			if 'recon - ac' in fn:
				dcm_paths[patient_id]['spect'] = spect_path + "\\" + fn
				break
				
	if not os.path.exists(dcm_paths[patient_id]["ct"]):
		spect_path = "Z:\\Isa\\spect\\"+str(patient_df.loc[patient_id, "SPECT"])
		for fn in os.listdir(spect_path):
			if 'y90 sirs' in fn and 'ac' not in fn:
				dcm_paths[patient_id]['ct'] = spect_path + "\\" + fn
				break
				
	if not os.path.exists(dcm_paths[patient_id]["fused"]):
		spect_path = "Z:\\Isa\\spect\\"+str(patient_df.loc[patient_id, "SPECT"])
		for fn in os.listdir(spect_path):
			if 'fused tran' in fn:
				dcm_paths[patient_id]['fused'] = spect_path + "\\" + fn
				break

def save_niis(patient_id, dcm_paths, nii_paths, overwrite=True):
	if (not overwrite) and os.path.exists(nii_p['spect']):
		return

	nii_p = nii_paths[patient_id]

	spect_img = hf.get_spect_series(dcm_paths[patient_id]['spect'])
	save_nii(spect_img, nii_p['spect'])

	fused_img = hf.get_spect_series(dcm_paths[patient_id]['fused'])
	save_nii(fused_img, nii_p['fused'])
	save_nii(fused_img[:,:,:,0], nii_p['fused-ch1'])
	save_nii(fused_img[:,:,:,1], nii_p['fused-ch2'])
	save_nii(fused_img[:,:,:,2], nii_p['fused-ch3'])

	ct_img, dims = hf.dcm_load(dcm_paths[patient_id]['ct'])
	save_nii(ct_img, nii_p['ct'], dims)

	blmri_art, dims = hf.dcm_load(dcm_paths[patient_id]['blmri-art'])
	save_nii(blmri_art, nii_p['blmri-art'], dims)

	blmri_pre, dims = hf.dcm_load(dcm_paths[patient_id]['blmri-pre'])
	save_nii(blmri_pre, nii_p['blmri-pre'], dims)

	fumri_art, dims = hf.dcm_load(dcm_paths[patient_id]['fumri-art'])
	save_nii(fumri_art, nii_p['fumri-art'], dims)

	fumri_pre, dims = hf.dcm_load(dcm_paths[patient_id]['fumri-pre'])
	save_nii(fumri_pre, nii_p['fumri-pre'], dims)


###########################
### SimpleITK
###########################

def reg_img(fixed_path, moving_path, out_transform_path, out_img_path, verbose=False):
	"""Assumes fixed and moving images are the same dimensions"""

	fixed = sitk.ReadImage(fixed_path, sitk.sitkFloat32)
	moving = sitk.ReadImage(moving_path, sitk.sitkFloat32)

	matcher = sitk.HistogramMatchingImageFilter()
	matcher.SetNumberOfHistogramLevels(1024)
	matcher.SetNumberOfMatchPoints(7)
	matcher.ThresholdAtMeanIntensityOn()
	moving = matcher.Execute(moving,fixed)

	demons = sitk.DemonsRegistrationFilter()
	demons.SetNumberOfIterations( 50 )
	demons.SetStandardDeviations( 1.0 )

	if verbose:
		def command_iteration(filter):
			print("{0:3} = {1:10.5f}".format(filter.GetElapsedIterations(), filter.GetMetric()))
		demons.AddCommand( sitk.sitkIterationEvent, lambda: command_iteration(demons) )

	displacementField = demons.Execute( fixed, moving )

	outTx = sitk.DisplacementFieldTransform( displacementField )
	sitk.WriteTransform(outTx, out_transform_path)
	sitk.WriteImage(moving, out_img_path)


def transform(moving_path, transform_path, target_path=None):
	"""Transforms without scaling image"""
	
	if target_path is None:
		target_path = moving_path
	
	moving = sitk.ReadImage(moving_path, sitk.sitkFloat32)
	tx = sitk.ReadTransform(transform_path)
	moving_reg = sitk.Resample(moving, tx)
	sitk.WriteImage(moving_reg, target_path)

###########################
### MISC
###########################

def add_to_filename(fn, addition):
	x = fn.find(".")
	return fn[:x] + addition + fn[x:]