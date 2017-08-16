from dicom2nifti.convert_dicom import dicom_series_to_nifti
import matplotlib
import matplotlib.pyplot as plt
import nibabel
import numpy as np
import os
import requests

###########################
### API REQUESTS
###########################

def get_reports(user='wangkcc4', pw='pass', mrns=None):
    if region = 'test':
        stem=None#'https://montage.ynhh.org/api/v1'
    elif region = 'prod':':
        stem='https://montage.ynhh.org/api/v1'
    else:
        
        
    query_str = '&'.join(['format=json', 'exam_type=mri'])
    
    if mrns is not None:
        query_str += '&patient_mrn__in='+','.join(mrns)
        
    url = stem + '/report/?' + query_str
    
    r = requests.get(url, auth=(user, pw))
    
    return r.json()

def get_dcms(studies=None, region='test', user='wangkcc4', pw='pass'):
    if region = 'test':
        host='vnatest1vt'
        port='8083'
    elif region = 'prod':
        host='vnatest1vt'
        port='8083'
    else:
        
        
    query_str = '&'.join(['requestType=WADO','contenttype=image/jpeg'])
    
    if studies is not None:
        query_str = '/studies'+studies
        query_str = '/studies
    #'AcuoREST/search='
    #'0x00080060=ct' modality
    #'0x00100010=a*' names beginning with a
    #'&studyUID='
    #'&seriesUID='
    #'&objectUID='
    url = ''.join(['http://', host, ':', port, '/AcuoREST', query_str])
    
    r = requests.get(url, auth=(user, pw))
    
    return r.json()


###########################
### IMAGE FEATURES
###########################

def get_vol(img, dims, dim_units):
    return np.sum(img>0) * np.prod(dims)

def get_hist(img):
    """Returns histogram in array and graphical forms."""
    h = plt.hist(flatten(img, times=2))
    plt.title("Histogram")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    
    #mean_intensity = np.mean(diff[img > 0])
    #std_intensity = np.std(diff[img > 0])

    return h, plt.gcf()


###########################
### IMAGE PREPROCESSING
###########################

def dcm_load(path2series):
    """
    Load a dcm series as a 3D array along with its dimensions.
    
    returns as a tuple:
    - the normalized (0-255) image
    - the spacing between pixels in cm
    """
    tmp_fn = './tmp.nii'
    while os.path.isfile(tmp_fn):
        tmp_fn += '+'

    dicom_series_to_nifti(path2series, tmp_fn)
    ret = ni_load(tmp_fn)
    os.remove(tmp_fn)

    return ret

def ni_load(filename):
    """
    Load a nifti image as a 3D array along with its dimensions.
    
    returns as a tuple:
    - the normalized (0-255) image
    - the spacing between pixels in cm
    """
    
    img = nibabel.load(filename)
    img = nibabel.as_closest_canonical(img) # make sure it is in the correct orientation
    
    dims = img.header['pixdim'][1:4]
    dim_units = img.header['xyzt_units']
    
    img = np.asarray(img.dataobj).astype(dtype='float64')
    img = img[::-1,:,:] # Orient the image along the same axis as the binary masks.
    img =  (255 * (img / np.amax(img))).astype(dtype='uint8')
    
    if dim_units == 2: #or np.sum(img) * dims[0] * dims[1] * dims[2] > 10000:
        dims = [d/10 for d in dims]
    
    return img, dims

def apply_mask(img, mask_file):
    """Apply the mask in mask_file to img and return the masked image."""
    
    with open(mask_file, 'rb') as f:
        mask = f.read()
        mask = np.fromstring(mask, dtype='uint8')
        mask = np.reshape(mask, img.shape, order='F')
        mask = mask[:,::-1,:]
    
    img[mask <= 0] = 0
    
    return img

def create_diff(art_img, pre_img):
    diff = art_img - pre_img # Calculate the subtracted image.
    diff[diff < 0] = 0  # The pre-contrast should never be more than the arterial. Clamp negative values to zero.

    draw_fig(diff, 'whole')
    
    return diff
    
    
    
#########################
### UTILITY
#########################

def flatten(l, times=1):
    for _ in range(times):
        l = [item for sublist in l for item in sublist]
    return l

def save_img_array(img, filename, binary=True):
    """Save a numpy array to a binary file in .npy file format.

    In the future, may include options for """
    if binary:
        np.save(filename, img)

def draw_slice(img, filename, slice=None):
    """Draw a slice of an image of type np array and save it to disk."""
    cnorm = matplotlib.colors.Normalize(vmin=0, vmax=np.amax(img))
    
    if slice is None and len(img.shape)>2:
        slice=img.shape[2]//2

    w = 20
    h = int(float(img.shape[1]) / img.shape[0] * w)
    img_slice = img[:,:,slice]

    img_slice = np.rot90(img_slice)

    fig = plt.figure(frameon=False)
    fig.set_size_inches(w,h)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(img_slice, interpolation='bilinear', norm=cnorm, cmap=plt.cm.gray, aspect='auto')

    filename += '.png'
    plt.savefig(filename)
    print('Slice saved as %s' % filename)
    fig.set_size_inches(w//3,h//3)
    plt.show()