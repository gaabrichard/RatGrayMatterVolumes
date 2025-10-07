#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
   This script is adapted from scilpy-1.3.0 to be used as standalone
   Converts Varian .fdf files in a directory to .nii files.
   Does not extract bval and bvec. If you want this information, use scilpy.
   
   ex: fdf2nii.py /data/fsems_01.img /data/nii/fsems_01 --scale 10
"""

import argparse
import glob
import logging
import os
import re
import struct

import nibabel as nib
import numpy as np
import scipy

def load_fdf(file_path):
    """
    Load a Varian FDF file.

    Parameters
    ----------
    file_path: Path to the fdf file or directory

    Return
    ------
    data: FDF raw data
    header: Dictionary of the fdf header info.
    """
    if os.path.isdir(file_path):
        data, header = read_directory(file_path)
        dir_path = file_path
    else:
        data, header = read_file(file_path)
        dir_path = os.path.dirname(os.path.abspath(file_path))
        
    fid_dir_path = dir_path[:-3] + 'fid'

    
    procpar_path = os.path.join(fid_dir_path, 'procpar')

    if os.path.exists(procpar_path):
        add_orientation_info(procpar_path, header)
    else:
        logging.warning('Could not find the procpar file {0}. \n'
                        'Image orientation may be incorrect'
                        'you need to have a procpar file in the fdf directory.')
                        

    return data, header
    
def add_orientation_info(procpar_path, header):
    """
    Extract orientation information from the procpar file
    and puts it in the header. Adapted from NWM2 code by Réjean Lebel
    and Jérémie Fouquet

    Parameters
    ----------
    procpar_path: Path to the procpar file.
    header: Header to add orientation info.

    Return
    ------
    None
    """
    
    #Read procpar
    with open(procpar_path, 'r') as procpar:
        for line in procpar:
            if line.startswith('acqdim'):
                dim = next(procpar)
                dim=[float(x) for x in dim.split()[1:]]
            if line.startswith('fn '):
                fn = next(procpar)
                fn=[float(x) for x in fn.split()[1:]]
            if line.startswith('fn1 '):
                fn1 = next(procpar)
                fn1=[float(x) for x in fn1.split()[1:]]
            if line.startswith('fn2 '):
                fn2 = next(procpar)
                fn2=[float(x) for x in fn2.split()[1:]]
            if line.startswith('gap '):
                gap = next(procpar)
                gap=[float(x) for x in gap.split()[1:]]
            if line.startswith('lpe '):
                lpe = next(procpar)
                lpe=[float(x) for x in lpe.split()[1:]]
            if line.startswith('lpe2 '):
                lpe2 = next(procpar)
                lpe2=[float(x) for x in lpe2.split()[1:]]
            if line.startswith('lro '):
                lro = next(procpar)
                lro=[float(x) for x in lro.split()[1:]]
            if line.startswith('ns '):
                ns = next(procpar)
                ns=[float(x) for x in ns.split()[1:]]
            if line.startswith('phi '):
                phi = next(procpar)
                phi=[float(x) for x in phi.split()[1:]]
            if line.startswith('ppe '):
                ppe = next(procpar)
                ppe=[float(x) for x in ppe.split()[1:]]
            if line.startswith('ppe2 '):
                ppe2 = next(procpar)
                ppe2=[float(x) for x in ppe2.split()[1:]]
            if line.startswith('pro '):
                pro = next(procpar)
                pro=[float(x) for x in pro.split()[1:]]
            if line.startswith('psi '):
                psi = next(procpar)
                psi=[float(x) for x in psi.split()[1:]]
            if line.startswith('pss '):
                pss = next(procpar)
                pss=[float(x) for x in pss.split()[1:]]
            if line.startswith('theta '):
                theta = next(procpar)
                theta=[float(x) for x in theta.split()[1:]]
            if line.startswith('thk '):
                thk = next(procpar)
                thk=[float(x) for x in thk.split()[1:]]
    
    #Get matrix size and voxel size           
    
    if dim[0] == 2:
       recon_voxel_size = [fn[0]/2, fn1[0]]
       recon_voxel_size = recon_voxel_size + ns
       voxel_d3 = thk[0] + gap[0]
       translation_d3 = min(pss)*10
    elif dim[0] == 3:
        recon_voxel_size = [fn[0]/2, fn1[0]/2]
        recon_voxel_size = recon_voxel_size + [fn2[0]/2]
        voxel_d3 = lpe2[0]*10/recon_voxel_size[2]
        translation_d3 = ppe2[0]*10 - (lpe2[0]*10 - voxel_d3)/2

        
    voxel_size = [lro[0]*10/recon_voxel_size[0], lpe[0]*10/recon_voxel_size[1], voxel_d3]
    

    
    #Compute translation in all directions
    translation = [pro[0]*10 - (lro[0]*10 - voxel_size[0])/2, \
                  ppe[0]*10 -(lpe[0]*10 - voxel_size[1])/2, \
                  translation_d3]
                  
    #Create scaling and translation matrix
    voxel_size.append(1)
    scaling_affine = np.diag(voxel_size)
    translation_affine = np.eye(4)
    translation_affine[0:3,3] = translation

         
    #Create rotation matrix
    rotation_affine = np.eye(4)
    
    psi = np.radians(psi[0])
    phi = np.radians(phi[0])
    theta = np.radians(theta[0])
    
    if dim[0] == 3:
       psi = psi + np.pi/2
       theta = theta - np.pi/2 #Not sure it works for all cases (maybe just a switch between psi and theta?)
    elif dim[0] == 2: #Again, grasping at straws here
       phi = phi - np.pi
       theta = theta + np.pi/2

       
    rot_theta = [[np.cos(theta), 0, np.sin(theta)], [0,1,0], [-np.sin(theta), 0, np.cos(theta)]]
    rot_psi = [[np.cos(psi), -np.sin(psi),0], [np.sin(psi), np.cos(psi), 0], [0,0,1]]
    rot_phi = [[1,0,0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]]
    
    
    rot_phi = np.array(rot_phi)
    rot_theta = np.array(rot_theta)
    rot_psi = np.array(rot_psi)
    
    
    
    rotation_affine[0:3,0:3] = np.dot(rot_psi, np.dot(rot_theta, rot_phi))
 
    
    affine_matrix = np.dot(rotation_affine, np.dot(translation_affine, scaling_affine))
    header['affine'] = affine_matrix
    
    #Rotate data to fit MATLAB scripts
    
    if dim[0]==3:
       header['dim']=3
    elif dim[0]==2:
       header['dim']=2
       
    
    return header
    
    
def read_file(file_path):
    """
    Read a single fdf file.

    Parameters
    ----------
    path: Path to the fdf file

    Return
    ------
    raw_header: Dictionary of header information
    data: Numpy array of the fdf data
    """
    raw_header = dict()

    raw_header['shape'] = [-1, -1, 1, 1]
    raw_header['endian'] = '>'

    # Extracts floating point numbers
    # ex: 'float  roi[] = {3.840000,3.840000,0.035000};'
    # returns ['3.840000', '3.840000', '0.035000']
    float_regex = '[-+]?[0-9]*[.]?[0-9]+'

    # Extracts value of a line of the type:
    # 'int    slice_no = 1;' would return '1'
    named_value_regex = '= *\"*(.*[^\"])\"* *;'

    # (tag_in_file, tag_in_header)
    find_values = (('echos', 'nechoes'),
                   ('echo_no', 'echo_no'),
                   ('nslices', 'nslices'),
                   ('slice_no', 'sl'),
                   ('bigendian', 'endian'),
                   ('array_dim', 'array_dim'),
                   ('bigendian', 'endian'),
                   ('studyid', 'studyid'))

    with open(file_path, 'rb') as fp:
        # Read entire file
        while True:
            line = fp.readline()
            line = line.decode()

            if line[0] == chr(12):
                break

            # Check line for tag, extract value with the regex then put it in
            # the header with the associated tag.
            for file_key, head_key in find_values:
                if line.find(file_key) > 0:
                    raw_header[head_key] = \
                        re.findall(named_value_regex, line)[0]
                    break

            if type(raw_header['endian']) is int:
                raw_header['endian'] = \
                    '>' if int(raw_header['endian']) != 0 else '<'

            if line.find('abscissa') > 0:
                # Extracts units in quotes.
                # ex: 'char  *abscissa[] = {"cm", "cm"}' returns
                # ["cm", "cm"]
                m = re.findall('\"[a-z]{2}\"', line.rstrip())

                unit = m[0].strip('"')

                # We convert everything in mm
                # Nifti doesn't support 'cm' anyway...
                if unit == 'cm':
                    unit = 'mm'

                raw_header['xyz_units'] = unit
                raw_header['t_units'] = 'unknown'

            elif line.find('roi') > 0:
                m = re.findall(float_regex, line.rstrip())
                raw_header['real_voxel_dim'] = \
                    np.array([float(x)*10 for x in m])

            elif line.find('orientation') > 0:
                m = re.findall(float_regex, line.rstrip())
                raw_header['orientation'] = np.array([float(x) for x in m])

            elif line.find('origin') > 0:
                m = re.findall(float_regex, line.rstrip())
                raw_header['origin'] = \
                    np.array([float(x) for x in m])

            elif line.find('matrix') > 0:
                # Extracts digits.
                # ex: 'float  matrix[] = {128, 128};'
                # returns ['128', '128']
                m = re.findall(r'\d+', line.rstrip())
                raw_header['shape'] = np.array([int(x) for x in m])

        # Total number of data pixels
        # nb_voxels = reduce(operator.mul, raw_header['shape'])
        nb_voxels = np.prod(raw_header['shape'])

        # Set how data is packed
        raw_header['fmt'] = "{}f".format(nb_voxels)
        if '<' in raw_header['endian']:
            raw_header['fmt'] = '<'+raw_header['fmt']
        else:
            raw_header['fmt'] = '>'+raw_header['fmt']

        # Go to the beginning of the data segment
        fp.seek(-nb_voxels * 4, 2)
        data = struct.unpack(raw_header['fmt'], fp.read(nb_voxels*4))

    # Get correct voxel dimensions in mm
    raw_header['voxel_dim'] = \
        [j/i for i, j in zip(raw_header['shape'],
                             raw_header['real_voxel_dim'])]

    correct_shape = raw_header['shape'][::-1]

    # Reshape the data according to image dimensions
    data = np.array(data).reshape(correct_shape).squeeze()

    if len(raw_header['shape']) != 2:
        data = data.transpose(2, 1, 0)

    data = np.rot90(data, 3)[:, ::-1]

    return raw_header, data


def read_directory(path):
    """
    Parameters
    ----------
    Read a directory containing multiple 2D ``.fdf`` files. The method
    should return ``None`` if the directory is empty.

    path: Path to the input directory

    Return
    ------
    data: Numpy array containing data
    final_header: Header information
    """
    files = glob.glob(os.path.join(path, '*.fdf'))
    files.sort()


    all_headers, all_data = zip(*[read_file(fl) for fl in files])

    if not all_headers:
        return None

    final_header = all_headers[0]

    # Fix data axis
    all_data = np.array(all_data)
    
    if len(all_data.shape) < 4:
        all_data = np.transpose(all_data, (1, 2, 0))

        # Set real shape
        final_header['shape'] = all_data.shape

        # Correct voxel dimensions to fit data shape
        if len(final_header['shape']) != len(final_header['voxel_dim']):
            final_header['voxel_dim'].extend(
                final_header['real_voxel_dim'][len(final_header['shape'])-1:])

        # Support for fourth dimension
        time = int(float(final_header['array_dim']))
        if time > 1:
            final_header['shape'] = (final_header['shape'][0],
                                     final_header['shape'][1],
                                     round(final_header['shape'][2] / time),
                                     time)

            final_header['voxel_dim'].append(1)

            all_data = all_data.reshape(final_header['shape'])
    else:
        all_data = np.transpose(all_data, (3, 2, 1, 0))
        all_data = np.transpose(all_data, (1, 0, 2, 3))
        final_header['shape'] = all_data.shape
        final_header['voxel_dim'].append(1.0)

    return all_data, final_header
    
def format_raw_header(header):
    """
    Format the header to a Nifti1Image format.

    Parameters
    ----------
    header: Raw dictionary of header information

    Return
    ------
    nifti1_header: Header to save in the nifti1 file.
    """
    if header is None:
        return header

    nifti1_header = nib.nifti1.Nifti1Header()

    nifti1_header.set_data_shape(header['shape'])
    nifti1_header.set_xyzt_units(header['xyz_units'], header['t_units'])
    nifti1_header.set_data_dtype('float32')
    

    return nifti1_header


def save_babel(data, header, out_path, scale, affine=None):
    """
    Save a loaded fdf file to nifti.

    Parameters
    ----------
    out_path: Path of the nifti file to be saved
    data: Raw data to be saved
    raw_header: Raw header from fdf files
    affine: Affine transformation to save with the data

    Return
    ------
    None
    """
    nifti1_header = format_raw_header(header)

    if 'affine' in header:    
        affine = header['affine']

    nifti1_header.set_data_shape(data.shape)

    img = nib.nifti1.Nifti1Image(dataobj=data,
                                 header=None,
                                 affine=affine)
    
    vox_dim = [round(num, 3)*scale for num in header['voxel_dim'][0:4]]
    img.header.set_zooms(vox_dim)
    qform = img.header.get_qform()
    #print(qform)

    if 'origin' in nifti1_header:
        qform[:len(nifti1_header['origin']), 3] = -nifti1_header['origin']

    img.header.set_qform(qform)
    img.update_header()
    img.to_filename(out_path)
    
    
def build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,)

    p.add_argument('in_path',
                   help='Path to the FDF file or folder to convert.')
    p.add_argument('out_path',
                   help='Path to the nifti file to write on disk.')
    p.add_argument('--scale', type=float, default=1.0,
                   help='Rescale voxel size by the specified factor default=1')

    return p


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    data, header = load_fdf(args.in_path)
    
    if header['dim']==3:
       data = np.rot90(data, axes=(1,2))
       data = np.rot90(data, k=2, axes=(0,1))

    save_babel(data, header,
               args.out_path,
               scale=args.scale)


if __name__ == "__main__":
    main()    
