import pandas as pd
import numpy as np
import os
import sys
import math
from enum import Enum
import igraph as ig

from typing import *

import warnings
from scipy.spatial.distance import pdist, squareform

from sklearn.preprocessing import MinMaxScaler
warnings.filterwarnings("ignore", category=RuntimeWarning)


def flatten_list(nested_list_arg: List[list]) -> List:
    """
    Flatten a nested list.
    turn [[1,2],[3,4]] to [1,2,3,4]
    """
    flat_list=[item for sublist in nested_list_arg for item in sublist]
    return flat_list

def split_strings(strings_list):
    split_list = []
    for string in strings_list:
        split_list.extend(string.split())
    return split_list

def get_df_from_file(filename,columns=['atom','x','y','z'],index=None):
    """
    Parameters
    ----------
    filename : str
        full file name to read.
    columns : str , optional
        list of column names for DataFrame. The default is None.
    splitter : str, optional
        input for [.split().] , for csv-',' for txt leave empty. The default is None.
    dtype : type, optional
        type of variables for dataframe. The default is None.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    with open(filename, 'r') as f:
        lines=f.readlines()[2:]
    splitted_lines=split_strings(lines)
    df=pd.DataFrame(np.array(splitted_lines).reshape(-1,4),columns=columns,index=index)
    df[['x','y','z']]=df[['x','y','z']].astype(float)
    return df


class XYZConstants(Enum):
    """
    Constants related to XYZ file processing
    """
    DF_COLUMNS=['atom','x','y','z']
    STERIMOL_INDEX = ['B1', 'B5', 'L', 'loc_B1','loc_B5']
    DIPOLE_COLUMNS = ['dip_x', 'dip_y', 'dip_z', 'total_dipole']
    RING_VIBRATION_COLUMNS = ['cross', 'cross_angle', 'para', 'para_angle']
    RING_VIBRATION_INDEX=['Product','Frequency','Sin_angle']
    VIBRATION_INDEX = ['Frequency', 'Amplitude']
    BONDED_COLUMNS = ['atom_1', 'atom_2', 'index_1', 'index_2']
    NOF_ATOMS = ['N', 'O', 'F']
    STERIC_PARAMETERS = ['B1', 'B5', 'L', 'loc_B1', 'loc_B5','RMSD']
    ELECTROSTATIC_PARAMETERS = ['dip_x', 'dip_y', 'dip_z', 'total_dipole','energy']

class GeneralConstants(Enum):
    """
    Holds constants for calculations and conversions
    1. covalent radii from Alvarez (2008) DOI: 10.1039/b801115j
    2. atomic numbers
    2. atomic weights
    """
    COVALENT_RADII= {
            'H': 0.31, 'He': 0.28, 'Li': 1.28,
            'Be': 0.96, 'B': 0.84, 'C': 0.76, 
            'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
            'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 
            'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
            'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 
            'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 
            'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 
            'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95,
            'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54,
            'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
            'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
            'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
            'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04,
            'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98,
            'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92,
            'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
            'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62,
            'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
            'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46,
            'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 
            'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
            'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87,
            'Am': 1.80, 'Cm': 1.69
    }
    
    BONDI_RADII={
        'H': 1.10, 'C': 1.70, 'F': 1.47,
        'S': 1.80, 'B': 1.92, 'I': 1.98, 
        'N': 1.55, 'O': 1.52, 'Co': 2.00, 
        'Br': 1.83, 'Si': 2.10,'Ni': 2.00,
        'P': 1.80, 'Cl': 1.75, 
    }
    CPK_RADII = {
    'C': 1.50,
    'C3': 1.60,
    'C6/N6': 1.70,
    'H': 1.00,
    'N': 1.50,
    'N4': 1.45,
    'O': 1.35,
    'O2': 1.35,
    'P': 1.40,
    'S': 1.70,
    'S1': 1.00,
    'F': 1.35,
    'Cl': 1.80,
    'S4': 1.40,
    'Br': 1.95,
    'I': 2.15,
    'X': 1.92,
    'F': 1.35
}
    # CPK_RADII={
    #     'C':1.50,   'H':1.00,   'S.O':1.70,  'Si':2.10,
    #     'C2':1.60,  'N':1.50,   'S1':1.00,   'Co':2.00,
    #     'C3':1.60,  'C66':1.70, 'F':1.35,    'Ni':2.00,
    #     'C4':1.50,  'N4':1.45,  'Cl':1.75,
    #     'C5/N5':1.70, 'O':1.35, 'S4':1.40,
    #     'C6/N6':1.70, 'O2':1.35, 'Br':1.95,
    #     'C7':1.70,    'P':1.40,  'I':2.15,
    #     'C8':1.50,    'S':1.70,  'B':1.92,
    
    # }

    REGULAR_BOND_TYPE = {

        'O.2': 'O', 'N.2': 'N', 'S.3': 'S',
        'O.3': 'O', 'N.1': 'N', 'S.O2': 'S',
        'O.co2': 'O', 'N.3': 'N', 'P.3': 'P',
        'C.1': 'C', 'N.ar': 'N',
        'C.2': 'C', 'N.am': 'N',
        "C.cat": 'C', 'N.pl3': 'N',
        'C.3': 'C', 'N.4': 'N',
        'C.ar': 'C', 'S.2': 'S',
    }

    BOND_TYPE={
        
        'O.2':'O2', 'N.2':'C6/N6','S.3':'S4',
        'O.3':'O', 'N.1':'N', 'S.O2':'S',
        'O.co2':'O', 'N.3':'C6/N6','P.3':'P',
        'C.1':'C', 'N.ar':'C6/N6',
        'C.2':'C3', 'N.am':'C6/N6',
        "C.cat":'C3', 'N.pl3':'C6/N6',
        'C.3':'C', 'N.4':'N4',
        'C.ar':'C6/N6', 'S.2':'S','H':'H' 
        }
    
    ATOMIC_NUMBERS ={
    '1':'H', '5':'B', '6':'C', '7':'N', '8':'O', '9':'F', '14':'Si',
             '15':'P', '16':'S', '17':'Cl', '35':'Br', '53':'I', '27':'Co', '28':'Ni'}
        

    ATOMIC_WEIGHTS = {
            'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,
            'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,
            'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,
            'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,
            'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,
            'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,
            'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,
            'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,
            'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,
            'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,
            'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,
            'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,
            'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,
            'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,
            'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,
            'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,
            'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,
            'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,
            'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,
            'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,
            'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,
            'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,
            'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,
            'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247
    }

import numpy.typing as npt


def adjust_indices(indices: npt.ArrayLike, adjustment_num: int=1) -> npt.ArrayLike:
    """
    adjust indices by adjustment_num
    """
    return np.array(indices)-adjustment_num

def adjust_indices_xyz(indices: npt.ArrayLike) -> npt.ArrayLike:
    """
    adjust indices by adjustment_num
    """
    return adjust_indices(indices, adjustment_num=1)

def calc_angle(p1: npt.ArrayLike, p2: npt.ArrayLike, degrees: bool=False) -> float: ###works, name in R: 'angle' , radians
    dot_product=np.dot(p1, p2)
    norm_p1=np.linalg.norm(p1)
    norm_p2=np.linalg.norm(p2)
    thetha=np.arccos(dot_product/(norm_p1*norm_p2))
    if degrees:
        thetha=np.degrees(thetha)   
    return thetha
    
def calc_new_base_atoms(coordinates_array: npt.ArrayLike, atom_indices: npt.ArrayLike):  #help function for calc_coordinates_transformation
    """
    a function that calculates the new base atoms for the transformation of the coordinates.
    optional: if the atom_indices is 4, the origin will be the middle of the first two atoms.
    """
    new_origin=coordinates_array[atom_indices[0]]
    if (len(atom_indices)==4):
        new_origin=(new_origin+coordinates_array[atom_indices[1]])/2
    new_y=(coordinates_array[atom_indices[-2]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[-2]]-new_origin))
    coplane=((coordinates_array[atom_indices[-1]]-new_origin)/np.linalg.norm((coordinates_array[atom_indices[-1]]-new_origin)+0.00000001))
    return (new_origin,new_y,coplane)

def np_cross_and_vstack(plane_1, plane_2):
    cross_plane=np.cross(plane_1, plane_2)
    united_results=np.vstack([plane_1, plane_2, cross_plane])
    return united_results

def calc_basis_vector(origin, y: npt.ArrayLike, coplane: npt.ArrayLike):#help function for calc_coordinates_transformation
    """
    origin: origin of the new basis
    y: y direction of the new basis
    coplane: a vector that is coplanar with the new y direction
    """
    # cross_y_plane=np.cross(coplane,y)
    # coef_mat=np.vstack([y, coplane, cross_y_plane])
    coef_mat=np_cross_and_vstack(coplane, y)
    angle_new_y_coplane=calc_angle(coplane,y)
    cop_ang_x=angle_new_y_coplane-(np.pi/2)
    # result_vector=[0,np.cos(cop_ang_x),0]
    result_vector=[np.cos(cop_ang_x), 0, 0]
    new_x,_,_,_=np.linalg.lstsq(coef_mat,result_vector,rcond=None)
    new_basis=np_cross_and_vstack(new_x, y)
    # new_z=np.cross(new_x,y)
    # new_basis=np.vstack([new_x, y, new_z])
    return new_basis

def transform_row(row_array, new_basis, new_origin, round_digits):
    translocated_row = row_array - new_origin
    return np.dot(new_basis, translocated_row).round(round_digits)



def calc_coordinates_transformation(coordinates_array: npt.ArrayLike, base_atoms_indices: npt.ArrayLike, round_digits:int=4 ,origin:npt.ArrayLike=None) -> npt.ArrayLike:#origin_atom, y_direction_atom, xy_plane_atom
    """
    a function that recives coordinates_array and new base_atoms_indices to transform the coordinates by
    and returns a dataframe with the shifted coordinates
    parameters:
    ----------
    coordinates_array: np.array
        xyz molecule array
    base_atoms_indices: list of nums
        indices of new atoms to shift coordinates by.
    origin: in case you want to change the origin of the new basis, middle of the ring for example. used in npa_df
    returns:
        transformed xyz molecule dataframe
    -------
        
    example:
    -------
    calc_coordinates_transformation(coordinates_array,[2,3,4])
    
    Output:
        atom       x       y       z
      0    H  0.3477 -0.5049 -1.3214
      1    B     0.0     0.0     0.0
      2    B    -0.0  1.5257     0.0
    """
    indices=adjust_indices_xyz(base_atoms_indices)
    new_basis=calc_basis_vector(*calc_new_base_atoms(coordinates_array,indices))    
    if origin is None:
        new_origin=coordinates_array[indices[0]]
    else:
        new_origin=origin

    transformed_coordinates = np.apply_along_axis(lambda x: transform_row(x, new_basis, new_origin, round_digits), 1,
                                                  coordinates_array)
    # transformed_coordinates=np.array([np.dot(new_basis,(row-new_origin)) for row in coordinates_array]).round(round_digits)
    return transformed_coordinates

def preform_coordination_transformation(xyz_df, indices=None):
    xyz_copy=xyz_df.copy()
    coordinates=np.array(xyz_copy[['x','y','z']].values)
    if indices is None:
        xyz_copy[['x','y','z']]=calc_coordinates_transformation(coordinates, [1,2,3])
    else:
        # print('indices', indices)
        xyz_copy[['x','y','z']]=calc_coordinates_transformation(coordinates, indices)
    # xyz_copy[['x','y','z']]=calc_coordinates_transformation(coordinates, get_indices([1,2,3])
    return xyz_copy

def calc_npa_charges(coordinates_array: npt.ArrayLike,charge_array: npt.ArrayLike,  geom_transform_indices=None):##added option for subunits
    """
    a function that recives coordinates and npa charges, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    
    Parameters
    ---------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    charge_array: np.array
        array of npa charges
    base_atoms_indices:list
        3/4 atom indices for coordinates transformation
        
    optional-sub_atoms:list
        calculate npa charges from a set of sub_atoms instead of all atoms.       
    Returns:
    -------
    dipole_df=calc_npa_charges(coordinates_array,charges,base_atoms_indices,sub_atoms)
    Output:
    dipole_df : pd.DataFrame
        output:            dip_x     dip_y     dip_z     total
                       0  0.097437 -0.611775  0.559625  0.834831
    """
    # what is NPA here??
    # indices=adjust_indices(base_atoms_indices)
    # transformed_coordinates=calc_coordinates_transformation(coordinates_array, indices)
    # if sub_atoms:
    #     atom_mask=sub_atoms
    # else:
    #     atom_mask=range(len(charge_array))
    # atom_mask=range(charge_array) if sub_atoms==None else sub_atoms
#TODO: Add option for sub_atoms!
    # Apply geometric transformation if specified
    # print(geom_transform_indices)
    if geom_transform_indices is not None:
        geometric_center = np.mean(coordinates_array[geom_transform_indices], axis=0)
        coordinates_array -= geometric_center

    dipole_xyz = np.vstack([(row[0] * row[1])for row in
                            list(zip(coordinates_array, charge_array))])
    dipole_vector=np.sum(dipole_xyz,axis=0)
    array_dipole=np.hstack([dipole_vector,np.linalg.norm(dipole_vector)])
    dipole_df=pd.DataFrame(array_dipole,index=XYZConstants.DIPOLE_COLUMNS.value).T
    # print(dipole_df)
    return dipole_df

def calc_dipole_gaussian(coordinates_array, gauss_dipole_array, base_atoms_indices ,geometric_transformation_indices=None):
    """
    a function that recives coordinates and gaussian dipole, transform the coordinates
    by the new base atoms and calculates the dipole in each axis
    """
    if geometric_transformation_indices:
        # Calculate the geometric center of specified indices
        geometric_center = np.mean(coordinates_array[geometric_transformation_indices], axis=0)
        # Translate all coordinates
        coordinates_array -= geometric_center

    indices=adjust_indices(base_atoms_indices)
    basis_vector=calc_basis_vector(*calc_new_base_atoms(coordinates_array, indices))
    gauss_dipole_array[0,0:3]=np.matmul(basis_vector,gauss_dipole_array[0,0:3])
    dipole_df=pd.DataFrame(gauss_dipole_array,columns=['dipole_x','dipole_y','dipole_z','total'])
    # print(geometric_transformation_indices, dipole_df)
    return dipole_df

def check_imaginary_frequency(info_df):##return True if no complex frequency, called ground.state in R
        bool_imaginary=not any([isinstance(frequency, complex) for frequency in info_df['Frequency']])
        return bool_imaginary



def indices_to_coordinates_vector(coordinates_array,indices):
    """
    a function that recives coordinates_array and indices of two atoms
    and returns the bond vector between them
    """

    if  isinstance(indices[0], tuple):
        bond_vector=[(coordinates_array[index[0]]-coordinates_array[index[1]]) for index in indices]
    else:
        bond_vector= coordinates_array[indices[0]]-coordinates_array[indices[1]]

    return bond_vector

def get_bonds_vector_for_calc_angle(coordinates_array,atoms_indices): ##for calc_angle_between_atoms

    indices=adjust_indices(atoms_indices)#three atoms-angle four atoms-dihedral
    augmented_indices=[indices[0],indices[1],indices[1],indices[2]]
    if len(indices)==4:
        augmented_indices.extend([indices[2],indices[3]])
    indices_pairs=list(zip(augmented_indices[::2],augmented_indices[1::2]))
  
    bond_vector=indices_to_coordinates_vector(coordinates_array,indices_pairs)
    return bond_vector
  

def calc_angle_between_atoms(coordinates_array,atoms_indices): #gets a list of atom indices
    """
    a function that gets 3/4 atom indices, and returns the angle between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    atoms_indices- list of ints
        a list of atom indices to calculate the angle between- [2,3,4]
   
    Returns
    -------
    angle: float
        the bond angle between the atoms
    """
    bonds_list=get_bonds_vector_for_calc_angle(coordinates_array,atoms_indices)
    if len(atoms_indices)==3:
      
        angle=calc_angle(bonds_list[0], bonds_list[1]*(-1), degrees=True)
    else:
        first_cross=np.cross(bonds_list[0],bonds_list[1]*(-1))
        second_cross=np.cross(bonds_list[2]*(-1),bonds_list[1]*(-1)) 
        angle=calc_angle(first_cross, second_cross, degrees=True)
    return angle

def get_angle_df(coordinates_array, atom_indices):
    """
    a function that gets a list of atom indices, and returns a dataframe of angles between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
    
    atom_indices- list of lists of ints
        a list of atom indices to calculate the angle between- [[2,3,4],[2,3,4,5]]
    """
 
    if isinstance(atom_indices, list) and all(isinstance(elem, list) for elem in atom_indices):
        indices_list=['angle_{}'.format(index) if len(index)==3 else 'dihedral_{}'.format(index) for index in atom_indices]
        angle_list=[calc_angle_between_atoms(coordinates_array,index) for index in atom_indices]
        return pd.DataFrame(angle_list,index=indices_list)
    else:
        indices_list=['angle_{}'.format(atom_indices) if len(atom_indices)==3 else 'dihedral_{}'.format(atom_indices)]
        angle=[calc_angle_between_atoms(coordinates_array,atom_indices)]
        return pd.DataFrame(angle,index=indices_list)


def calc_single_bond_length(coordinates_array,atom_indices):
    """
    a function that gets 2 atom indices, and returns the distance between thos atoms.
    Parameters
    ----------
    coordinates_array: np.array
        contains x y z atom coordinates
        
    atom_indices- list of ints
        a list of atom indices to calculate the distance between- [2,3]
   
    Returns
    -------
    distance: float
        the bond distance between the atoms
    """
    indices=adjust_indices(atom_indices)
    distance=np.linalg.norm(coordinates_array[indices[0]]-coordinates_array[indices[1]])
    return distance

def calc_bonds_length(coordinates_array,atom_pairs): 
    """
    a function that calculates the distance between each pair of atoms.
    help function for molecule class
    
    Parameters
    ----------
    coordinates_array: np.array
        xyz coordinates 
        
    atom_pairs : iterable
        list containing atom pairs-(([2,3],[4,5]))
        
    Returns
    -------
    pairs_df : dataframe
        distance between each pair 
        
        Output:
                               0
    bond length[2, 3]            1.525692
    bond length[4, 5]            2.881145
    
    """
    # Check the order of bond list and pairs
    bond_list=[calc_single_bond_length(coordinates_array,pair) for pair in atom_pairs]
    pairs=adjust_indices(atom_pairs)
    index=[('bond_length')+str(pair) for pair in pairs]
    pairs_df=pd.DataFrame(bond_list,index=index)
    return pairs_df



def direction_atoms_for_sterimol(bonds_df,base_atoms)->list: #help function for sterinol
    """
    a function that return the base atom indices for coordination transformation according to the bonded atoms.
    you can insert two atom indicess-[1,2] output [1,2,8] or the second bonded atom
    if the first one repeats-[1,2,1] output [1,2,3]
    """
    
    base_atoms_copy=base_atoms[0:2]
    origin,direction=base_atoms[0],base_atoms[1]
    bonds_df = bonds_df[~((bonds_df[0] == origin) & (bonds_df[1] == direction)) & 
                              ~((bonds_df[0] == direction) & (bonds_df[1] == origin))]
    
    try :
        base_atoms[2]==origin
        if(any(bonds_df[0]==direction)):
            # take the second atom in the bond where the first equeal to the direction, second option
            base_atoms_copy[2]=int(bonds_df[(bonds_df[0]==direction)][1].iloc[1])
        else:
            # take the first atom in the bond where the first equeal to the direction, second option
            base_atoms_copy[2]=int(bonds_df[(bonds_df[1]==direction)][0].iloc[1])
    except: 
        # if (any(bonds_df[0]==direction)):
        #     print(bonds_df[(bonds_df[0]==direction)][1].iloc[0],'h')
        #     base_atoms_copy.append(int(bonds_df[(bonds_df[0]==direction)][1].iloc[0])) if int(bonds_df[(bonds_df[0]==direction)][1].iloc[0])!=origin else base_atoms_copy.append(int(bonds_df[(bonds_df[0]==direction)][0].iloc[0]))
        # else:
        #     print(bonds_df[(bonds_df[1]==direction)][0].iloc[0],'w')
        #     base_atoms_copy.append(int(bonds_df[(bonds_df[1]==direction)][0].iloc[0])) if int(bonds_df[(bonds_df[1]==direction)][0].iloc[0])!=origin else base_atoms_copy.append(int(bonds_df[(bonds_df[1]==direction)][1].iloc[0]))
        for _, row in bonds_df.iterrows():
            if row[0] == direction:
                base_atoms_copy.append(row[1])
                break
            elif row[1] == direction:
                base_atoms_copy.append(row[0])
                break
    return base_atoms_copy

def get_molecule_connections(bonds_df,source,direction):
    graph=ig.Graph.DataFrame(edges=bonds_df,directed=True)
    paths=graph.get_all_simple_paths(v=source,mode='all')
    with_direction=[path for path in paths if (direction in path)]
    longest_path=np.unique(flatten_list(with_direction))
    return longest_path




def get_specific_bonded_atoms_df(bonds_df,longest_path,coordinates_df):
    """
    a function that returns a dataframe of the atoms that are bonded in the longest path.
    bonded_atoms_df: dataframe
        the atom type and the index of each bond.

       atom_1 atom_2 index_1 index_2
0       C      N       1       4
1       C      N       1       5
2       C      N       2       3
    """
    if longest_path is not None:
        edited_bonds_df=bonds_df[(bonds_df.isin(longest_path))].dropna().reset_index(drop=True)
    else:
        edited_bonds_df=bonds_df
    bonds_array=(np.array(edited_bonds_df)-1).astype(int) # adjust indices? 
    atom_bonds=np.vstack([(coordinates_df.iloc[bond]['atom'].values) for bond in bonds_array]).reshape(-1,2)
    bonded_atoms_df=(pd.concat([pd.DataFrame(atom_bonds),edited_bonds_df],axis=1))
    bonded_atoms_df.columns=[XYZConstants.BONDED_COLUMNS.value]
    return bonded_atoms_df


def remove_atom_bonds(bonded_atoms_df,atom_remove='H'):
    atom_bonds_array=np.array(bonded_atoms_df)
    delete_rows_left=np.where(atom_bonds_array[:,0]==atom_remove)[0] #itterrow [0] is index [1] are the values
    delete_rows_right=np.where(atom_bonds_array[:,1]==atom_remove)[0]
    atoms_to_delete=np.concatenate((delete_rows_left,delete_rows_right))
    new_bonded_atoms_df=bonded_atoms_df.drop((atoms_to_delete),axis=0)
    return new_bonded_atoms_df



def extract_connectivity(xyz_df, threshhold_distance=1.82):
    coordinates=np.array(xyz_df[['x','y','z']].values)
    atoms_symbol=np.array(xyz_df['atom'].values)
    # compute the pairwise distances between the points
    distances = pdist(coordinates)
    # convert the flat array of distances into a distance matrix
    dist_matrix = squareform(distances)
    dist_df=pd.DataFrame(dist_matrix).stack().reset_index()
    dist_df.columns = ['a1', 'a2', 'value']
    dist_df['first_atom']=[atoms_symbol[i] for i in dist_df['a1']]
    dist_df['second_atom']=[atoms_symbol[i] for i in dist_df['a2']]
    remove_list=[]
    dist_array=np.array(dist_df)
    for idx,row in enumerate(dist_array):
        if ((row[3]=='H') & (row[4] not in XYZConstants.NOF_ATOMS.value)):
            remove_list.append(idx)
        if ((row[3] == 'H') & (row[4] == 'H')):
            remove_list.append(idx)
        if (((row[3] == 'H') | (row[4] == 'H')) & (row[2]>=1.5) ):
            remove_list.append(idx)
        if ((row[2]>=threshhold_distance) | (row[2]==0)):
            remove_list.append(idx)
    dist_df=dist_df.drop(remove_list)
    dist_df=dist_df.drop_duplicates(subset=['value'])
    dist_array=np.array(dist_df[['a1','a2']])+1
    return pd.DataFrame(dist_array)

def get_center_of_mass(xyz_df):
    coordinates=np.array(xyz_df[['x','y','z']].values,dtype=float)
    atoms_symbol=np.array(xyz_df['atom'].values)
    masses=np.array([GeneralConstants.ATOMIC_WEIGHTS.value[symbol] for symbol in atoms_symbol])
    center_of_mass=np.sum(coordinates*masses[:,None],axis=0)/np.sum(masses)
    return center_of_mass

def get_closest_atom_to_center(xyz_df,center_of_mass):
    distances = np.sqrt((xyz_df['x'] - center_of_mass[0]) ** 2 + (xyz_df['y'] - center_of_mass[1]) ** 2 + (xyz_df['z'] - center_of_mass[2]) ** 2)
    idx_closest = np.argmin(distances)
    center_atom = xyz_df.loc[idx_closest]
    return center_atom

def get_sterimol_base_atoms(center_atom, bonds_df):
    # print(center_atom)
    center_atom_id=int(center_atom.name)+1
    base_atoms = [center_atom_id]
    if (any(bonds_df[0] == center_atom_id)):
        base_atoms.append(int(bonds_df[(bonds_df[0]==center_atom_id)][1].iloc[0]))
    else:
        base_atoms.append(int(bonds_df[(bonds_df[1]==center_atom_id)][0].iloc[0]))
    return base_atoms

def center_substructure(coordinates_array,atom_indices):
    atom_indices=adjust_indices(atom_indices)
    substructure=coordinates_array[atom_indices]
    center_substructure=np.mean(substructure,axis=0)
    return center_substructure

def nob_atype(xyz_df, bonds_df):
    
    symbols = xyz_df['atom'].values
    
    list_results=[]
    for index,symbol in enumerate(symbols):
        index+=1
        nob = bonds_df[(bonds_df[0] == index) | (bonds_df[1] == index)].shape[0]
        if symbol == 'H':
            result = 'H'
        elif symbol == 'F':
            result = 'F'
        elif symbol == 'P':
            result = 'P'
        elif symbol == 'Cl':
            result = 'Cl'
        elif symbol == 'Br':
            result = 'Br'
        elif symbol == 'I':
            result = 'I'
        elif symbol == 'O':
            if nob < 1.5:
                result = 'O2'
            elif nob > 1.5:
                result = 'O'
        elif symbol == 'S':
            if nob < 2.5:
                result = 'S'
            elif 2.5 < nob < 5.5:
                result = 'S4'
            elif nob > 5.5:
                result = 'S1'
        elif symbol == 'N':
            if nob < 2.5:
                result = 'C6/N6'
            elif nob > 2.5:
                result = 'N'
        elif symbol == 'C':
            if nob < 2.5:
                result = 'C3'
            elif 2.5 < nob < 3.5:
                result = 'C6/N6'
            elif nob > 3.5:
                result = 'C'
        else:
            result = 'X'

        list_results.append(result)

    return list_results

def get_sterimol_indices(coordinates,bonds_df):
    center = get_center_of_mass(coordinates)
    center_atom = get_closest_atom_to_center(coordinates,center)
    base_atoms = get_sterimol_base_atoms(center_atom,bonds_df)
    return base_atoms

def filter_atoms_for_sterimol(bonded_atoms_df, coordinates_df):
    """
    a function that filter out NOF bonds and H bonds and returns
     a dataframe of the molecule coordinates without them.
    """
    allowed_bonds_indices= pd.concat([bonded_atoms_df['index_1'],bonded_atoms_df['index_2']],axis=1).reset_index(drop=True)
    atom_filter=adjust_indices(np.unique([atom for sublist in allowed_bonds_indices.values.tolist() for atom in sublist]))
    edited_coordinates_df=coordinates_df.loc[atom_filter].reset_index(drop=True)
   
    return edited_coordinates_df


def get_extended_df_for_sterimol(coordinates_df, bonds_df, radii='CPK'):
    """
    A function that adds information to the regular coordinates_df

    Parameters
    ----------
    coordinates_df : dataframe
    bond_type : str
        The bond type of the molecule
    radii : str, optional
        The type of radii to use ('bondi' or 'CPK'), by default 'bondi'

    Returns
    -------
    dataframe
        The extended dataframe with additional columns

    """
    
    bond_type_map = GeneralConstants.BOND_TYPE.value
    ## if radius is cpk mapping should be done on atype, else on atom
    radii_map = GeneralConstants.CPK_RADII.value if radii == 'CPK' else GeneralConstants.BONDI_RADII.value
    
    df = coordinates_df.copy()  # make a copy of the dataframe to avoid modifying the original
    
    df['atype']=nob_atype(coordinates_df, bonds_df)
    
    
    # df['atype'] = df['atom'].map(bond_type_map).fillna(bond_type)
    df['magnitude'] = calc_magnitude_from_coordinates_array(df[['x', 'z']].astype(float))
    
    df['radius'] = df['atom'].map(radii_map)
    df['B5'] = df['radius'] + df['magnitude']
    df['L'] = df['y'] + df['radius']
    return df


def get_transfomed_plane_for_sterimol(plane,degree):
    """
    a function that gets a plane and rotates it by a given degree
    in the case of sterimol the plane is the x,z plane.
    Parameters:
    ----------
    plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
    degree : float
    """
    # print(degree,plane)
    cos_deg=np.cos(degree*(np.pi/180))
    sin_deg=np.sin(degree*(np.pi/180))
    rot_matrix=np.array([[cos_deg,-1*sin_deg],[sin_deg,cos_deg]])
    transformed_plane=np.vstack([np.matmul(rot_matrix,row) for row in plane]).round(4)
   # print("transformed plane")
    #print(transformed_plane)
    return transformed_plane


def calc_B1(transformed_plane,avs,edited_coordinates_df,column_index):
    """
    Parameters
    ----------
    transformed_plane : np.array
        [x,z] plane of the molecule coordinates.
        example:
            [-0.6868 -0.4964]
            [-0.7384 -0.5135]
            [-0.3759 -0.271 ]
            [-1.1046 -0.8966]
            [ 0.6763  0.5885]
    avs : list
        the max & min of the [x,z] columns from the transformed_plane.
        example:[0.6763, -1.1046, 0.5885, -0.8966
                 ]
    edited_coordinates_df : TYPE
        DESCRIPTION.
    column_index : int
        0 or 1 depending- being used for transformed plane.
    """
    
    ## get the index of the min value in the column compared to the avs.min
    idx=np.where(np.isclose(np.abs(transformed_plane[:,column_index]),(avs.min()).round(4)))[0][0]
    if transformed_plane[idx,column_index]<0:
        new_idx=np.where(np.isclose(transformed_plane[:,column_index],transformed_plane[:,column_index].min()))[0][0]
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[new_idx,column_index],
                                 transformed_plane[:,column_index]<=transformed_plane[new_idx,column_index]+1)
        
        transformed_plane[:,column_index]=-transformed_plane[:,column_index]
    else:
        bool_list=np.logical_and(transformed_plane[:,column_index]>=transformed_plane[idx,column_index]-1,
                                 transformed_plane[:,column_index]<=transformed_plane[idx,column_index])
        
    against,against_loc=[],[]
    B1,B1_loc=[],[]
    for i in range(1,transformed_plane.shape[0]): 
        if bool_list[i]:
            against.append(np.array(transformed_plane[i,column_index]+edited_coordinates_df['radius'].iloc[i]))
            against_loc.append(edited_coordinates_df['L'].iloc[i])
        if len(against)>0:
            B1.append(max(against))
            B1_loc.append(against_loc[against.index(max(against))])
            
        else:
            B1.append(np.abs(transformed_plane[idx,column_index]+edited_coordinates_df['radius'].iloc[idx]))
            B1_loc.append(edited_coordinates_df['radius'].iloc[idx])
            
    # print(f'B1: {B1}, B1_loc: {B1_loc}')      
    return [B1,B1_loc]

#find the direction

def b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane):
    
    """
    a function that gets a plane transform it and calculate the b1s for each degree.
    checks if the plane is in the x or z axis and calculates the b1s accordingly.
    Parameters:
    ----------
    extended_df : pd.DataFrame
    b1s : list
    b1s_loc : list
    degree_list : list
    plane : np.array
    """
    degree=[]
    for degree in degree_list:
        transformed_plane=get_transfomed_plane_for_sterimol(plane, degree)
        
        avs=np.abs([max(transformed_plane[:,0]),min(transformed_plane[:,0]), 
                    max(transformed_plane[:,1]),min(transformed_plane[:,1])])
        
        if min(avs) == 0:
            min_avs_indices = np.where(avs == min(avs))[0]
            if any(index in [0, 1] for index in min_avs_indices):
                tc = np.round(transformed_plane, 1)
                B1 = max(extended_df['radius'].iloc[np.where(tc[:, 0] == 0)])
                B1_loc = extended_df['L'].iloc[np.argmax(extended_df['radius'].iloc[np.where(tc[:, 0] == 0)])]
                b1s.append(B1)
                b1s_loc.append(B1_loc)
                continue  # Skip the rest of the loop

            elif any(index in [2, 3] for index in min_avs_indices):
                tc = np.round(transformed_plane, 1)
                B1 = max(extended_df['radius'].iloc[np.where(tc[:, 1] == 0)])
                B1_loc = extended_df['L'].iloc[np.argmax(extended_df['radius'].iloc[np.where(tc[:, 1] == 0)])]
                b1s.append(B1)
                b1s_loc.append(B1_loc)
                continue

        if np.where(avs==avs.min())[0][0] in [0,1]:
            B1,B1_loc=calc_B1(transformed_plane,avs,extended_df,0)
            
  
        elif np.where(avs==avs.min())[0][0] in [2,3]:
            B1,B1_loc=calc_B1(transformed_plane,avs,extended_df,1)
             
        
        b1s.append(np.unique(np.vstack(B1)).max())####check
        b1s_loc.append(np.unique(np.vstack(B1_loc)).max())
    

def get_b1s_list(extended_df, scans=90//5):
    
    b1s,b1s_loc=[],[]
    scans=scans
    degree_list=list(range(18,108,scans))
    plane=np.array(extended_df[['x','z']].astype(float))
    b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane)
    
    if b1s:
        try:
            back_ang=degree_list[np.where(b1s==min(b1s))[0][0]]-scans   
            front_ang=degree_list[np.where(b1s==min(b1s))[0][0]]+scans
            degree_list=range(back_ang,front_ang+1)
        except:
            print(np.where(np.isclose(b1s, min(b1s), atol=1e-8)))
            back_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]-scans
            front_ang=degree_list[np.where(np.isclose(b1s, min(b1s), atol=1e-8))[0][0]]+scans
            degree_list=range(back_ang,front_ang+1)
    else:
        print('no b1s found')
        return [np.array(b1s),np.array(b1s_loc)]
    # print(f'specific degree list: {degree_list}')
    b1s_for_loop_function(extended_df, b1s, b1s_loc, degree_list, plane)
    # print(f'b1 arrays: {[np.array(b1s),np.array(b1s_loc)]}')
    return [np.array(b1s),np.array(b1s_loc)]

def calc_sterimol(bonded_atoms_df,extended_df):
    edited_coordinates_df=filter_atoms_for_sterimol(bonded_atoms_df,extended_df)
    b1s,b1s_loc=get_b1s_list(edited_coordinates_df)
    B1=min(b1s[b1s>=0])
    loc_B1=max(b1s_loc[np.where(b1s[b1s>=0]==min(b1s[b1s>=0]))])
    B5=max(edited_coordinates_df['B5'].values)
    L=max(edited_coordinates_df['L'].values)
    loc_B5 = min(edited_coordinates_df['y'].iloc[np.where(edited_coordinates_df['B5'].values == B5)[0]])
    sterimol_df = pd.DataFrame([B1, B5, L, loc_B1,loc_B5], index=XYZConstants.STERIMOL_INDEX.value)
    return sterimol_df.T


def get_sterimol_df(coordinates_df, bonds_df, base_atoms, connected_from_direction, radii='bondi', sub_structure=True):

    bonds_direction = direction_atoms_for_sterimol(bonds_df, base_atoms)
    new_coordinates_df = preform_coordination_transformation(coordinates_df, bonds_direction)
    if sub_structure:
        if connected_from_direction is None:
            connected_from_direction = get_molecule_connections(bonds_df, base_atoms[0], base_atoms[1])
        else:
            connected_from_direction = connected_from_direction
    else:
        connected_from_direction = None
    bonded_atoms_df = get_specific_bonded_atoms_df(bonds_df, connected_from_direction, new_coordinates_df)
    
    extended_df = get_extended_df_for_sterimol(new_coordinates_df, bonds_df, radii)
    # print(f'extended dataframe: {extended_df}')
    ###calculations
   # print("new coords df")
   # print(new_coordinates_df)
    sterimol_df = calc_sterimol(bonded_atoms_df, extended_df)
    sterimol_df= sterimol_df.rename(index={0: str(base_atoms[0]) + '-' + str(base_atoms[1])})
    sterimol_df=sterimol_df.round(2)
    return sterimol_df


def calc_magnitude_from_coordinates_array(coordinates_array: npt.ArrayLike) -> List[float]:
    """
    Calculates the magnitudes of each row in the given coordinates array.

    Parameters
    ----------
    coordinates_array: np.ndarray
        A nx3 array representing the x, y, z coordinates of n atoms.

    Returns
    -------
    magnitudes: List[float]
        A list of magnitudes corresponding to the rows of the input array.

    """
    magnitude = np.linalg.norm(coordinates_array, axis=1)
    return magnitude



class Molecule():

    def __init__(self, molecule_xyz_filename):

   
        self.molecule_name = molecule_xyz_filename.split('.')[0]
        self.molecule_path = os.path.dirname(os.path.abspath(molecule_xyz_filename))
        os.chdir(self.molecule_path)
        
        
        self.xyz_df = get_df_from_file(molecule_xyz_filename)
        
        self.coordinates_array = np.array(self.xyz_df[['x', 'y', 'z']].astype(float))
        self.bonds_df = extract_connectivity(self.xyz_df)
        
        
            



    def process_sterimol_atom_group(self, atoms, radii):

        connected = get_molecule_connections(self.bonds_df, atoms[0], atoms[1])
        return get_sterimol_df(self.xyz_df, self.bonds_df, atoms, connected, radii)

    def get_sterimol(self, base_atoms: Union[None, Tuple[int, int]] = None, radii: str = 'bondi') -> pd.DataFrame:
        """
        Returns a DataFrame with the Sterimol parameters calculated based on the specified base atoms and radii.

        Args:
            base_atoms (Union[None, Tuple[int, int]], optional): The indices of the base atoms to use for the Sterimol calculation. Defaults to None.
            radii (str, optional): The radii to use for the Sterimol calculation. Defaults to 'bondi'.

        Returns:
            pd.DataFrame: A DataFrame with the Sterimol parameters.
            
            to add
            - only_sub- sterimol of only one part - i only have that.
            - drop some atoms.
        """
        if base_atoms is None:
            base_atoms = get_sterimol_indices(self.xyz_df, self.bonds_df)

        if isinstance(base_atoms[0], list):
            # If base_atoms is a list of lists, process each group individually and concatenate the results
            sterimol_list = [self.process_sterimol_atom_group(atoms, radii) for atoms in base_atoms]
            sterimol_df = pd.concat(sterimol_list, axis=0)

        else:
            # If base_atoms is a single group, just process that group
            sterimol_df = self.process_sterimol_atom_group(base_atoms, radii)
        return sterimol_df


    def swap_atom_pair(self, pair_indices: Tuple[int, int]) -> pd.DataFrame:
        """
        Swaps the positions of two atoms in the molecule and returns a new DataFrame with the updated coordinates.

        Args:
            pair_indices (Tuple[int, int]): The indices of the atoms to swap.

        Returns:
            pd.DataFrame: A new DataFrame with the updated coordinates.
        """
        pairs = adjust_indices(pair_indices)
        xyz_df = self.xyz_df
        temp = xyz_df.iloc[pairs[0]].copy()
        xyz_df.iloc[pairs[0]] = self.coordinates_array[pairs[1]]
        xyz_df.iloc[pairs[1]] = temp
        return xyz_df

  
 
    def get_coordination_transformation_df(self, base_atoms_indices: List[int]) -> pd.DataFrame:
        """
        Returns a new DataFrame with the coordinates transformed based on the specified base atoms.

        Args:
            base_atoms_indices (List[int]): The indices of the base atoms to use for the transformation.

        Returns:
            pd.DataFrame: A new DataFrame with the transformed coordinates.
        """
        new_coordinates_df = preform_coordination_transformation(self.xyz_df, base_atoms_indices)
        return new_coordinates_df

    ## not working after renumbering for some reason
    
    

    def get_bond_angle(self, atom_indices: List[int]) -> pd.DataFrame:
        """
        Returns a DataFrame with the bond angles calculated based on the specified atom indices.

        Args:
            atom_indices (List[int]): The indices of the atoms to use for the bond angle calculation.

        Returns:
            pd.DataFrame: A DataFrame with the bond angles.
        """
        return get_angle_df(self.coordinates_array, atom_indices)
    
    def get_bond_length_single(self, atom_pair):
        bond_length = calc_single_bond_length(self.coordinates_array, atom_pair)
        bond_length_df = pd.DataFrame([bond_length], index=[f'bond_length_{atom_pair[0]}-{atom_pair[1]}'])
        return bond_length_df

    def get_bond_length(self, atom_pairs):
        """
            Returns a DataFrame with the bond lengths calculated based on the specified atom pairs.

            Args:
                atom_pairs (Union[List[Tuple[int, int]], Tuple[int, int]]): The pairs of atoms to use for the bond length calculation.

            Returns:
                pd.DataFrame: A DataFrame with the bond lengths.
            """
        if isinstance(atom_pairs[0], list):
            # If atom_pairs is a list of lists, process each pair individually and concatenate the results
            bond_length_list = [self.get_bond_length_single(pair) for pair in atom_pairs]
            bond_df = pd.concat(bond_length_list, axis=0)
        else:
            # If atom_pairs is a single pair, just process that pair
            bond_df = self.get_bond_length_single(atom_pairs)
        return bond_df


    

class Molecules():
    
    def __init__(self,molecules_dir_name, renumber=False):
        self.molecules_path=os.path.abspath(molecules_dir_name)
        os.chdir(self.molecules_path) 
        self.molecules=[]
        self.failed_molecules=[]
        for file in os.listdir(): 
            if file.endswith('.xyz'):
                try:
                    self.molecules.append(Molecule(file))
                except:
                    self.failed_molecules.append(file)
                    print(f'Error: {file} could not be processed')
                    failed_file = file.rsplit('.feather', 1)[0] + '.feather_fail'
                # self.molecules.append(Molecule(log_file))
       
    


    def filter_molecules(self, indices):
        self.molecules = [self.molecules[i] for i in indices]
        self.molecules_names = [self.molecules_names[i] for i in indices]

    def get_sterimol_dict(self,base_atoms):
        sterimol_dict={}
        for molecule in self.molecules:
            try:
                sterimol_dict[molecule.molecule_name]=molecule.get_sterimol(base_atoms)
            except:
                print(f'failed to calculate for {molecule.molecule_name}')
                sterimol_dict[molecule.molecule_name]=np.nan

        return sterimol_dict
   

def main():
    import pprint
    
    os.chdir(r'/home/nati/Roni/Roni Alpha/Optimized_structures_xyz')
    mols=Molecules(os.getcwd())
    x=mols.get_sterimol_dict([4,7])
    pprint.pprint(x)
    
if __name__ == "__main__":
    main()