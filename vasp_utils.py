import os
import numpy as np
from vasp import Vasp
from vasp.vasprc import VASPRC
from ase.io import read

VASPRC['mode'] = 'queue'
VASPRC['scheduler'] = 'OGE'
VASPRC['queue.command'] = 'qsub'
# Those variable are set for the demonstration only
# You may want to change them to suit your system
VASPRC['queue.q'] = 'tw'
VASPRC['queue.pe'] = 'orte'
VASPRC['queue.nprocs'] = '12'
VASPRC['queue.slot'] = 1
VASPRC['script']='vasp_script.sh'
VASPRC['validate'] = False
os.environ['VASP_PP_PATH']=os.path.abspath('examples/fakePotentials/')+'/'

def show_charge(calc):
    calc.bader()
    v=show_atoms(calc.atoms)
    v.add_label(color='black',scale=0.5,labelType='text',
                labelText=['%.2f'%c for c in calc.get_charges()],
                zOffset=1.2,attachment='middle_center')
    return v

def show_magmom(calc):
    v=show_atoms(calc.atoms)
    v.add_label(color='black',scale=0.5,labelType='text',
                labelText=['%.2f'%c for c in calc.get_magnetic_moments()],
                zOffset=1.2,attachment='middle_center')
    return v

def show_atoms(atoms):
    import nglview
    v = nglview.show_ase(atoms)
    v.clear_representations()
    v.add_unitcell()
    v.add_spacefill(radius_type='vdw',scale=0.5,
                    roughness=1,metalness=0)

    v.parameters = dict(clipDist=-100, sampleLevel=2)
    v.control.spin([1,0,0],3.14*1.5)
    return v
