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


def show_traj(traj):
    # traj can be list of atoms
    import nglview
    v = nglview.show_asetraj(traj)
    v.clear_representations()
    # it seems that nglview traj cannot handle unit cell right now
    #v.add_unitcell()
    v.add_spacefill(radius_type='vdw',scale=0.5,
                    roughness=1,metalness=0)
    v.parameters = dict(clipDist=-100, sampleLevel=2)
    return v

def show_vib(calc,mode=-1,amplitude=1,info=True):
    freqs,modes = calc.get_vibrational_modes()
    displace = np.concatenate((
        np.linspace(0, -1, 5, endpoint=False),
        np.linspace(-1, 1, 10, endpoint=False),
        np.linspace(1, 0, 5, endpoint=False)
        ))
    traj=[]
    n_img = sum([np.iscomplexobj(freq) for freq in freqs])
    n_mode = len(freqs)
    n_show=(mode+n_mode)%n_mode+1

    msg='%i modes in total, %i imaginary, showing the %ith,'%(n_mode,n_img,n_show)
    if np.iscomplexobj(freqs[mode]):
        msg+='freq is %.2fi(meV).'%(freqs[mode].real*1000)
    else:
        msg+='freq is %.2f(meV).'%(freqs[mode]*1000)
    if info:
        print(msg)
    for i in displace:
        image = calc.atoms.copy()
        image.set_positions(image.positions+i*modes[mode]*amplitude)
        traj.append(image)
    return show_traj(traj)


def form_traj(atoms,mode,amplitude):
    traj=[]
    displace = np.concatenate((
        np.linspace(0, -1, 5, endpoint=False),
        np.linspace(-1, 1, 10, endpoint=False),
        np.linspace(1, 0, 5, endpoint=False)))
    for i in displace:
        image = atoms.copy()
        image.set_positions(image.positions+i*mode*amplitude)
        traj.append(image)
    return traj

def show_vibs(calc,amplitude=1):
    import traitlets
    import nglview
    import ipywidgets as widgets
    from ipywidgets import Layout,HBox,VBox

    trajs=[]
    freqs,modes=calc.get_vibrational_modes()
    freq_list = ['-%.2f'%(f.real*1000) if np.iscomplexobj(f)
             else '%.2f'%(f*1000) for f in freqs]

    for n,mode in enumerate(modes):
            trajs.append([a.positions for a in form_traj(calc.atoms,mode,amplitude)])

    v=nglview.show_asetraj([calc.atoms])
    v.clear_representations()
    v.add_spacefill(radius_type='vdw',scale=0.5,
                    roughness=1,metalness=0)
    v._set_size('450px','300px')
    v.camera='orthographic'
    v.parameters=dict(clipDist=0,sampleLevel=1)

    select = widgets.Select(
        options=freq_list,
        value=freq_list[-1],
        layout=Layout(width='100px',height='270px'),
        disabled=False
    )

    play = widgets.Play(value=0, min=-10, max=10, step=1, _repeat=True,
                        description="Press play", disabled=False)

    slider = widgets.IntSlider(value=0, min=-10, max=10,)

    class FnScope:
        traj=trajs[-1]

    def handle_slider_change(change):
        v.set_coordinates({0:FnScope.traj[change.new]})

    def handle_select_change(change):
        FnScope.traj = trajs[change.new]

    slider.observe(handle_slider_change, names='value')
    select.observe(handle_select_change, names='index')
    traitlets.link((play,'value'),(slider,'value'))
    vib_viewer = HBox([VBox([v,HBox([play,slider])]),select])
    display(vib_viewer)
