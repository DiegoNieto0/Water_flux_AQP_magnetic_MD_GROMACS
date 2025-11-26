import pickle
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis 

path = os.getcwd()
pathimages = 'PATH/TO/OUTPUT'
patron = 'trajectories_names'
tprfile = f'{path}/{patron}.tpr'
pdbfile = f'{path}/{patron}.pdb'
xtcfile = f'{path}/{patron}.xtc'
waterbox = ' and (prop z < 54) and (prop z > 37) and (prop x < 99) and (prop x > 50) and (prop y < 99) and (prop y > 50)' # channel delimitation
acceptor = False
selection='(id 1:16264)'
treatment='8b'
name = 'Acceptor' if acceptor else 'Donor'

def run_hydrogen_bond_analysis(tprfile, pdbfile, xtcfile, waterbox,acceptor,selection):
    #Run hydrogen bond analysis on MD trajectory.

    #Parameters:
    #- tprfile (str): Path to the .tpr file.
    #- pdbfile (str): Path to the .pdb file.
    #- xtcfile (str): Path to the .xtc file.
    #- waterbox (str): Selection criteria for water box.
    #- acceptor (bool): If Tytrue, acceptor atoms from protein
    #- selection (str): atoms selection from protein
    #
    if acceptor == True:
        u = mda.Universe(tprfile, pdbfile, xtcfile)
        hbonds = HydrogenBondAnalysis(
            universe=u,
            hydrogens_sel=f'(name HW1 or name HW2){waterbox} ',
            acceptors_sel=f'{selection} and (type O or type N)',
        )
    else:
        u = mda.Universe(tprfile, pdbfile, xtcfile)
        hbonds = HydrogenBondAnalysis(
            universe=u,
            hydrogens_sel=f'type H and {selection}',
            acceptors_sel=f'name OW{waterbox}',
        )
        
    hbonds.run(
        start=None,
        stop=None,
        step=1,
        verbose=True
    )
    tau_max = 25
    window_step = 1

# bins in z for the histogram
    bin_edges = np.linspace(38, 55, 20)
    bin_centers = bin_edges[:-1] + 0.5
    counts = np.full(bin_centers.size, fill_value=0.0)
    for frame, donor_ix, *_ in hbonds.results.hbonds:

        u.trajectory[frame.astype(int)]
        donor = u.atoms[donor_ix.astype(int)]

        zpos = donor.position[2]
        hist, *_ = np.histogram(zpos, bins=bin_edges)
        counts += hist * 2  # multiply by two as each hydrogen bond involves two water molecules


    plt.plot(bin_centers, counts, lw=2)

#plt.title(r"Number of hydrogen bonds as a funcion of height in $z$", weight="bold")
    plt.xlabel(r"$z\ \rm (\AA)$")
    plt.ylabel(r"$N_{HB}$")

    plt.show()

# Crear un DataFrame con bin_centers y counts
    df = pd.DataFrame({'bin_centers': bin_centers, 'counts': counts})

# Agregar más información al DataFrame si es necesario

# Guardar el DataFrame en un archivo CSV
    df.to_csv(f'{pathimages}/hbonds_{name}_{treatment}.csv', index=False)

# Print donor, hydrogen, acceptor and count info for these hbonds
    counts = hbonds.count_by_ids()
    lines = []
    for donor, hydrogen, acceptor, count in counts[:10]:
        d, h, a = u.atoms[donor], u.atoms[hydrogen], u.atoms[acceptor]
        lines.append(f"{d.resname}-{d.resid}-{d.name}\t{h.name}\t{a.resname}-{a.resid}-{a.name}\tcount={count}")
    for line in sorted(lines):
        print(line)
    


