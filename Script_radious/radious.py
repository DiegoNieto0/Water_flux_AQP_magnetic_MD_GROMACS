from MDAnalysis.analysis import hole2 as hole
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

path = os.getcwd()
pathimages = '/PATH/TO/OUTPUT'
patron = 'trajectories_names'
pdbfile = f'{path}/{patron}.pdb'
xtcfile = f'{path}/{patron}.xtc'
acceptor = True
selection='(id 1:4066)'
treatment='0a'
mon = 1

# Cargar el archivo de simulación
u = mda.Universe(f'{pdbfile}', f'{xtcfile}')

# Seleccionar los monómeros
monomer1 = u.select_atoms(f'protein and {selection}')
center1 = [69.49, 56.31, 42.14]
# Crear objetos de análisis de poro
ha = hole.HoleAnalysis(u, select=f'protein and {selection}',
                       cpoint=center1,
                       executable='~/hole2/bin/hole',
                       cvect=[0, 0, 1]
                       )
# Iniciar los bucles de tiempo y cálculo de poro
ha.run(random_seed=3145, step=1, start=1, stop=2500)
#ha.create_vmd_surface(filename='surf_0Ta.vmd')
# Bucle para eliminar archivos holeXXX.out y holeXXX.sph
n = 2500

for i in range(1, n + 1):
    # Crear el nombre de archivo con el formato holeXXX.out
    out_file = f'hole{i:03d}.out'
    
    # Eliminar el archivo si existe
    if os.path.exists(out_file):
        os.remove(out_file)

    # Crear el nombre de archivo con el formato holeXXX.sph
    sph_file = f'hole{i:03d}.sph'
    
    # Eliminar el archivo si existe
    if os.path.exists(sph_file):
        os.remove(sph_file)
means, edges = ha.histogram_radii(bins=100, range=None,
                                  aggregator=np.mean)
                                  
# Obtener datos de radio promedio y desviación estándar
std_devs, _ = ha.histogram_radii(bins=100, range=None, aggregator=np.std)

midpoints = 0.5*(edges[1:]+edges[:-1])

# Crear un DataFrame
df = pd.DataFrame({'Radio_Promedio': means,
                   'Desviacion_Estandar': std_devs,
                   'Valor_Z': midpoints})
                  
# Guardar DataFrame en un archivo CSV
#df.to_csv(f'{pathimages}/radius_{treatment}_{mon}.csv', index=False)

# Plotear los resultados
ha.plot_mean_profile(bins=100, n_std=1, color='red', fill_alpha=0.2, legend=True)
#plt.yticks(0,15)
plt.savefig('mat.png',dpi=800)
plt.show()
ha.plot()



