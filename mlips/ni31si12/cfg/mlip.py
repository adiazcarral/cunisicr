import os
from glob import glob
import sys

import ase
from ase import Atom, Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.io import read as ase_read

#from mlippy.mtp import Relax



def convert_to_db(filename, format='vasp', db_path="cfgs.db"):
    """
    If the filenames follow a pattern, eg., POSCAR1 POSCAR2 POSCAR3, then you will 
    write the filename as `POSCAR`

    Args: 
        filename (str): Filename (or pattern) of the structures that you want to convert
    """
    print(filename)
    input_files = glob(filename + '*')
    print(input_files)
    
    print(sorted(input_files))
    atoms = []


    for _struct in range(len(input_files)):
        struct = _struct + 1
        print(struct)
        try:
            filepath = os.path.join('/work/angel/cunisicr_mlips/phases/ni31si12/cfg/dir-OUTCAR' + str(struct), 'OUTCAR')
            print(filepath)
            _atom = ase_read(filepath, format=format)
            atoms.append(_atom)
            print(_atom)
        except:
            filepath = os.path.join('/work/angel/cunisicr_mlips/phases/ni31si12/cfg/POSCAR-' + str(struct) , 'POSCAR')
            print(filepath)
            _atom = ase_read(filepath, format='vasp')
            atoms.append(_atom)
            print(_atom)

        # Write structures to a database for recordkeeping
    with connect(db_path) as db:
        for struct in atoms:
            db.write(struct)


def cfg_to_db(cfg_path, species_map, system_yml=None):
    """
    Reads in an mlip style .cfg file, parses each structure into an ASE atoms object.
    It currently reads the file in line by line which is very memory efficient.
    However, all the structures are 
    """
#    db = connect('relaxed.db')

#    print(species_map)

    mlip_data = []
    extra_data = []
    
    with open(cfg_path, 'r') as infile:
        # Chunk parameters for pushing data to the database object
        counter = 0
        chunksize = 12000
        crystal_set = []
        

        features = {}

        for line in infile:
            if 'BEGIN_CFG' in line:
                continue
                    
            if 'Size' in line:
                cell_size = int(next(infile))
                continue

            if 'supercell' in line.lower():
                cell_vecs = []
                cell_vecs.append([float(i) for i in next(infile).split()])
                cell_vecs.append([float(i) for i in next(infile).split()])
                cell_vecs.append([float(i) for i in next(infile).split()])
                continue

            if 'AtomData' in line:
                atoms = []
                forces = []
                for i in range(cell_size):
                    adata = [float(f) for f in next(infile).split()]
                    num = int(adata[0])
                    mlip_atom_type = int(adata[1])
                    species = species_map[mlip_atom_type]
                    position = adata[2:5]
                    forces.append(adata[5:8])
#                    magmom = system_yml['magmom'][mlip_atom_type]
                    atoms.append(Atom(species, position, index=num))#, magmom=magmom))
                    
                    continue

            if 'Energy' in line:
                energy = float(next(infile))
                continue

            if 'PlusStress' in line:
                plusstress = [float(i) for i in next(infile).split()]
                continue

            if 'conf_id' in line:
                conf_id = int(line.split()[2])
                continue

            if 'Feature' in line and 'relaxation' not in line:
                feat_type = line.split()[1]
                feat_val = line.split()[2]
                #             print(feat_type, feat_val)
                features[feat_type] = feat_val
                continue

            if 'relaxation' in line:
                relax = line.split()[2].split('_')
                relax_state = relax[1]
                features[relax[0]] = relax_state
                continue

            if 'END_CFG' in line:
                """ Create the Atoms object from the array `atoms` created earlier.
                All MLIP tags are stored in `info`

                """

                # INITIALIZE CRYSTAL PROPERTIES
                crystal = Atoms(atoms, pbc=True)#, calculator=calc)
                crystal.set_cell(cell_vecs)
                
                crystal.set_calculator(
                    calc=SinglePointCalculator(
                        atoms=crystal,
                        energy=energy,
                        forces=forces,
                        stress=plusstress
                    )
                )
                


                crystal_set.append(crystal)
#                print(len(crystal_set), chunksize)
#                print(crystal_set)
                counter += 1
#                print("Now reading:", counter)


                if len(crystal_set) == chunksize:
                    #print(len(crystal_set), chunksize)
                    counter = 0
                    
                    ### Print the remainder structures to the database
                    with connect('relaxed.db') as db: 
                        for struct in crystal_set:
                            db.write(struct)
                    crystal_set = []   
        ### Print the remainder structures to the database
        with connect('relaxed.db') as db: 
            for struct in crystal_set:
                db.write(struct)




def write_to_cfg(struct, mapping):
    map_backwards = {mapping[key]:key for key in mapping}
#    print(map_backwards)
    cfg_lines = []
#    print(struct)
    cfg_lines.append("BEGIN_CFG")
    
    # Set the size of the cell
    cfg_lines.append(" Size")
    cfg_lines.append('   ' + str(struct.natoms))
    
    # Set the supercell lattice vectors
    cfg_lines.append(" Supercell")
    for vec in struct.cell:
        cfg_lines.append('   ' + '\t'.join([str(val) for val in vec]))

    # Set the AtomData values

    atom_num = 1
    try:
        struct.forces
        cfg_lines.append(" AtomData:  id type   cartes_x  cartes_y  cartes_z  fx  fy  fz")
    except:
        cfg_lines.append(" AtomData:  id type   cartes_x  cartes_y  cartes_z")

    for atom in range(len(struct.numbers)): # loop through the 
        atom_num = str(atom + 1) # correcting for the base-0 python index
#        print(struct.numbers, mapping)
        species = '  ' + str(map_backwards[struct.numbers[atom]]) # have to map the value to the key
        position = '  '.join([str(val) for val in struct.positions[atom]])
        try:        
            
            forces = '  '.join([str(val) for val in struct.forces[atom]])
            _atom_line = '\t' + atom_num + '\t' + species + '\t' + position + '\t' + forces
        except:
            _atom_line = '\t' + atom_num + '\t' + species + '\t' + position
            pass


        atom_line = ''.join([str(val) for val in _atom_line])
        cfg_lines.append(atom_line)
    
    # Set energy
    try:
        energy = struct.energy
        cfg_lines.append(" Energy")    
        cfg_lines.append('\t' + str(energy))

        cfg_lines.append(" PlusStress: xx   yy   zz   yz   xz   xy")
        cfg_lines.append('\t' + '\t'.join([str(f) for f in struct.stress]))

        
    except:
        pass
    
    cfg_lines.append(" Feature    EFS_by      VASP")
    cfg_lines.append("END_CFG\n")
    
 #   print(cfg_lines)
    return cfg_lines


    
def set_mapping(yaml_path):
    import yaml
    params = yaml.safe_load(yaml_path)
    with open(yaml_path, 'r') as y:
        yml = yaml.safe_load(y)
    species = yml['species']

    # List comprehension for writing a mapping between chemical numbers to MLIP species ID
    mapping = {s:ase.symbols.symbols2numbers(species[s])[0] for s in species}
    
    return mapping


def db_to_cfg(cfg_out_path, db_path='cfgs.db', yaml_path='system.yaml'):
    db = connect(db_path)
    

    mapping = set_mapping(yaml_path)

    with open(cfg_out_path, 'w') as f:
        for struct in db.select():
            lines = write_to_cfg(struct, mapping)

            for line in lines:
                print(line, file=f)





#if __name__ == "__main__":
#    yaml_path = sys.argv[1]
#    mapping = set_mapping(yaml_path)
#    db_to_cfg("_out.cfg", db_path="./cfgs.db", yaml_path=yaml_path)



def extract_training_error(filepath):


    with open(filepath, 'r') as infile:
        error_report = {}
        error_find = False 

        for line in infile:
            if '___Errors report____' in line: 
                error_find = True
                
            if 'Energy:' in line:
                n_configs = int(next(infile).split()[3])
                max_error = float(next(infile).split()[4])
                average_error = float(next(infile).split()[4])
                energy_rms_error = float(next(infile).split()[4])
                
                error_report['n_configs'] = n_configs
                error_report['Energy'] = {'max':max_error,
                                          'avg':average_error,
                                          'rms':energy_rms_error
                }

            else:
                pass
                

        
            if 'Energy per atom:' in line:
                next(infile)
                max_error = float(next(infile).split()[4])
                average_error = float(next(infile).split()[4])
                energy_rms_error = float(next(infile).split()[4])
                
                error_report['Energy_per_atom'] = {'max':max_error,
                                                   'avg':average_error,
                                                   'rms':energy_rms_error
                }

            else:
                pass

            if 'Forces:' in line:
                n_atoms = int(next(infile).split()[3])
                max_error = float(next(infile).split()[4])
                average_error = float(next(infile).split()[4])
                energy_rms_error = float(next(infile).split()[4])
                max_stress_diff = float(next(infile).split()[4])
                rms_stress_diff = float(next(infile).split()[4])

                error_report['Energy_per_atom'] = {'max':max_error,
                                                   'avg':average_error,
                                                   'rms':energy_rms_error,
                                                   'max_stress_diff':max_stress_diff,
                                                   'rms_stress_diff':rms_stress_diff
                }

            else:
                pass


            if 'Stresses:' in line:
                next(infile)
                max_error = float(next(infile).split()[4])
                average_error = float(next(infile).split()[4])
                energy_rms_error = float(next(infile).split()[4])
                max_stress_diff = float(next(infile).split()[4])
                rms_stress_diff = float(next(infile).split()[4])

                error_report['Energy_per_atom'] = {'max':max_error,
                                                   'avg':average_error,
                                                   'rms':energy_rms_error,
                                                   'max_stress_diff':max_stress_diff,
                                                   'rms_stress_diff':rms_stress_diff
                }

            else:
                pass


            if 'Virial stresses:' in line:
                max_error = float(next(infile).split()[4])
                average_error = float(next(infile).split()[4])
                energy_rms_error = float(next(infile).split()[4])
                max_stress_diff = float(next(infile).split()[4])
                rms_stress_diff = float(next(infile).split()[4])

                error_report['Virial_stresses'] = {'max':max_error,
                                                   'avg':average_error,
                                                   'rms':energy_rms_error,
                                                   'max_stress_diff':max_stress_diff,
                                                   'rms_stress_diff':max_stress_diff
                }

            else:
                pass

    return error_report
