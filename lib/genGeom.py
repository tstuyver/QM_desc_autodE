import autode as ade
from pebble import ProcessPool

class Molecule:

    def __init__(self, id, smiles, solvent):
        self.smiles = smiles
        self.id = id
        if solvent:
            self.mol = ade.Molecule(name = f'geometry_{self.id}', smiles=f'{self.smiles}', solvent_name = solvent)
        else:
            self.mol = ade.Molecule(name = f'geometry_{self.id}',smiles=f'{self.smiles}')

        self._xyz_file = None

    @property
    def xyz_file(self):
        return self._xyz_file

    @xyz_file.setter
    def xyz_file(self, filename):
        if isinstance(filename,str):
            self._xyz_file = filename

    def generate_structure(self):
        self.mol.find_lowest_energy_conformer(lmethod=ade.methods.XTB())
        self.mol.print_xyz_file()

        self.xyz_file = f'geometry_{self.id}.xyz'


def get_geometry(mol):
    mol.generate_structure()

    return mol, mol.id


def get_opt_molecules(args, molecule_list, logger):
    opt_molecules = []

    with ProcessPool(max_workers=args.xtb_n_procs) as pool:
        future = pool.map(get_geometry, molecule_list, timeout=900)

        iterator = future.result()
        
        while True:
            try: 
                result, id = next(iterator)
                opt_molecules.append(result)

                logger.info(f'optimization of {id} completed')

            except StopIteration:
                break
            except TimeoutError as error:
                logger.info(f'get_geometry call took more than {error.args} seconds')
                raise
            except ValueError as error:
                logger.info(f'get_geometry call failed due to ValueError: {error.args}')
                raise
            except Exception as error:
                logger.info(f'get_geometry call failed due to an unexpected error: {error.args}')
                pass

        pool.close()
        pool.join()

    return opt_molecules
