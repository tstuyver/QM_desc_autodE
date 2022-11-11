import autode as ade
from pebble import ProcessPool
import os
import shutil


class Molecule:
    """Class to perform autodE geometry generation for individual molecules.

    Attributes:
        id (str): molecule ID
        smiles (str): molecule SMILES
        solvent (str): Solvent in which the molecule is present.
    """

    def __init__(self, id, smiles, solvent):
        self.smiles = smiles
        self.id = id
        self.conformer_folder = str(self.id)
        self.xtb_folder = os.getcwd()
        self.xyz_file = None
        if solvent:
            self.mol = ade.Molecule(
                name=f"geometry_{self.id}",
                smiles=f"{self.smiles}",
                solvent_name=solvent,
            )
        else:
            self.mol = ade.Molecule(name=f"geometry_{self.id}", smiles=f"{self.smiles}")

        self._xyz_file = None

    @property
    def xyz_file(self):
        return self._xyz_file

    @xyz_file.setter
    def xyz_file(self, filename):
        if isinstance(filename, str):
            self._xyz_file = filename

    def generate_structure(self):
        dir_path = os.path.join(self.xtb_folder, str(self.id))

        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

        os.chdir(dir_path)
        try:
            self.mol.find_lowest_energy_conformer(lmethod=ade.methods.XTB())
            self.mol.print_xyz_file()

            self.xyz_file = f"geometry_{self.id}.xyz"
            shutil.copyfile(
                os.path.join(dir_path, self.xyz_file),
                os.path.join(self.xtb_folder, self.xyz_file),
            )
            os.chdir(self.xtb_folder)
        except Exception:
            os.chdir(self.xtb_folder)
        shutil.rmtree(dir_path)


def get_geometry(mol):
    """Get geometry of single molecule."""
    mol.generate_structure()
    if mol.xyz_file is not None:
        return mol, mol.id
    else:
        return None, mol.id


def get_opt_molecules(args, molecule_list, logger):
    """_summary_

    Args:
        args (Namespace): command line arguments
        molecule_list (List[Molecule]): list of molecules for which a geometry needs to be determined
        logger (Logging.Logger): logger

    Returns:
        opt_molecules (List[Molecule, str]): list of molecules with optimized geometry and the associated ID
    """
    opt_molecules = []

    with ProcessPool(max_workers=args.xtb_n_procs) as pool:
        future = pool.map(get_geometry, molecule_list, timeout=900)

        iterator = future.result()

        while True:
            try:
                result, id = next(iterator)
                if result != None:
                    opt_molecules.append(result)
                    logger.info(f"optimization of {id} completed")
                else:
                    logger.info(f"optimization for {id} failed")

            except StopIteration:
                break
            except TimeoutError as error:
                logger.info(f"get_geometry call took more than {error.args} seconds")
                raise
            except ValueError as error:
                logger.info(f"get_geometry call failed due to ValueError: {error.args}")
                raise
            except Exception as error:
                logger.info(
                    f"get_geometry call failed due to an unexpected error: {error.args}"
                )
                pass

        pool.close()
        pool.join()

    # opt_molecules = list(filter(None, opt_molecules))

    return opt_molecules
