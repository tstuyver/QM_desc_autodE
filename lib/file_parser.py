from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd

def xyz2com(xyz, head, footer, comfile, charge=0, mult=1):
    title = xyz.splitlines()[1]
    coords = [x+'\n' for x in xyz.splitlines()[2:]]

    with open(comfile, 'w') as com:
        com.write(head)
        com.write('\n')
        com.write(title+'\n')
        com.write('\n')
        com.write('{} {}\n'.format(charge, mult))
        com.writelines(coords)
        com.write('\n')
        com.write(footer)
        com.write('\n\n\n')


def NBO2csv(out, acsv):
    #TODO bcsv
    with open(out, 'r') as outhandle:
        txt = outhandle.readlines()

    txt = [x.strip() for x in txt]
    df,txt = _GetNPACharge(txt)

    return df    
             

def _GetNPACharge(txt):
    columns = 'Atom No    Charge        Core      Valence    Rydberg      Total'
    start_id = txt.index(columns)+2
    end_id = start_id + txt[start_id:].index('====================================================================')

    NPACharge = txt[start_id:end_id] 

    NPACharge = [x.split() for x in NPACharge]

    df = pd.DataFrame(NPACharge, columns=columns.split())
    
    return df, txt[end_id:]
