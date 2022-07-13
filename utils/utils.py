from rdkit import Chem

def smiles_concat(first_smiles: list, second_smiles: list):
    """
    Given two list of smiles, returns the their list with the concatenated version. 
    Examples: A='CO' and B='CN' returns CO.CN
    """
    if bool(first_smiles) is False:
        raise ValueError("You should know not to pass an empty list, dumbass.")
    if bool(second_smiles) is False:
        raise ValueError("You should know not to pass an empty list, dumbass.")

    full_smiles = []
    for (smile_1, smile_2) in zip(first_smiles, second_smiles):
        instance = type(Chem.rdchem.Mol())
        if (isinstance(Chem.MolFromSmiles(smile_1), instance), isinstance(Chem.MolFromSmiles(smile_2), instance)) == (True, True):
            full_smiles.append(smile_1 + "." + smile_2)
        else:
           raise ValueError(f"Invalid smiles on on {smiles_1} or {smiles_2}.")
        
    return full_smiles
