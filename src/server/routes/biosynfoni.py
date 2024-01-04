import random 

from flask import Blueprint, Response, request 
from rdkit import Chem

from .common import Status, ResponseData

blueprint_predict_biosynthetic_class = Blueprint("predict_biosynthetic_class", __name__)
@blueprint_predict_biosynthetic_class.route("/api/predict_biosynthetic_class", methods=["POST"])
def predict_biosynthetic_class() -> Response:
    """
    API endpoint for predicting the biosynthetic class of a molecule.

    :param str smiles: The SMILES string.
    :return: The predicted biosynthetic classes with their probabilities.
    :rtype: ty.Dict[str, float]
    """
    # Parse request data.
    data = request.get_json()
    
    # Retrieve the SMILES string.
    smiles = data.get("smiles", None)
    
    if smiles is None:
        msg = "No SMILES string provided!"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    else:
        # Parse the SMILES string.
        try:
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                raise ValueError()
            
            # Predict biosynthetic classes for the molecule. 
            # For now, get random floats between 0 and 1.
            predictions = {
                "Polyketide":   random.random(),
                "NRP":          random.random(),
                "Saccharide":   random.random(),
                "Terpene":      random.random(),
                "Alkaloid":     random.random(),
                "Nucleoside":   random.random(),
                "RiPP":         random.random(),
                "Other":        random.random()
            }
            payload = {"predictions": predictions}

            msg = "Successfully predicted biosynthetic class!"
            return ResponseData(Status.Success, payload, msg).to_dict()
        
        except Exception:
            msg = "Could not parse SMILES string!"
            return ResponseData(Status.Failure, message=msg).to_dict()
