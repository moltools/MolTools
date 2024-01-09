import random 
import os

import joblib
from flask import Blueprint, Response, request 
from rdkit import Chem

from .common import Status, ResponseData

from biosynfoni import overlapped_fp

# Load the predictor. If it fails, set it to None.
try:
    absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    PREDICTOR = joblib.load(os.path.join(absolute_path, "models/model.pkl"))
    PREDICTOR.set_params(n_jobs=1)
except Exception as e:
    print(e)
    PREDICTOR = None

blueprint_predict_biosynthetic_class = Blueprint("predict_biosynthetic_class", __name__)
@blueprint_predict_biosynthetic_class.route("/api/predict_biosynthetic_class", methods=["POST"])
def predict_biosynthetic_class() -> Response:
    """
    API endpoint for predicting the biosynthetic class of a molecule.

    :param str smiles: The SMILES string.
    :return: The predicted biosynthetic classes with their probabilities.
    :rtype: ResponseData
    """
    # Parse request data.
    data = request.get_json()

    # Check if the predictor is loaded.
    if PREDICTOR is None:
        msg = "No predictor loaded!"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
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
            
            # Define labels.
            pathway_labels = [
                "Alkaloid",
                "Amino acid",
                "Carbohydrate",
                "Fatty acid",
                "Isoprenoid",
                "Phenylpropanoid",
                "Polyketide"
            ]

            # Get the overlapped fingerprint.
            fp = overlapped_fp(mol)
            preds = PREDICTOR.predict_proba([fp])
            preds = [pred[0][1] for pred in preds]
            labeled_preds = list(zip(pathway_labels, preds))
            labeled_preds = {label: pred for label, pred in labeled_preds}

            # Return the response.
            payload = {"predictions": labeled_preds}

            msg = "Successfully predicted biosynthetic class!"
            return ResponseData(Status.Success, payload, msg).to_dict()
        
        except Exception as e:
            print(e)
            msg = "Could not parse SMILES string!"
            return ResponseData(Status.Failure, message=msg).to_dict()
