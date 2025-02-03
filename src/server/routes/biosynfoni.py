import os
import subprocess

import joblib
from flask import Blueprint, Response, request 
from rdkit import Chem

from .common import fail, success

import biosynfoni
from biosynfoni import overlapped_fp

# Load the predictor. If it fails, set it to None.
PREDICTOR_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models/bsf_model.joblib")

blueprint_predict_biosynthetic_class = Blueprint("predict_biosynthetic_class", __name__)
@blueprint_predict_biosynthetic_class.route("/api/predict_biosynthetic_class", methods=["POST"])
def predict_biosynthetic_class() -> Response:
    """
    API endpoint for predicting the biosynthetic class of a molecule.
    """
    # Parse request data.
    data = request.get_json()

    # load model
    try:
        PREDICTOR = joblib.load(PREDICTOR_PATH)
        PREDICTOR.set_params(n_jobs=1)
    except Exception as e:
        print(e)
        PREDICTOR = None

    # Check if the predictor is loaded.
    if PREDICTOR is None:
        msg = "No predictor loaded!"
        return fail(msg)
    
    # Retrieve the SMILES string.
    smiles = data.get("smiles", None)
    
    if smiles is None:
        msg = "No SMILES string provided!"
        del PREDICTOR
        fail(msg)
    
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
            del PREDICTOR
            return success(msg, payload)
        
        except Exception as e:
            print(e)
            msg = "Could not parse SMILES string!"
            del PREDICTOR
            return fail(msg)
        
# end point for getting biosynfoni version of installed package
blueprint_biosynfoni_version = Blueprint("biosynfoni_version", __name__)
@blueprint_biosynfoni_version.route("/api/biosynfoni_version", methods=["GET"])
def biosynfoni_version() -> Response:
    """
    API endpoint for getting the version of the installed biosynfoni package.
    """
    version = "0.0.0"
    try:
        output = subprocess.run(
            ["pip", "show", "biosynfoni"], capture_output=True, text=True
        )
        for line in output.stdout.split("\n"):
            if line.startswith("Version:"):
                version = line.split(":", 1)[1].strip()
    except Exception as e:
        pass

    # version = biosynfoni.__version__
    payload = {"version": version}

    msg = "Successfully retrieved biosynfoni version!"
    return success(msg, payload)