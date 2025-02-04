import os
import subprocess

import joblib
import shap
import numpy as np
from flask import Blueprint, Response, request 
from rdkit import Chem

from .common import fail, success

from biosynfoni import overlapped_fp
from biosynfoni.subkeys.versionfonis import fpVersions
from biosynfoni.subkeys.biosmartfonis import substructureSmarts

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

            # get feature importances
            subkeys_names = fpVersions.get("full_1103", None)
            subkeys = [substructureSmarts.get(subkey, None) for subkey in subkeys_names]

            if subkeys is None or any([subkey is None for subkey in subkeys]):
                msg = "Unable to retrieve feature importances"
                del PREDICTOR
                return fail(msg)
            
            # add tags as istopes to input molecule
            for atom in mol.GetAtoms():
                index = atom.GetIdx()
                atom.SetIsotope(index + 1)
            
            # match most important subkeys to compound
            highlights = {}
            for subkey in subkeys:
                name = subkey["name"]
                smarts = subkey["smarts"]
                pattern = Chem.MolFromSmarts(smarts)
                matching_atoms = set()
                matches = mol.GetSubstructMatches(pattern)
                for match in matches:
                    for atom_idx in match:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        matching_atoms.add(atom.GetIsotope())
                highlights[name] = list(matching_atoms)

            # get feature insights
            explainer = shap.TreeExplainer(PREDICTOR)
            shap_values = explainer.shap_values(np.array(fp))

            # print num classes predictor
            shap_instance = shap_values.T[1]

            class_names = [(pathway_labels[i], 1 + (i * 2)) for i in range(len(pathway_labels))]
            feature_names = [subkey["name"] for subkey in subkeys]

            # rank which features in the fingerprint are most important for the prediction of every class
            feature_importance = {}
            for i, (class_name, idx) in enumerate(class_names):
                shap_values_class = shap_values.T[idx]
                
                # get all features that contributed positively and all negatively to the prediction
                all_pos = np.where(shap_values_class >= 0)[0]
                all_neg = np.where(shap_values_class < 0)[0]

                # top_pos = np.argsort(shap_values_class)[-3:][::-1]
                # top_neg = np.argsort(shap_values_class)[:3]

                feature_importance[class_name] = {
                    "top_pos": [feature_names[j] for j in all_pos],
                    "top_neg": [feature_names[j] for j in all_neg]
                }
            
            fp_keys = []
            for i, key in enumerate(fp):
                fp_keys.append((subkeys[i]["name"], key))

            # Return the response.
            payload = {
                "tagged_smiles": Chem.MolToSmiles(mol),
                "predictions": labeled_preds,
                "feature_importance": feature_importance,
                "highlights": highlights,
                "fingerprint": fp_keys
            }

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