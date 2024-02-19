from flask import Blueprint, Response, request 
import os 
import json 
import typing as ty

import time
import timeout_decorator

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.parsing import parse_reaction_rules, parse_molecular_patterns, parse_mol

from .common import Status, ResponseData

try:
    absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rules = json.load(open(os.path.join(absolute_path, "data/rules.json"), "r"))
    REACTIONS = parse_reaction_rules(json.dumps(rules["reactions"]))
    MONOMERS = parse_molecular_patterns(json.dumps(rules["monomers"]))
except Exception as e:
    print(e)
    REACTIONS = []
    MONOMERS = []

Record = ty.Tuple[Molecule, ty.List[ReactionRule], ty.List[MolecularPattern]]

blueprint_parse_retromol = Blueprint("parse_retromol", __name__)
@blueprint_parse_retromol.route("/api/parse_retromol", methods=["POST"])
def parse_retromol() -> Response:
    # Parse request data.
    data = request.get_json()
    
    # Retrieve the SMILES string.
    smiles = data.get("smiles", None)

    if smiles is None:
        msg = "No SMILES string provided!"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    else:
        try:
            mol = Molecule("input", smiles)
            result = parse_mol(mol, REACTIONS, MONOMERS)

            # TODO: create easy to parse graph structure for client
            print(result.monomer_mapping)
            print(result.monomer_graph)

            payload = {"monomer_graph": json.loads(result.to_json())["monomer_graph"]}
            msg = "Successfully parsed molecule!"
            return ResponseData(Status.Success, payload=payload, message=msg).to_dict()

        except Exception as e:
            print(e)
            msg = "Failed to parse molecule!"
            return ResponseData(Status.Failure, message=msg).to_dict()
