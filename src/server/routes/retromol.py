from flask import Blueprint, Response, request 
import os 
import json 
import typing as ty

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.parsing import parse_reaction_rules, parse_molecular_patterns, parse_mol
from retromol_sequencing.sequencing import parse_modular_natural_product

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
    data = request.get_json()

    smiles = data.get("smiles", None)
    if smiles is None:
        msg = "No SMILES string provided!"

        return ResponseData(Status.Failure, message=msg).to_dict()
    
    else:
        mol = Molecule("input", smiles)
        result = parse_mol(mol, REACTIONS, MONOMERS)
        sequences = parse_modular_natural_product(result)
        payload = {"sequences": sequences}
        message = "Successfully parsed molecule!"

        return ResponseData(Status.Success, payload=payload, message=message).to_dict()
    
blueprint_embed_retromol = Blueprint("embed_retromol", __name__)
@blueprint_embed_retromol.route("/api/embed_retromol", methods=["POST"])
def embed_retromol() -> Response:
    data = request.get_json()

    return ResponseData(Status.Failure, message="Matching not inplemented!").to_dict()