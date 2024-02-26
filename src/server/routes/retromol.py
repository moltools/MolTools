from flask import Blueprint, Response, request 
import os 
import json 
import typing as ty
# import re
# import joblib

# import umap 
# import rdkit
# from rdkit import Chem
from rdkit.Chem import AllChem 
# from sklearn.neighbors import KDTree

import numpy as np

from retromol.chem import Molecule, MolecularPattern, ReactionRule
# from retromol_sequencing.fingerprint import get_amino_acid_fingerprint, amino_acid_class_to_label
from retromol.parsing import parse_reaction_rules, parse_molecular_patterns, parse_mol
# from retromol_sequencing.primary_sequence import resolve_biosynthetic_sequence
from retromol_sequencing.fingerprint import get_biosynthetic_fingerprint

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

# try:
#     absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#     PREDICTOR = joblib.load(os.path.join(absolute_path, "models/aa_classifier.joblib"))
#     PREDICTOR.set_params(n_jobs=1)
# except Exception as e:
#     print(e)
#     PREDICTOR = None

# try:
#     TREE = KDTree(np.load(os.path.join(absolute_path, "models/embedding.npy")))
#     EMBEDDING_LABELS = np.loadtxt(os.path.join(absolute_path, "models/embedding_identifiers.txt"), dtype=str)
#     EMBEDDER = joblib.load(os.path.join(absolute_path, "models/embedder.joblib"))
# except Exception as e:
#     print(e)
#     TREE = None
#     EMBEDDING_LABELS = None
#     EMBEDDER = None

# def is_polyketide(name: str, motif: Chem.Mol) -> bool:
#     """
#     Check if motif is a polyketide.

#     :param str name: Name of motif.
#     :param Chem.Mol motif: Motif.
#     :returns: True if motif is a polyketide.
#     :rtype: bool
#     """
#     if re.match(r"^[A-D]\d+$", name):
#         return True
    
#     return False

# def is_amino_acid(name: str, motif: Chem.Mol) -> bool:
#     """
#     Check if motif is an amino acid.

#     :param str name: Name of motif.
#     :param Chem.Mol motif: Motif.
#     :returns: True if motif is an amino acid.
#     :rtype: bool
#     """
#     smarts_strings = [
#         r"[NH1,NH2][C][C](=[O])[O]", # Alpha amino acid
#         r"[NH1,NH2][C][C][C](=[O])[O]" # Beta amino acid
#     ]

#     if any([motif.HasSubstructMatch(Chem.MolFromSmarts(x)) for x in smarts_strings]):
#         return True
    
#     return False

# def is_primary_motif(motif: Chem.Mol) -> ty.Optional[str]:
#     """
#     Check if motif is a primary motif.

#     :param Chem.Mol motif: Motif.
#     :returns: Primary motif identity or None.
#     :rtype: ty.Optional[str]
#     """
#     for monomer in MONOMERS:
#         if motif.HasSubstructMatch(monomer.compiled):
#             name = monomer.name

#             if is_polyketide(name, motif):
#                 return name 
#             elif is_amino_acid(name, motif):
#                 return name
        
#     return False

# def classify(name: str, mol: Chem.Mol) -> str:
#     """
#     Classify monomer.
    
#     :param str name: Name of monomer.
#     :param Chem.Mol mol: Molecule representing monomer.
#     :returns: Classification.
#     :rtype: str
#     """
#     if is_polyketide(name, mol):
#         return name[0] # Any of A, B, C, D.
#     else:
#         fp = get_amino_acid_fingerprint(mol)
#         pred = PREDICTOR.predict([fp])[0]
#         label = amino_acid_class_to_label(pred)
#         return label

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

            # embed 2d coords in mol object
            AllChem.Compute2DCoords(mol.compiled)
            conf = mol.compiled.GetConformer()

            result = parse_mol(mol, REACTIONS, MONOMERS)

            primary_seq_monomer_ids = []
            primary_seq = []

            # if result.success == True:
            #     seq = resolve_biosynthetic_sequence(result, is_primary_motif, return_monomer_id=True)

            #     for i, _ in enumerate(seq): # monomer ID 
            #         # get edges
            #         if i > 0:
            #             primary_seq_monomer_ids.append([seq[i-1][2], seq[i][2]])

            #     if PREDICTOR is not None:
            #         seq = [(t[0], classify(t[0], t[1])) for t in seq]
            #         for item in seq:
            #             primary_seq.append(item[1])
            #     else:
            #         seq = []
            # else:
            #     return ResponseData(Status.Failure, message="Failed to parse molecule!").to_dict()

            monomer_to_nodes = {}
            for sub_id, (node_id, name) in result.monomer_mapping.items():
                sub_mol = result.reaction_mapping[sub_id]
                isotope_numbers = [atom.GetIsotope() - 1 for atom in sub_mol.GetAtoms() if atom.GetIsotope() > 0]
                monomer_to_nodes[node_id] = isotope_numbers

            atom_coords = {}
            for atom in mol.compiled.GetAtoms():
                atom_coords[atom.GetIdx()] = np.array([conf.GetAtomPosition(atom.GetIdx()).x, conf.GetAtomPosition(atom.GetIdx()).y])
            
            atom_colors = {
                "C": "black",
                "O": "red",
                "N": "blue",
                "F": "green",
                "Cl": "green",
                "Br": "green",
                "I": "green",
                "P": "purple",
                "S": "yellow",
            }

            molecule_graph_data = {"nodes": [], "links": []}
            for atom in mol.compiled.GetAtoms():
                molecule_graph_data["nodes"].append({
                    "id": atom.GetIdx(),
                    "name": atom.GetSymbol(),
                    "color": atom_colors.get(atom.GetSymbol(), "grey"),
                    "label": "",
                    "x": conf.GetAtomPosition(atom.GetIdx()).x,
                    "y": conf.GetAtomPosition(atom.GetIdx()).y
                })
            for bond in mol.compiled.GetBonds():
                molecule_graph_data["links"].append({
                    "bondtype": bond.GetBondTypeAsDouble(),
                    "source": bond.GetBeginAtom().GetIdx(),
                    "target": bond.GetEndAtom().GetIdx()
                })

            mapping = {}
            for _, (node_id, name) in result.monomer_mapping.items():
                mapping[node_id] = name
            monomer_graph_data = {"nodes": [], "links": []}
            for node_id in result.monomer_graph.nodes: 
                name = mapping.get(node_id, None)
                if name is not None: color = "orange"
                else: color = "grey"
                name = name if name is not None else "None"

                # calc centroid 
                atom_ids = monomer_to_nodes.get(node_id, [])
                if len(atom_ids) > 0:
                    centroid = np.mean([atom_coords[atom_id] for atom_id in atom_ids], axis=0)
                else:
                    centroid = [0, 0]

                monomer_graph_data["nodes"].append({
                    "id": node_id,
                    "name": name,
                    "color": color,
                    "atom_ids": monomer_to_nodes.get(node_id, []),
                    "x": centroid[0],
                    "y": centroid[1]
                })
            for edge in result.monomer_graph.edges:
                monomer_graph_data["links"].append({
                    "source": edge[0],
                    "target": edge[1]
                })

            payload = {
                "molecule_graph_data": molecule_graph_data,
                "monomer_graph_data": monomer_graph_data,
                "primary_seq": primary_seq,
                "primary_seq_monomer_ids": primary_seq_monomer_ids
            }

            msg = "Successfully parsed molecule!"
            return ResponseData(Status.Success, payload=payload, message=msg).to_dict()

        except Exception as e:
            print(e)
            msg = "Failed to parse molecule!"
            return ResponseData(Status.Failure, message=msg).to_dict()

blueprint_embed_retromol = Blueprint("embed_retromol", __name__)
@blueprint_embed_retromol.route("/api/embed_retromol", methods=["POST"])
def parse_retromol() -> Response:
    # Parse request data.
    data = request.get_json()

    # if seq := data.get("primary_seq", None):
        
        # fp = get_biosynthetic_fingerprint(seq)
        
        # if (EMBEDDER is not None) and (TREE is not None) and (EMBEDDING_LABELS is not None):
        #     query = EMBEDDER.transform([fp])
        #     dists, ind = TREE.query(query, k=10)
            
        #     matches = []    
        #     for i, index in enumerate(ind[0]):
        #         label = EMBEDDING_LABELS[index]
        #         dist = dists[0][i]

        #         if label.startswith("NPA"):
        #             url = f"https://www.npatlas.org/explore/compounds/{label}"
        #         else:
        #             url = None 

        #         matches.append({
        #             "index": i,
        #             "label": label.replace("_", " "),
        #             "distance": round(dist, 3),
        #             "link": url
        #         })

        #     payload = {"matches": matches}
        #     return ResponseData(Status.Success, payload=payload, message="Retrieved matches!").to_dict()
        # else:
    #     return ResponseData(Status.Success, message="Could not load embedding!").to_dict()
    
    # else:
    #     msg = "No primary sequence provided!"
    #     return ResponseData(Status.Failure, message=msg).to_dict()

    return ResponseData(Status.Failure, message="Matching not inplemented!").to_dict()