from flask import Blueprint, Response, request 
import os 
import json 
import typing as ty

from tqdm import tqdm

from rdkit import Chem

import neo4j

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.parsing import parse_reaction_rules, parse_molecular_patterns, parse_mol
from retromol_sequencing.sequencing import parse_modular_natural_product
from retromol_sequencing.alignment import ModuleSequence, parse_primary_sequence, MultipleSequenceAlignment, PolyketideMotif, PeptideMotif, Gap

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

blueprint_bioactivity_labels = Blueprint("bioactivity_labels", __name__)
@blueprint_bioactivity_labels.route("/api/bioactivity_labels", methods=["GET"])
def bioactivity_labels() -> Response:
    biaoctivity_labels = []

    driver = neo4j.GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))

    with driver.session() as session:
        query = """
        MATCH (b:Bioactivity)
        RETURN b
        """
        result = session.run(query)

        bioactivity_labels = []
        for record in result:
            bioactivity_labels.append(record["b"]["name"])

    # sort bioactivity labels
    bioactivity_labels = list(set(bioactivity_labels))
    bioactivity_labels.sort()

    payload = {"bioactivityLabels": bioactivity_labels}
    msg = "Successfully retrieved bioactivity labels!"
    return ResponseData(Status.Success, payload=payload, message=msg).to_dict()

blueprint_parse_retromol = Blueprint("parse_retromol", __name__)
@blueprint_parse_retromol.route("/api/parse_retromol", methods=["POST"])
def parse_retromol() -> Response:
    data = request.get_json()
    data = data["data"]

    smiles = data.get("smiles", None)
    if smiles is None:
        msg = "No SMILES string provided!"

        return ResponseData(Status.Failure, message=msg).to_dict()
    
    else:
        mol = Molecule("input", smiles)
        result = parse_mol(mol, REACTIONS, MONOMERS)

        input_smiles = Chem.MolToSmiles(result.mol)

        # Get all monomer amns
        monomer_ids = []
        for k, props in result.monomer_graph.items():
            if props["identity"] is not None:
                monomer_ids.append(props["reaction_tree_id"])
        subs = [Chem.MolFromSmiles(result.reaction_tree[x]["smiles"]) for x in monomer_ids]

        monomer_amns = set()
        for sub in subs:
            for atom in sub.GetAtoms():
                amn = atom.GetAtomMapNum()
                if amn > 0: monomer_amns.add(amn)
        
        # filter reaction tree for nodes that contain all amns 
        mols = []
        for k, props in result.reaction_tree.items():
            mol = Chem.MolFromSmiles(props["smiles"])
            amns = set()
            for atom in mol.GetAtoms():
                amn = atom.GetAtomMapNum()
                if amn > 0: amns.add(amn)
            if monomer_amns.issubset(amns):
                mols.append(mol)

        # get num of cycles in mol, keep ones with smallest num of cycles
        mols_with_num_cycles = []
        for mol in mols:
            ssr = Chem.GetSymmSSSR(mol)
            num_cycles = len(ssr)
            mols_with_num_cycles.append((mol, num_cycles))
        min_num_cycles = min([x[1] for x in mols_with_num_cycles])
        mols = [x[0] for x in mols_with_num_cycles if x[1] == min_num_cycles]
        linearized_smiles = Chem.MolToSmiles(mols[0])
        


        sequences = parse_modular_natural_product(result)

        # 
        input_mol = Chem.MolFromSmiles(input_smiles)
        for atom in input_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        input_smiles = Chem.MolToSmiles(input_mol)

        # 
        linearized_mol = Chem.MolFromSmiles(linearized_smiles)
        for atom in linearized_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        linearized_smiles = Chem.MolToSmiles(linearized_mol)

        payload = {
            "sequences": sequences,
            "input_smiles": input_smiles,
            "linearized_smiles": linearized_smiles
        }
        message = "Successfully parsed molecule!"

        return ResponseData(Status.Success, payload=payload, message=message).to_dict()
    
blueprint_embed_retromol = Blueprint("embed_retromol", __name__)
@blueprint_embed_retromol.route("/api/embed_retromol", methods=["POST"])
def embed_retromol() -> Response:
    data = request.get_json()

    return ResponseData(Status.Failure, message="Matching not inplemented!").to_dict()

def retrieve_primary_sequence(session: neo4j.Session, identifier: str) -> ty.List[ty.Dict[str, ty.Any]]:
    query = """
    MATCH (:PrimarySequence {identifier: $identifier})-[:START]->(startUnit)
    MATCH s = (startUnit)-[:NEXT*]->(endUnit)
    WITH nodes(s) AS nodes, length(s) AS path_length
    ORDER BY path_length DESC
    LIMIT 1
    RETURN nodes
    """
    result = session.run(query, identifier=identifier)
    record = result.single()

    seq = []
    if record is not None:
        if nodes := record["nodes"]:
            for node in nodes:
                props = dict(node.items())
                seq.append(props)

    return seq 

blueprint_find_matches = Blueprint("find_matches", __name__)
@blueprint_find_matches.route("/api/find_matches", methods=["POST"])
def find_matches() -> Response:
    data = request.get_json()
    data = data["data"]

    payload = {"matches": []}

    try:
        # driver = neo4j.GraphDatabase.driver("bolt://database:7687")
        driver = neo4j.GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))

        query = [] 
        for item in data["matchItems"]:
            props = item["properties"]
            props["identifier"] = item["identifier"]
            query.append(props)
        seq_a = ModuleSequence("Query", parse_primary_sequence(query))

        gap_cost = 2
        end_gap_cost = 1
        
        with driver.session() as session:
            query = """
            MATCH (b:PrimarySequence)
            RETURN b
            """
            result = session.run(query, fetch_size=1)

            top_10 = []
            for record in tqdm(result):

                id_b = record["b"]["identifier"]
                seq_b = retrieve_primary_sequence(session, id_b)

                if len(seq_b) != 0:
                    seq_b = ModuleSequence(id_b, parse_primary_sequence(seq_b))
                    score = seq_a.optimal_alignment(seq_b, gap_cost, end_gap_cost).score

                    # Keep 10 best scores.
                    if len(top_10) < 50:
                        top_10.append((id_b, score, seq_b))
                        top_10.sort(key=lambda x: x[1], reverse=True)
                    else:
                        if score > top_10[-1][1]:
                            top_10.pop()
                            top_10.append((id_b, score, seq_b))
                            top_10.sort(key=lambda x: x[1], reverse=True)

                else:
                    continue

            # get top 10 sequences and do msa 
            msa = MultipleSequenceAlignment([seq_a] + [x[2] for x in top_10], gap_cost, end_gap_cost).get_alignment()

            matches = []
            for seq in msa:

                seq_repr = []
                for x, tag in seq._seq: 
                    if isinstance(x, PolyketideMotif):
                        x_repr = x.type.name + (str(x.decoration_type) if x.decoration_type is not None else "")
                    elif isinstance(x, PeptideMotif):
                        x_repr = "AA" + str(x.type.value)
                    elif isinstance(x, Gap):
                        x_repr = "GAP"
                    else:
                        x_repr = "???"
                    seq_repr.append(x_repr)

                if seq.name.startswith("NPA"):
                    url = "https://www.npatlas.org/explore/compounds/" + seq.name
                else:
                    url = None 

                bioactivities = []
                query = """
                MATCH (c:Compound {identifier: $identifier})-[:HAS_BIOACTIVITY]->(b:Bioactivity)
                RETURN b
                """
                result = session.run(query, identifier=seq.name)
                for record in result:
                    bioactivities.append(record["b"]["name"])

                matches.append({
                    "identifier": seq.name,
                    "bioactivity": list(set(bioactivities)),
                    "sequence": seq_repr,
                    "url": url
                })

            payload["matches"] = matches

        driver.close()
        
        return ResponseData(Status.Success, payload=payload, message="Matching completed!").to_dict()
    
    except Exception as e:
        return ResponseData(Status.Failure, message=f"{e}").to_dict()
    
blueprint_query_database = Blueprint("query_database", __name__)
@blueprint_query_database.route("/api/query_database", methods=["POST"])
def query_database() -> Response:
    data = request.get_json()

    payload = {"matches": []}

    try:
        # driver = neo4j.GraphDatabase.driver("bolt://database:7687")
        driver = neo4j.GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))

        gap_cost = 2
        end_gap_cost = 1
        
        with driver.session() as session:
            query = """
            MATCH (b:PrimarySequence)-[:START]->(u1)
            MATCH path = (u1)-[:NEXT*]->(u2)-[:NEXT]->(u3)-[:NEXT]->(u4)-[:NEXT]->(u5)
            WHERE
            (u2.identifier = 'polyketide' AND u2.accessory_domains = ['KR', 'DH'] AND u2.decoration_type = '4') AND
            (u3.identifier = 'polyketide' AND u3.accessory_domains = ['KR', 'DH'] AND u3.decoration_type = '1') AND
            (u4.identifier = 'polyketide' AND u4.accessory_domains = ['KR']) AND
            (u5.identifier = 'polyketide' AND u5.accessory_domains = ['KR', 'DH'] AND u5.decoration_type = '1')
            AND NOT (u5)-[:NEXT]->()
            RETURN DISTINCT b
            """
            result = session.run(query, fetch_size=1)

            first_10 = []
            for record in tqdm(result):

                id_b = record["b"]["identifier"]
                seq_b = retrieve_primary_sequence(session, id_b)

                if len(seq_b) != 0:
                    seq_b = ModuleSequence(id_b, parse_primary_sequence(seq_b))
                    first_10.append((id_b, seq_b))

                else:
                    continue

            print(first_10)

            # get top 10 sequences and do msa 
            msa = MultipleSequenceAlignment([x[1] for x in first_10], gap_cost, end_gap_cost).get_alignment()

            matches = []
            for seq in msa:

                seq_repr = []
                for x, tag in seq._seq: 
                    if isinstance(x, PolyketideMotif):
                        x_repr = x.type.name + (str(x.decoration_type) if x.decoration_type is not None else "")
                    elif isinstance(x, PeptideMotif):
                        x_repr = "AA" + str(x.type.value)
                    elif isinstance(x, Gap):
                        x_repr = "GAP"
                    else:
                        x_repr = "???"
                    seq_repr.append(x_repr)

                if seq.name.startswith("NPA"):
                    url = "https://www.npatlas.org/explore/compounds/" + seq.name
                else:
                    url = None 

                bioactivities = []
                query = """
                MATCH (c:Compound {identifier: $identifier})-[:HAS_BIOACTIVITY]->(b:Bioactivity)
                RETURN b
                """
                result = session.run(query, identifier=seq.name)
                for record in result:
                    bioactivities.append(record["b"]["name"])

                matches.append({
                    "identifier": seq.name,
                    "bioactivity": list(set(bioactivities)),
                    "sequence": seq_repr,
                    "url": url
                })

            payload["matches"] = matches

        driver.close()
        
        return ResponseData(Status.Success, payload=payload, message="Matching completed!").to_dict()
    
    except Exception as e:
        return ResponseData(Status.Failure, message=f"{e}").to_dict()