"""This module contains the API endpoints for RetroMol."""
import os
import json
import typing as ty

import neo4j
from flask import Blueprint, Response, request

from retromol.chem import Molecule, MolecularPattern, ReactionRule
from retromol.parsing import parse_reaction_rules, parse_molecular_patterns, parse_mol
from retromol.sequencing import parse_modular_natural_product
from retromol.alignment import (
    ModuleSequence,
    MultipleSequenceAlignment,
    PolyketideMotif,
    PeptideMotif,
    Gap,
    parse_primary_sequence,
)

from .common import Status, ResponseData

try:
    absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rules = json.load(open(os.path.join(absolute_path, "data/rules.json"), "r", encoding="utf-8"))
    REACTIONS = parse_reaction_rules(json.dumps(rules["reactions"]))
    MONOMERS = parse_molecular_patterns(json.dumps(rules["monomers"]))
except Exception:
    REACTIONS = []
    MONOMERS = []

Record = ty.Tuple[Molecule, ty.List[ReactionRule], ty.List[MolecularPattern]]

# ======================================================================================================================
#
# Retrieving data from the database.
#
# ======================================================================================================================

blueprint_bioactivity_labels = Blueprint("bioactivity_labels", __name__)
@blueprint_bioactivity_labels.route("/api/bioactivity_labels", methods=["GET"])
def bioactivity_labels() -> Response:
    """API endpoint for retrieving all bioactivity labels.
    
    :return: The bioactivity labels.
    :rtype: ResponseData
    """
    driver = neo4j.GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))

    # Retrieve all bioactivity labels.
    labels = []
    with driver.session() as session:
        query = """
        MATCH (b:Bioactivity)
        RETURN b
        """
        result = session.run(query)
        for record in result:
            labels.append(record["b"]["name"])

    labels = list(set(labels))
    labels.sort()
    payload = {"bioactivityLabels": labels}
    message = "Successfully retrieved bioactivity labels!"
    return ResponseData(Status.Success, payload=payload, message=message).to_dict()

# ======================================================================================================================
#
# Parsing a SMILES string.
#
# ======================================================================================================================

blueprint_parse_smiles = Blueprint("parse_smiles", __name__)
@blueprint_parse_smiles.route("/api/parse_smiles", methods=["POST"])
def parse_smiles() -> Response:
    """API endpoint for parsing a SMILES string.
    
    :return: Parsed biosynthetic fingerprints.
    :rtype: ResponseData
    """
    data = request.get_json()
    data = data["data"]

    smiles = data.get("smiles", None)
    if smiles is None:
        message = "No SMILES string provided!"
        return ResponseData(Status.Failure, message=message).to_dict()
    else:
        mol = Molecule("input", smiles)
        result = parse_mol(mol, REACTIONS, MONOMERS)
        sequences = parse_modular_natural_product(result)
        payload = {"sequences": sequences}
        message = "Successfully parsed molecule!"
        return ResponseData(Status.Success, payload=payload, message=message).to_dict()

# ======================================================================================================================
#
# Parsing a proto-cluster.
#
# ======================================================================================================================

blueprint_parse_proto_cluster = Blueprint("parse_proto_cluster", __name__)
@blueprint_parse_proto_cluster.route("/api/parse_proto_cluster", methods=["POST"])
def parse_proto_cluster() -> Response:
    """API endpoint for parsing a proto=cluster.
    
    :return: Parsed biosynthetic fingerprints.
    :rtype: ResponseData
    """
    data = request.get_json()
    data = data["data"]

    try:
        selected_input_type = data["selectedInputType"] # jobId or jsonSrc
        job_id = data["jobId"]
        ncbiAccession = data["ncbiAccession"]
        jsonSrc = data["jsonSrc"]
    except KeyError as err:
        message = f"Key not present in submission: {err}"
        return ResponseData(Status.Failure, message=message).to_dict()

    message = "This endpoint is not implemented yet!"
    return ResponseData(Status.Warning, message=message).to_dict()

# ======================================================================================================================
#
# Matching and querying the database with a primary sequence.
#
# ======================================================================================================================

def retrieve_primary_sequence(session: neo4j.Session, identifier: str) -> ty.List[ty.Dict[str, ty.Any]]:
    """Retrieve the primary sequence of a compound.

    :param session: The Neo4j session.
    :type session: neo4j.Session
    :param identifier: The identifier of the compound.
    :type identifier: str
    :return: The primary sequence of the compound.
    :rtype: ty.List[ty.Dict[str, ty.Any]]
    """
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

def match_exact(
    session: neo4j.Session,
    match_items: ty.List[ty.Dict[str, ty.Any]],
    top_n: int,
    match_against_molecules: bool,
    match_against_proto_clusters: bool,
    selected_bioactivity_labels: ty.List[str],
    gap_cost: int,
    end_gap_cost: int
) -> ty.Tuple[ty.List[ModuleSequence], ty.Dict[str, ty.List[str]]]:
    """Match the primary sequence exactly against the database.
    
    :param session: The Neo4j session.
    :type session: neo4j.Session
    :param match_items: The items to match against.
    :type match_items: ty.List[ty.Dict[str, ty.Any]]
    :param top_n: The number of top matches to return.
    :type top_n: int
    :param match_against_molecules: Whether to match against molecules.
    :type match_against_molecules: bool
    :param match_against_proto_clusters: Whether to match against proto-clusters.
    :type match_against_proto_clusters: bool
    :param selected_bioactivity_labels: The selected bioactivity labels.
    :type selected_bioactivity_labels: ty.List[str]
    :param gap_cost: The gap cost.
    :type gap_cost: int
    :param end_gap_cost: The end gap cost.
    :type end_gap_cost: int
    :return: The top sequences and their bioactivities.
    :rtype: ty.Tuple[ty.List[ModuleSequence], ty.Dict[str, ty.List[str]]]
    :raises ValueError: If no matches are found.
    """
    # Make sure the query is non-ambiguous.
    assert all([len(x) == 1 for x in match_items]), "Query is too ambiguous for exact matching!"

    # Exact matching.
    query = []
    for items in match_items:
        item = items[0]
        props = item["properties"]
        props["identifier"] = item["identifier"]
        query.append(props)
    seq_a = ModuleSequence("Query", parse_primary_sequence(query))

    # Retrieve all primary sequences, or only those coming from molecules or proto-clusters.
    if match_against_molecules and match_against_proto_clusters:
        result = session.run("""
            MATCH (b:PrimarySequence)
            WHERE (b)<-[:HAS_PRIMARY_SEQUENCE]-(:Compound) 
            OR (b)<-[:HAS_PRIMARY_SEQUENCE]-(:ProtoCluster)
            RETURN b""",
            fetch_size=1
        )
    elif match_against_molecules:
        result = session.run("MATCH (:Compound)-[:HAS_PRIMARY_SEQUENCE]->(b:PrimarySequence) RETURN b", fetch_size=1)
    elif match_against_proto_clusters:
        result = session.run(
            "MATCH (:ProtoCluster)-[:HAS_PRIMARY_SEQUENCE]->(b:PrimarySequence) RETURN b", 
            fetch_size=1
    )
    else:
        raise ValueError("No database to match against!")

    top = []
    for record in result: # Loop over all primary sequences in the database.
        id_b = record["b"]["identifier"]

        # Get bioactivity labels for id_b.
        bioactivities = []
        query = "MATCH (c:Compound {identifier: $identifier})-[:HAS_BIOACTIVITY]->(b:Bioactivity) RETURN b"
        result = session.run(query, identifier=id_b)
        for record in result:
            bioactivities.append(record["b"]["name"])

        # If selected bioactivity labels is empty, we are not going to filter by bioactivity.
        if len(selected_bioactivity_labels) != 0:
            if not all([x in bioactivities for x in selected_bioactivity_labels]):
                continue

        # Retrieve primary sequence.
        seq_b = retrieve_primary_sequence(session, id_b)

        # If the sequence is not empty, do the alignment and get the score.
        if len(seq_b) != 0:
            seq_b = ModuleSequence(id_b, parse_primary_sequence(seq_b))
            score = seq_a.optimal_alignment(seq_b, gap_cost, end_gap_cost).score

            # Keep top_n best scores.
            if len(top) < top_n:
                top.append((id_b, score, seq_b, bioactivities))
                top.sort(key=lambda x: x[1], reverse=True)
            else:
                if score > top[-1][1]:
                    top.pop()
                    top.append((id_b, score, seq_b, bioactivities))
                    top.sort(key=lambda x: x[1], reverse=True)
        else:
            continue

    # Check if top is empty.
    if len(top) == 0:
        raise ValueError("No matches found!")

    # Compile the bioactivities.
    bioactivites = {x[0]: x[3] for x in top}

    # Sort top on name.
    top.sort(key=lambda x: x[0]) # Exact organization of MSA is sensitive to order...

    # Compile the alignment.
    seqs = [seq_a] + [x[2] for x in top]

    return seqs, bioactivites

def match_pattern(
    session: neo4j.Session,
    match_items: ty.List[ty.Dict[str, ty.Any]],
    top_n: int,
    match_against_molecules: bool,
    match_against_proto_clusters: bool,
    selected_bioactivity_labels: ty.List[str],
    has_no_leading_modules: bool,
    has_no_trailing_modules: bool
) -> ty.Tuple[ty.List[ModuleSequence], ty.Dict[str, ty.List[str]]]:
    """Match the primary sequence pattern against the database.
    
    :param session: The Neo4j session.
    :type session: neo4j.Session
    :param match_items: The items to match against.
    :type match_items: ty.List[ty.Dict[str, ty.Any]]
    :param top_n: The number of top matches to return.
    :type top_n: int
    :param match_against_molecules: Whether to match against molecules.
    :type match_against_molecules: bool
    :param match_against_proto_clusters: Whether to match against proto-clusters.
    :type match_against_proto_clusters: bool
    :param selected_bioactivity_labels: The selected bioactivity labels.
    :type selected_bioactivity_labels: ty.List[str]
    :param has_no_leading_modules: Whether the sequence has no leading modules.
    :type has_no_leading_modules: bool
    :param has_no_trailing_modules: Whether the sequence has no trailing modules.
    :type has_no_trailing_modules: bool
    :return: The top sequences and their bioactivities.
    :rtype: ty.Tuple[ty.List[ModuleSequence], ty.Dict[str, ty.List[str]]]
    :raises ValueError: If no matches are found.
    """
    def compile_path_query(has_leading_modules: bool) -> str:
        """Compile the path query.
        
        :param has_leading_modules: Whether the sequence has leading modules.
        :type has_leading_modules: bool
        :return: The compiled path query.
        :rtype: str
        """
        query = []

        if has_leading_modules:
            module_range = range(2, len(match_items) + 1)
        else:
            module_range = range(1, len(match_items))

        for i in module_range:
            query.append(f"(u{i})-[:NEXT]->")

        if has_leading_modules:
            query.append(f"(u{len(match_items) + 1})")
        else:
            query.append(f"(u{len(match_items)})")

        query = "".join(query)

        return query

    def match_item_to_query(index: int, match_item: ty.Dict[str, ty.Any]) -> str:
        """Convert a match item to a query.
        
        :param index: The index of the match item.
        :type index: int
        :param match_item: The match item.
        :type match_item: ty.Dict[str, ty.Any]
        :return: The query addendum.
        :rtype: str
        """
        subqueries = []
        for option in match_item:
            if option["identifier"] == "polyketide":
                domains = option["properties"]["accessory_domains"]
                decoration = option["properties"]["decoration_type"]
                subquery = \
                    (f"(u{index}.identifier = 'polyketide'") + \
                    (f" AND u{index}.accessory_domains = {domains}" if domains else "") + \
                    (f" AND u{index}.decoration_type = '{decoration}')" if decoration is not None else ")")
                subqueries.append(subquery)
            elif option["identifier"] == "peptide":
                classification = option["properties"]["classification"]
                subquery = \
                    (f" (u{index}.identifier = 'peptide'") + \
                    (f" AND u{index}.classification = '{classification}')" if classification else ")")
                subqueries.append(subquery)
            else:
                pass

        subquery = " OR ".join(subqueries)
        return " AND (" + subquery + ")" if subquery != "" else ""

    if has_no_leading_modules:
        query = "MATCH (b:PrimarySequence)-[:START]->(u1) MATCH path = " + compile_path_query(False)
    else:
        query = "MATCH (b:PrimarySequence)-[:START]->(u1) MATCH path = (u1)-[:NEXT*]->" + compile_path_query(True)

    # Retrieve all primary sequences, or only those coming from molecules or proto-clusters.
    if match_against_molecules and match_against_proto_clusters:
        query  += " WHERE ((b)<-[:HAS_PRIMARY_SEQUENCE]-(:Compound) OR (b)<-[:HAS_PRIMARY_SEQUENCE]-(:ProtoCluster))"
    elif match_against_molecules:
        query  += " WHERE ((b)<-[:HAS_PRIMARY_SEQUENCE]-(:Compound))"
    elif match_against_proto_clusters:
        query  += " WHERE ((b)<-[:HAS_PRIMARY_SEQUENCE]-(:ProtoCluster))"
    else:
        raise ValueError("No database to match against!")

    if has_no_leading_modules and has_no_trailing_modules:
        query += f" AND (NOT ()-[:NEXT]->(u1) AND NOT (u{len(match_items)})-[:NEXT]->())"
    elif has_no_leading_modules:
        query += " AND NOT ()-[:NEXT]->(u1)"
    elif has_no_trailing_modules:
        query += f" AND NOT (u{len(match_items) + 1})-[:NEXT]->()"

    for i, match_item in enumerate(match_items):
        if has_no_leading_modules:
            index = i + 1
        else:
            index = i + 2

        query_addendum = match_item_to_query(index, match_item)
        if query_addendum != "":
            query += query_addendum

    query += " RETURN DISTINCT b"

    print(query)

    result = session.run(query, fetch_size=1)

    top = []
    for record in result:
        id_b = record["b"]["identifier"]

        # Get bioactivity labels for id_b.
        bioactivities = []
        query = "MATCH (c:Compound {identifier: $identifier})-[:HAS_BIOACTIVITY]->(b:Bioactivity) RETURN b"
        result = session.run(query, identifier=id_b)
        for record in result:
            bioactivities.append(record["b"]["name"])

        # If selected bioactivity labels is empty, we are not going to filter by bioactivity.
        if len(selected_bioactivity_labels) != 0:
            if not all([x in bioactivities for x in selected_bioactivity_labels]):
                continue

        seq_b = retrieve_primary_sequence(session, id_b)
        if len(seq_b) != 0:
            seq_b = ModuleSequence(id_b, parse_primary_sequence(seq_b))
            top.append((id_b, seq_b, bioactivities))
        if len(top) >= top_n: # Returns first top_n matches found.
            break

    # Check if top is empty.
    if len(top) == 0:
        raise ValueError("No matches found!")

    # Compile the bioactivities.
    bioactivites = {x[0]: x[2] for x in top}

    return [x[1] for x in top], bioactivites

blueprint_match_database = Blueprint("match_database", __name__)
@blueprint_match_database.route("/api/match_database", methods=["POST"])
def match_database() -> Response:
    """API endpoint for finding matches in the database.
    
    :return: The matches.
    :rtype: ResponseData
    """
    data = request.get_json()
    data = data["data"]

    # Unpack data.
    try:
        match_items = data["matchItems"]
        exact_matching = data["ambiguousNotAllowed"]
        selected_bioactivity_labels = data["selectedBioactivityLabels"]
        match_against_molecules = data["matchAgainstMolecules"]
        match_against_proto_clusters = data["matchAgainstProtoClusters"]
        top_n = data["numMatchesToReturn"]
        has_no_leading_modules = data["hasNoLeadingModules"]
        has_no_trailing_modules = data["hasNoTrailingModules"]
    except KeyError as err:
        message = f"Key not present in submission: {err}"
        return ResponseData(Status.Failure, message=message).to_dict()

    # Mount driver.
    try:
        # driver = neo4j.GraphDatabase.driver("bolt://database:7687")
        driver = neo4j.GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))
    except Exception as err:
        message = f"Could not mount neo4j driver: {err}"
        return ResponseData(Status.Failure, message=message).to_dict()

    # Mount driver and perform matching.
    with driver.session() as session:
        # Set gap costs for pairwise alignment and multiple sequence alignment.
        gap_cost = 2
        end_gap_cost = 1

        if exact_matching:
            # Exact matching was picked.
            try:
                seqs, bioactivities = match_exact(
                    session,
                    match_items,
                    top_n,
                    match_against_molecules,
                    match_against_proto_clusters,
                    selected_bioactivity_labels,
                    gap_cost,
                    end_gap_cost
                )
                assert len(seqs) != 0, "No matches found!"
                msa = MultipleSequenceAlignment(seqs, gap_cost, end_gap_cost).get_alignment()
            except Exception as err:
                message = f"Error during exact matching: {err}"
                return ResponseData(Status.Failure, message=message).to_dict()
        else:
            # Querying was picked.
            try:
                seqs, bioactivities = match_pattern(
                    session,
                    match_items,
                    top_n,
                    match_against_molecules,
                    match_against_proto_clusters,
                    selected_bioactivity_labels,
                    has_no_leading_modules,
                    has_no_trailing_modules
                )
                assert len(seqs) != 0, "No matches found!"
                msa = MultipleSequenceAlignment(seqs, gap_cost, end_gap_cost).get_alignment()
            except Exception as err:
                message = f"Error during querying: {err}"
                return ResponseData(Status.Failure, message=message).to_dict()

        # Compile matches.
        matches = []
        for seq in msa:

            # Compile sequence representation.
            try:
                seq_repr = []
                for x, _ in seq._seq:
                    if isinstance(x, PolyketideMotif):
                        x_repr = x.type.name + (str(x.decoration_type) if x.decoration_type is not None else "")
                    elif isinstance(x, PeptideMotif):
                        x_repr = "AA" + str(x.type.value)
                    elif x is Gap:
                        x_repr = "GAP"
                    else:
                        x_repr = "???"
                    seq_repr.append(x_repr)
            except Exception as err:
                message = f"Error during constructing sequence representation: {err}"
                return ResponseData(Status.Failure, message=message).to_dict()

            # Compile URL.
            try:
                if seq.name.startswith("NPA"):
                    url = "https://www.npatlas.org/explore/compounds/" + seq.name
                else:
                    url = None
            except Exception as err:
                message = f"Error during constructing URL: {err}"
                return ResponseData(Status.Failure, message=message).to_dict()

            # Compile result for this sequence.
            try:
                matches.append({
                    "identifier": seq.name,
                    "bioactivity": list(set(bioactivities[seq.name])) if seq.name != "Query" else [],
                    "sequence": seq_repr,
                    "url": url
                })
            except Exception as err:
                message = f"Error during compiling result: {err}"
                return ResponseData(Status.Failure, message=message).to_dict()

    return ResponseData(
        Status.Success,
        payload={"matches": matches},
        message="Matching completed!"
    ).to_dict()
