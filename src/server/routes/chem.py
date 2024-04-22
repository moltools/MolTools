"""This module contains the API endpoints for general chemistry-related tasks."""
import typing as ty

from flask import Blueprint, Response, request
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from .common import Status, ResponseData

def draw_mol(
    mol: Chem.Mol,
    window_size: ty.Tuple[int, int] = (800, 800),
    background_color: ty.Optional[str] = None
) -> str:
    """
    Draw a molecule with its substructures highlighted.

    :param mol: RDKit molecule object.
    :type mol: Chem.Mol
    :param window_size: Size of the window to draw the molecule in.
    :type window_size: ty.Tuple[int, int]
    :param background_color: Background color of the window.
    :type background_color: ty.Optional[str]
    :return: SVG representation of the molecule.
    :rtype: str
    """
    drawing = rdMolDraw2D.MolDraw2DSVG(*window_size)

    options = drawing.drawOptions()
    if background_color is not None:
        options.setBackgroundColour(background_color)
    options.useBWAtomPalette()

    drawing.DrawMolecule(mol)

    drawing.FinishDrawing()
    svg_str = drawing.GetDrawingText().replace("svg:", "")

    return svg_str

blueprint_smiles_to_svg = Blueprint("smiles_to_svg", __name__)
@blueprint_smiles_to_svg.route("/api/smiles_to_svg", methods=["POST"])
def parse_smiles_to_svg() -> Response:
    """API endpoint for converting SMILES to SVG.
    
    :return: SVG of the input SMILES string.
    :rtype: ResponseData
    """
    data = request.get_json()
    data = data["data"]

    # Unpack the SMILES string from the request data.
    try:
        smiles = data["smiles"]
        height = data["height"]
        width = data["width"]
    except KeyError as err:
        message = f"Key not present in submission: {err}"
        return ResponseData(Status.Failure, message=message).to_dict()

    # Attempt to parse the SMILES string.
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as err:
        message = f"Failed to parse SMILES string: {err}"
        return ResponseData(Status.Failure, message=message).to_dict()

    # Draw the molecule and return the SVG.
    svg = draw_mol(mol, window_size=(height, width))
    payload = {"svg": svg}
    message = "Successfully converted SMILES to SVG."
    return ResponseData(Status.Success, payload=payload, message=message).to_dict()
