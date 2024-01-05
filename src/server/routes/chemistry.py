import typing as ty

from flask import Blueprint, Response, request 
from rdkit import Chem 
from rdkit.Chem.Draw import rdMolDraw2D

from .common import Status, ResponseData

def mol_to_svg(
    mol: Chem.Mol,
    window_size: ty.Tuple[int, int] = (450, 350),
    background_color: ty.Optional[str] = None
) -> str:
    """
    Draw a molecule with its substructures highlighted.

    :param Chem.Mol mol: The molecule to draw.
    :param ty.Tuple[int, int] window_size: The size of the window.
    :param ty.Optional[str] background_color: The background color.
    :return: The SVG string.
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

blueprint_draw_smiles = Blueprint("draw_smiles", __name__)
@blueprint_draw_smiles.route("/api/draw_smiles", methods=["POST"])
def draw_smiles() -> Response:
    """
    API endpoint for drawing a SMILES string and returning the SVG string.

    :param str smiles: The SMILES string.
    :return: The SVG string.
    :rtype: ResponseData
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
            
            # Draw the molecule.
            svg_string = mol_to_svg(mol)
            payload = {"svg_string": svg_string}

            msg = "Successfully parsed SMILES string!"
            return ResponseData(Status.Success, payload, msg).to_dict()
        
        except Exception:
            msg = "Could not parse SMILES string!"
            return ResponseData(Status.Failure, message=msg).to_dict()