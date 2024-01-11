from flask import Blueprint, Response, request 
from pikachu.general import read_smiles, svg_string_from_structure
from pikachu.drawing.drawing import Drawer, Options

from .common import Status, ResponseData

blueprint_draw_smiles_with_pikachu = Blueprint("draw_smiles_with_pikachu", __name__)
@blueprint_draw_smiles_with_pikachu.route("/api/draw_smiles_with_pikachu", methods=["POST"])
def draw_smiles_with_pikachu() -> Response:
    """
    API endpoint for drawing a SMILES string with PIKAChU and returning the SVG string.

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
            structure = read_smiles(smiles)
            drawer = Drawer(structure, options=Options(), coords_only=True)

            x_pos = [atom.draw.position.x for atom in drawer.drawn_atoms]
            y_pos = [atom.draw.position.y for atom in drawer.drawn_atoms]
            margin = 25
            x_min = min(x_pos) - margin
            y_min = min(y_pos) - margin
            x_max = max(x_pos) + margin
            y_max = max(y_pos) + margin
            width = max(0, x_max - x_min)
            height = max(0, y_max - y_min)
            
            svg_string = f"""<svg width="{400}" height="{400}" viewBox="{x_min} {y_min} {width} {height}" xmlns="http://www.w3.org/2000/svg">"""
            svg_string += drawer.draw_svg()
            svg_string += "</svg>"

            # svg_string = svg_from_smiles(smiles)
            payload = {"svg_string": svg_string}
            msg = "Successfully parsed SMILES string!"
            return ResponseData(Status.Success, payload, msg).to_dict()
        
        except Exception as e:
            msg = "Failed to parse SMILES string!"
            return ResponseData(Status.Failure, message=msg).to_dict()
