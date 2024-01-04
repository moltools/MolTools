from cinemol.parsers import parse_sdf
from cinemol.chemistry import Style, Look, draw_molecule
from flask import Blueprint, Response, request 

from .common import Status, ResponseData

blueprint_draw_model = Blueprint("draw_model", __name__)
@blueprint_draw_model.route("/api/draw_model", methods=["POST"])
def draw_smiles() -> Response:
    """
    API endpoint for drawing a 3D model and returning the SVG string.

    :param str smiles: SDF string.
    :return: The SVG string.
    :rtype: str
    """
    # Parse request data.
    data = request.get_json()
    
    # Get SDF string.
    sdf_str = data.get("sdfString", None)
    
    if sdf_str is None:
        msg = "No SDF string provided!"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    else:
        try:
            # Parse SDF string; return the first molecule, if any.
            atoms, bonds = parse_sdf(sdf_str)

            svg_str = draw_molecule(atoms, bonds, Style.Tube, Look.Glossy, 25)

            return ResponseData(Status.Success, svg_str).to_dict()
        
        except Exception as e:
            msg = "Failed to parse SDF string!"
            return ResponseData(Status.Failure, message=msg).to_dict()