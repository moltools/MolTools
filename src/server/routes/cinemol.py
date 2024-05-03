"""This module contains the API endpoints for CineMol."""
from cinemol.parsers import parse_sdf
from cinemol.api import Style, Look, draw_molecule
# from cinemol.version import version
from flask import Blueprint, Response, request

from .common import Status, ResponseData

blueprint_draw_model = Blueprint("draw_model", __name__)
@blueprint_draw_model.route("/api/draw_model", methods=["POST"])
def draw_model() -> Response:
    """API endpoint for drawing a 3D model and returning the SVG string.

    :param str smiles: SDF string.
    :return: The SVG string.
    :rtype: ResponseData
    """
    # Parse request data.
    data = request.get_json()

    # Get SDF string.
    sdf_str = data.get("sdf_string", None)

    # Get parameters.
    style = {
        "space-filling":    Style.SpaceFilling,
        "ball-and-stick":   Style.BallAndStick,
        "tube":             Style.Tube,
        "wireframe":        Style.Wireframe
    }[data.get("style", "tube")]

    look = {
        "cartoon":           Look.Cartoon,
        "glossy":            Look.Glossy,           
    }[data.get("look", "glossy")]

    include_hydrogens = data.get("include_hydrogens", False)
    resolution = int(data.get("resolution", 25))
    rotation_x = float(data.get("rotation_x", 0.0))
    rotation_y = float(data.get("rotation_y", 0.0))
    rotation_z = float(data.get("rotation_z", 0.0))

    view_box = data.get("view_box", None)
    if view_box is not None:
        view_box = (
            float(view_box["min_x"]),
            float(view_box["min_y"]),
            float(view_box["width"]),
            float(view_box["height"])
        )

    # Get window size.
    width = data.get("width", None)
    height = data.get("height", None)

    if width is not None and height is not None:
        window = (int(width), int(height))
    else:
        window = None

    # Draw molecule.
    if sdf_str is None:
        msg = "No SDF string provided!"
        return ResponseData(Status.Failure, message=msg).to_dict()

    else:
        try:
            # Parse SDF string; return the first molecule, if any.
            atoms, bonds = parse_sdf(sdf_str)

            # Draw molecule.
            svg = draw_molecule(
                atoms=atoms,
                bonds=bonds,
                style=style,
                look=look,
                resolution=resolution,
                window=window,
                view_box=view_box,
                rotation_over_x_axis=rotation_x,
                rotation_over_y_axis=rotation_y,
                rotation_over_z_axis=rotation_z,
                scale=1.0,
                focal_length=100.0,
                exclude_atoms=None if include_hydrogens else ["H"],
                verbose=False
            )

            svg_str = svg.to_svg()
            vb = svg.view_box

            # Return the SVG string.
            payload = {
                "svg_string": svg_str, 
                "view_box": {
                    "min_x": vb.min_x, 
                    "min_y": vb.min_y, 
                    "width": vb.width, 
                    "height": vb.height
                }
            }

            return ResponseData(Status.Success, payload).to_dict()

        except Exception as e:
            msg = f"Failed to parse SDF string: {e}"
            return ResponseData(Status.Failure, message=msg).to_dict()

blueprint_fetch_cinemol_version = Blueprint("fetch_cinemol_version", __name__)
@blueprint_fetch_cinemol_version.route("/api/fetch_cinemol_version", methods=["GET"])
def fetch_cinemol_version() -> Response:
    """API endpoint for fetching the version of Cinemol.

    :return: The version of Cinemol.
    :rtype: ResponseData
    """
    # Get version.
    payload = {"version": "0.1.0-dev"}

    # Return the version.
    msg = "Successfully fetched Cinemol version!"

    return ResponseData(Status.Success, payload, msg).to_dict()
