import typing as ty

from cinemol.geometry import Point3D, Line3D, Sphere, Cylinder, CylinderCapType, get_perpendicular_lines
from cinemol.style import Color, Cartoon, Glossy, CoreyPaulingKoltungAtomColor as CPK, PubChemAtomRadius
from cinemol.model import ModelCylinder, ModelSphere, ModelWire, Scene 
from cinemol.parsers import parse_sdf
from cinemol.chemistry import Atom, Bond, Style, Look
from flask import Blueprint, Response, request 

from .common import Status, ResponseData

def draw_molecule(
    atoms: ty.List[Atom],
    bonds: ty.List[Bond],
    style: Style,
    look: Look,
    resolution: int, 
    rotation_over_x_axis: float = 0.0,
    rotation_over_y_axis: float = 0.0,
    rotation_over_z_axis: float = 0.0,
    verbose: bool = False,
    scale: float = 1.0,
    focal_length: float = 1.0
) -> str:
    """
    Draw a molecule using the given atoms and bonds.
    
    :param ty.List[Atom] atoms: The atoms to draw.
    :param ty.List[Bond] bonds: The bonds to draw.
    :param Style style: The style of the depiction.
    :param Look look: The look of the depiction.
    :param int resolution: The resolution of the SVG model.
    :param float rotation_over_x_axis: The rotation over the x-axis.
    :param float rotation_over_y_axis: The rotation over the y-axis.
    :param float rotation_over_z_axis: The rotation over the z-axis.
    :param bool verbose: Whether to print verbose output.
    :param float scale: The scale of the model and viewbox.
    :param float focal_length: The focal length of the camera.
    :return: The drawn molecule as an SVG string.
    :rtype: str
    """
    scene = Scene()
    scene.nodes = []
    atom_map = {atom.index: atom for atom in atoms}

    # Default settings for drawing.
    include_spheres = False
    include_cylinders = False
    include_wires = False
    calculate_sphere_sphere_intersections = False
    calculate_sphere_cylinder_intersections = False
    calculate_cylinder_sphere_intersections = False
    calculate_cylinder_cylinder_intersections = False
    
    # Default settings for cartoon look.
    stroke_color = Color(0, 0, 0)
    stroke_width = 0.05

    # ==========================================================================
    # Wireframe style
    # ==========================================================================

    # Wire style has a separate implementation that is faster than for geometric shapes.
    if style == Style.Wireframe:
        include_wires = True 
        wire_width = 0.05

        for bond in bonds:
            # Get start atom specifications.
            start_atom = atom_map[bond.start_index]
            start_symbol = start_atom.symbol
            start_color = CPK().get_color(start_symbol)
            start_pos = Point3D(*start_atom.coordinates)

            # Get end atom specifications.
            end_atom = atom_map[bond.end_index]
            end_symbol = end_atom.symbol
            end_color = CPK().get_color(end_symbol)
            end_pos = Point3D(*end_atom.coordinates)
            
            # Determine color of bond.
            if bond.color is not None:
                start_color = Color(*bond.color)
                end_color = Color(*bond.color)
            else:
                if start_atom.color is not None: start_color = Color(*start_atom.color)
                if end_atom.color is not None: end_color = Color(*end_atom.color)
            
            # If colors are not the same we split the bond down the middle, and 
            # draw two separate wires to represent the bond.
            if start_color != end_color:
                middle_pos = Point3D(
                    (start_pos.x + end_pos.x) / 2,
                    (start_pos.y + end_pos.y) / 2,
                    (start_pos.z + end_pos.z) / 2
                )
                scene.add_node(ModelWire(Line3D(start_pos, middle_pos), start_color, wire_width, bond.opacity))
                scene.add_node(ModelWire(Line3D(middle_pos, end_pos), end_color, wire_width, bond.opacity))
            else:
                scene.add_node(ModelWire(Line3D(start_pos, end_pos), start_color, wire_width, bond.opacity))

    # ==========================================================================
    # Space-filling style
    # ==========================================================================

    elif style == Style.SpaceFilling:
        include_spheres = True
        calculate_sphere_sphere_intersections = True
        
        for atom in atoms:
            # Get atom specifications.
            atom_symbol = atom.symbol
            atom_color = CPK().get_color(atom_symbol) if atom.color is None else Color(*atom.color)
            atom_radius = PubChemAtomRadius().to_angstrom(atom_symbol)
            atom_pos = Point3D(*atom.coordinates)

            # Determine atom look.
            if look == Look.Cartoon: depiction = Cartoon(atom_color, stroke_color, stroke_width, atom.opacity)
            elif look == Look.Glossy: depiction = Glossy(atom_color, atom.opacity)
            else: raise ValueError(f"Unknown look: '{look}'")

            # Add atom to scene.
            scene.add_node(ModelSphere(Sphere(atom_pos, atom_radius), depiction))

    # ==========================================================================
    # Ball-and-stick and tube styles
    # ==========================================================================

    elif style == Style.BallAndStick or style == Style.Tube:

        if style == Style.BallAndStick:
            include_spheres = True
            include_cylinders = True
            calculate_cylinder_sphere_intersections = True
            cap_type = CylinderCapType.NoCap 

            # Atoms are only drawn for the ball-and-stick style.
            for atom in atoms:
                # Get atom specifications.
                atom_symbol = atom.symbol
                atom_color = CPK().get_color(atom_symbol) if atom.color is None else Color(*atom.color)
                atom_radius = PubChemAtomRadius().to_angstrom(atom_symbol) / 3.0
                atom_pos = Point3D(*atom.coordinates)

                # Determine atom look.
                if look == Look.Cartoon: depiction = Cartoon(atom_color, stroke_color, stroke_width, atom.opacity)
                elif look == Look.Glossy: depiction = Glossy(atom_color, atom.opacity)
                else: raise ValueError(f"Unknown look: '{look}'")

                # Add atom to scene.
                scene.add_node(ModelSphere(Sphere(atom_pos, atom_radius), depiction))
                
        elif style == Style.Tube:
            include_cylinders = True
            calculate_cylinder_cylinder_intersections = True
            cap_type = CylinderCapType.Round

        for bond in bonds:
            # Get start atom specifications.
            start_atom = atom_map[bond.start_index]
            start_symbol = start_atom.symbol
            start_color = CPK().get_color(start_symbol)
            start_pos = Point3D(*start_atom.coordinates)

            # Get end atom specifications.
            end_atom = atom_map[bond.end_index]
            end_symbol = end_atom.symbol
            end_color = CPK().get_color(end_symbol)
            end_pos = Point3D(*end_atom.coordinates)

            # Determine color of bond.
            if bond.color is not None:
                start_color = Color(*bond.color)
                end_color = Color(*bond.color)
            else:
                if start_atom.color is not None: start_color = Color(*start_atom.color)
                if end_atom.color is not None: end_color = Color(*end_atom.color)

            # Determine number of cylinders to draw for each bond.
            bond_order = bond.order if style == Style.BallAndStick else 1 # Bond order is only used for ball-and-stick style.
            bond_radius = bond.radius if bond.radius is not None else 0.2
            temp_bond_radius = bond_radius / bond_order
            line = Line3D(start_pos, end_pos)
            lines = get_perpendicular_lines(line, temp_bond_radius * (bond_order + 1), bond_order)

            # Add cylinders to scene.
            for line in lines:
                if start_color != end_color:
                    middel_pos = Point3D(
                        (line.start.x + line.end.x) / 2,
                        (line.start.y + line.end.y) / 2,
                        (line.start.z + line.end.z) / 2
                    )
                    
                    # First part of bond.
                    if look == Look.Cartoon: depiction = Cartoon(start_color, stroke_color, stroke_width, bond.opacity)
                    elif look == Look.Glossy: depiction = Glossy(start_color, bond.opacity)
                    else: raise ValueError(f"Unknown look: '{look}'")
                    scene.add_node(ModelCylinder(Cylinder(line.start, middel_pos, temp_bond_radius, cap_type), depiction))

                    # Second part of bond.
                    if look == Look.Cartoon: depiction = Cartoon(end_color, stroke_color, stroke_width, bond.opacity)
                    elif look == Look.Glossy: depiction = Glossy(end_color, bond.opacity)
                    else: raise ValueError(f"Unknown look: '{look}'")
                    scene.add_node(ModelCylinder(Cylinder(middel_pos, line.end, temp_bond_radius, cap_type), depiction))

                else:
                    # Bond as a whole.
                    if look == Look.Cartoon: depiction = Cartoon(start_color, stroke_color, stroke_width, bond.opacity)
                    elif look == Look.Glossy: depiction = Glossy(start_color, bond.opacity)
                    else: raise ValueError(f"Unknown look: '{look}'")
                    scene.add_node(ModelCylinder(Cylinder(line.start, line.end, temp_bond_radius, cap_type), depiction))

    else:
        raise ValueError(f"Unknown style: '{style}'")

    # Draw scene.
    return scene.draw(
        resolution=resolution, 
        verbose=verbose, 
        rotation_over_x_axis=rotation_over_x_axis,
        rotation_over_y_axis=rotation_over_y_axis,
        rotation_over_z_axis=rotation_over_z_axis,
        include_spheres=include_spheres,
        include_cylinders=include_cylinders,
        include_wires=include_wires,
        calculate_sphere_sphere_intersections=calculate_sphere_sphere_intersections,
        calculate_sphere_cylinder_intersections=calculate_sphere_cylinder_intersections,
        calculate_cylinder_sphere_intersections=calculate_cylinder_sphere_intersections,
        calculate_cylinder_cylinder_intersections=calculate_cylinder_cylinder_intersections,
        scale=scale,
        focal_length=focal_length
    )

blueprint_draw_model = Blueprint("draw_model", __name__)
@blueprint_draw_model.route("/api/draw_model", methods=["POST"])
def draw_model() -> Response:
    """
    API endpoint for drawing a 3D model and returning the SVG string.

    :param str smiles: SDF string.
    :return: The SVG string.
    :rtype: str
    """
    # Parse request data.
    data = request.get_json()

    # Get SDF string.
    sdf_str = data.get("sdf_string", None)

    # Get parameters.
    style = {
        "Space-filling":    Style.SpaceFilling,
        "Ball-and-stick":   Style.BallAndStick,
        "Tube":             Style.Tube,
        "Wireframe":        Style.Wireframe
    }[data.get("style", "tube")]

    look = {
        "Cartoon":           Look.Cartoon,
        "Glossy":            Look.Glossy,           
    }[data.get("look", "glossy")]

    include_hydrogens = data.get("include_hydrogens", False)
    resolution = int(data.get("resolution", 25))
    rotation_x = float(data.get("rotation_x", 0.0))
    rotation_y = float(data.get("rotation_y", 0.0))
    rotation_z = float(data.get("rotation_z", 0.0))
    
    # Draw molecule.
    if sdf_str is None:
        msg = "No SDF string provided!"
        return ResponseData(Status.Failure, message=msg).to_dict()
    
    else:
        try:
            # Parse SDF string; return the first molecule, if any.
            atoms, bonds = parse_sdf(sdf_str)
            
            # Filter out hydrogens, if requested.
            if not include_hydrogens:
                filtered_atoms, filtered_bonds, inds_with_hs = [], [], []
                for atom in atoms:
                    if atom.symbol == "H": inds_with_hs.append(atom.index)
                    else: filtered_atoms.append(atom)

                for bond in bonds:
                    if bond.start_index not in inds_with_hs and bond.end_index not in inds_with_hs:
                        filtered_bonds.append(bond)
                        
                atoms, bonds = filtered_atoms, filtered_bonds

            svg_str = draw_molecule(
                atoms, 
                bonds, 
                style, 
                look, 
                resolution, 
                scale=1.0, 
                focal_length=100.0,
                rotation_over_x_axis=rotation_x,
                rotation_over_y_axis=rotation_y,
                rotation_over_z_axis=rotation_z
            )
            
            # Return the SVG string.
            payload = {"svg_string": svg_str}

            return ResponseData(Status.Success, payload).to_dict()
        
        except Exception as e:
            msg = "Failed to parse SDF string!"
            return ResponseData(Status.Failure, message=msg).to_dict()