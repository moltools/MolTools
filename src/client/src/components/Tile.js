import React from "react";
import { Link } from "react-router-dom";
import { FiImage } from "react-icons/fi";

const Tile = ({ name, icon, color, description, path }) => {

    const iconSize = 50;

    // Use the provided color or default to gray.
    const tileStyle = {backgroundColor: color || "#ccc"};

    return (
        <Link className="tile" to={path} style={tileStyle}>
            {typeof icon === 'string' ? (
                <img src={icon} alt={name} width={iconSize} height={iconSize} />
            ) : (
                icon ? (
                    React.cloneElement(icon, { size: iconSize })
                ) : (
                    <FiImage size={iconSize} />
                )
            )}
            <h3>{name}</h3>
            <p>{description}</p>
        </Link>
    );
};

export default Tile;