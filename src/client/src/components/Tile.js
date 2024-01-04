import React from "react";
import { Link } from "react-router-dom";
import { FiImage } from "react-icons/fi";

const Tile = ({ name, icon, description, path }) => {
    const iconSize = 30;

    return (
        <Link
            className="tile"
            to={path}
        >
            {icon ? (
                React.cloneElement(icon, { size: iconSize })
            ) : (
                <FiImage size={iconSize} />
            )}
            <h3>{name}</h3>
            <p>{description}</p>
        </Link>
    );
};

export default Tile;