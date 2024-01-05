import React from "react";
import { Link } from "react-router-dom";
import { FiImage } from "react-icons/fi";

/**
 * Tile component that represents a clickable tile with an icon, name, and description.
 * @param {Object} props - The props for the Tile component.
 * @param {string} props.name - The name of the tile.
 * @param {string | React.ReactNode} props.icon - The icon for the tile, either as a string (image URL) or a React component.
 * @param {string} props.color - The background color of the tile (optional, defaults to gray).
 * @param {string} props.description - The description of the tile.
 * @param {string} props.path - The URL path the tile should navigate to when clicked.
 * @returns {JSX.Element} - The rendered Tile component.
 */
const Tile = (props) => {
    // Destructure the props.
    const { name, icon, color, description, path } = props;

    // Define the size for the icon.
    const iconSize = 50;

    // Use the provided color or default to gray.
    const tileStyle = { backgroundColor: color || "#ccc" };

    return (
        <Link className="tile" to={path} style={tileStyle}>
            {typeof icon === "string" ? (
                <img src={icon} alt={name} width={iconSize} height={iconSize} />
            ) : (
                icon ? (
                    // Clone the provided React icon element and set the size.
                    React.cloneElement(icon, { size: iconSize })
                ) : (
                    // Use the default FiImage icon with the specified size.
                    <FiImage size={iconSize} />
                )
            )}
            <h3>{name}</h3>
            <p>{description}</p>
        </Link>
    );
};

// Export the Tile component as the default export.
export default Tile;