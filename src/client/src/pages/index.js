import React from "react";
import Tile from "../components/Tile";

// =====================================================================================================================
// Home component.
// =====================================================================================================================

/**
 * Home Component
 *
 * This component represents the home page of the MolTools application, displaying a collection of tiles based on the
 * provided 'tilesData' prop.
 *
 * @param {Object} props - The props for the Home component.
 * @param {Array} props.tilesData - An array of tile data objects to be displayed.
 * @returns {JSX.Element} The rendered Home component displaying the tiles.
 */
const Home = (props) => {
    // Destructure the props.
    const { tilesData } = props;

    return (
        <div>
            <div className="tile-container">
                {tilesData.length === 0 ? (
                    <p>No app tiles are available.</p>
                ) : (
                        tilesData.map((tile, index) => (
                            <Tile
                                key={index}
                                name={tile.name}
                                icon={tile.icon}
                                color={tile.color}
                                description={tile.description}
                                path={tile.path}
                            />
                        ))
                )}
            </div>
        </div>
    );
};

// Export the Home component as the default export.
export default Home;