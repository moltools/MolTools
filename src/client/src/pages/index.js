import React from "react";
import Tile from "../components/Tile";

const Home = ({ tilesData }) => {
    return (
        <div>
            <h1>Home</h1>
            <div className="tile-container">
                {tilesData.length === 0 ? (
                    <p>No app tiles are available.</p>
                ) : (
                        tilesData.map((tile, index) => (
                            <Tile
                                key={index}
                                name={tile.name}
                                icon={tile.icon}
                                description={tile.description}
                                path={tile.path}
                            />
                        )   )
                )}
            </div>
        </div>
    );
};

export default Home;