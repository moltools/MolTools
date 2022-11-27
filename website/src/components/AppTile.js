import React from "react";

export const AppTile = ({ name, logo, url }) => {
    return(
        <div className="apptile">
            <a href={url}>
                <img
                    className="img-fluid"
                    src={logo}
                    alt={name}
                />
            </a>
        </div>
    );
};
