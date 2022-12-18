import React from "react";
import { AppTile } from "../components/AppTile";
import "../index.css";

import undecided_logo from "../logos/undecided.png";
import retromol_logo from "../logos/retromol.png";
import cinemol_logo from "../logos/cinemol.png";

function Home() {
    return (
        <div>
            <div className="apptiles">
                <AppTile name="retromol" logo={retromol_logo} url="https://retromol.bioinformatics.nl/"/>
                <AppTile name="cinemol" logo={cinemol_logo} url="/cinemol"/>
                <AppTile name="smilescorrector" logo={undecided_logo} url="/smilescorrector"/>
                <AppTile name="activitymol" logo={undecided_logo} url="/activitymol"/>
            </div>
        </div>
    );
};

export default Home;
