import React, { useState } from "react";

import Overview from "../components/retromol/Overview";
import ParseMolecule from "../components/retromol/ParseMolecule";
import ParseProtoCluster from "../components/retromol/ParseProtoCluster";
import QueryDatabase from "../components/retromol/QueryDatabase";

const RetroMol = () => {

    // Define tabs.
    const tabs = [
        "Overview", 
        "Parse molecule", 
        "Parse proto-cluster", 
        "Query database"
    ];

    // Define state for selected tab.
    const [selectedTab, setSelectedTab] = useState("Overview");

    return (
        <div>
            <div className="tabs is-boxed" style={{ padding: "20px" }}>
                <ul>
                    {tabs.map((tab, index) => (
                        <li 
                            key={index}
                            className={selectedTab === tab ? "is-active" : ""}
                            onClick={() => setSelectedTab(tab)}

                        >
                            <a>{tab}</a>
                        </li>
                    ))}
                </ul>
            </div>
            <div style={{ padding: "20px", paddingTop: "0px"}}>
                {selectedTab === "Overview" && <Overview />}
                {selectedTab === "Parse molecule" && <ParseMolecule />}
                {selectedTab === "Parse proto-cluster" && <ParseProtoCluster />}
                {selectedTab === "Query database" && <QueryDatabase />}
            </div>
        </div>
    );
};

export default RetroMol;