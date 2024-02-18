import React, { useState } from "react";
import { toast } from "react-toastify";

const Overview = () => {
    return (
        <div>
            <h1 class="title">RetroMol</h1>
            <p>
                RetroMol is a tool for the generation of retrosynthetic trees for a given molecule, with the purpose
                of identifying potential precursors for the molecule. The scope of the tool is to provide a user-friendly
                interface for the discovery of retrobiosynthesis of modular natural products. Modular natural products
                currently supported by RetroMol include polyketides, non-ribosomal peptides, and polyketide-nonribosomal
                peptide hybrids.
            </p>
            <br />
            <p>
                You can use RetroMol from two starting inputs:
                <ol style={{paddingLeft: "40px"}}>
                    <li>A molecule represented by a SMILES string</li>
                    <li>An output file from <i>antiSMASH</i> containing mined proto-clusters</li>
                </ol>
            </p>
            <br />
            <p>
                Select the appropriate tab to get started. After parsing the molecule or proto-cluster, you can query
                the database to retrieve similar molecules and their retrosynthetic trees.
            </p>
        </div>
    );
};

const ParseMolecule = () => {
    const [smiles, setSmiles] = useState("");
    const [monomerGraph, setMonomerGraph] = useState({});

    const parseMolecule = async () => {
        try {
            const response = await fetch("/api/parse_molecule", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ "smiles": smiles })
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();
            
            if (json.status === "success") {
                setMonomerGraph(json.payload.monomer_graph);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            toast.error(error.message);
        }
    };

    return (
        <div class="column is-full">
            <div class="field">
                <div 
                    class="control"
                >
                    <input class="input" type="text" placeholder="Enter SMILES" onChange={(e) => setSmiles(e.target.value)} />
                </div>
            </div>
            <div class="control">
                <button class="button is-link is-light" style={{marginRight: "5px"}} onClick={parseMolecule}>Parse</button>
                <button class="button is-link is-light" onClick={() => setMonomerGraph({})}>Clear</button>
            </div>
            <div>
                
            </div>
        </div>
    );
};

const ParseProtoCluster = () => {
    return (
        <div></div>
    );
};

const EmbedResults = () => {
    return (
        <div></div>
    );
};

const RetroMol = () => {
    const tabs = ["Overview", "Parse molecule", "Parse proto-cluster", "Query database"];
    const [selectedTab, setSelectedTab] = useState("Overview");

    return (
        <div>
            <div class="tabs is-boxed" style={{padding: "20px"}}>
                <ul>
                    {tabs.map((tab, index) => (
                        <li 
                            key={index}
                            class={selectedTab === tab ? "is-active" : ""}
                            onClick={() => setSelectedTab(tab)}
                        >
                            <a>{tab}</a>
                        </li>
                    ))}
                </ul>
            </div>
            <div class="container">
                {selectedTab === "Overview" && <Overview />}
                {selectedTab === "Parse molecule" && <ParseMolecule />}
                {selectedTab === "Parse proto-cluster" && <ParseProtoCluster />}
                {selectedTab === "Embed results" && <EmbedResults />}
            </div>
        </div>
    );
};

export default RetroMol;