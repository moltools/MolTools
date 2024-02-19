import React, { useState } from "react";
import { toast } from "react-toastify";

const Overview = () => {
    return (
        <div>
            <h1 className="title">RetroMol</h1>
            <p>
                RetroMol is a tool for the generation of retrosynthetic trees for a given molecule, with the purpose
                of identifying potential precursors for the molecule. The scope of the tool is to provide a user-friendly
                interface for the discovery of retrobiosynthesis of modular natural products. Modular natural products
                currently supported by RetroMol include polyketides, non-ribosomal peptides, and polyketide-nonribosomal
                peptide hybrids.
            </p>
            <br />
            <div>
                You can use RetroMol from two starting inputs:
                <ol style={{paddingLeft: "40px"}}>
                    <li>A molecule represented by a SMILES string</li>
                    <li>An output file from <i>antiSMASH</i> containing mined proto-clusters</li>
                </ol>
            </div>
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
    const [svgString, setSvgString] = useState("");
    const [monomerGraph, setMonomerGraph] = useState({});
    const [isLoading, setIsLoading] = useState(false);

    React.useEffect(() => {
        const drawSmiles = async () => {
            try {
                const response = await fetch("/api/draw_smiles", {
                    method: "POST",
                    headers: {"Content-Type": "application/json"},
                    body: JSON.stringify({ "smiles": smiles })
                });

                if (!response.ok) {
                    throw new Error("Network response was not ok!");
                };

                const json = await response.json();
                
                // Unpack response.
                if (json.status === "success") {
                    setSvgString(json.payload.svg_string);
                } else if (json.status === "warning") {
                    toast.warn(json.message);
                } else if (json.status === "failure") {
                    toast.error(json.message);
                };
        
            } catch (error) {
                const msg = "Could not draw molecule!";
                toast.error(msg, { autoClose: true });
                console.error(error);
            };
        };

        if (smiles) {
            drawSmiles();
        };
    }, [smiles]);

    const parseMolecule = async () => {
        setIsLoading(true);
        try {
            const response = await fetch("/api/parse_retromol", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ "smiles": smiles })
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();
            
            if (json.status === "success") {
                toast.success(json.message);
                console.log(json.payload.monomer_graph)
                setMonomerGraph(json.payload.monomer_graph);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            toast.error(error.message);
        }
        setIsLoading(false);
    };

    // Load example -- sets the SMILES input to a default value
    const loadExample = () => {
        setSmiles("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O");
    };

    const clear = () => {
        setSmiles("");
        setMonomerGraph({});
        setSvgString("");
    };

    return (
        <div className="column is-full">
            <div className="field">
                <div 
                    className="control"
                >
                    <input className="input" value={smiles} type="text" placeholder="Enter SMILES" onChange={(e) => setSmiles(e.target.value)} disabled={isLoading} />
                </div>
            </div>
            <div className="control">
                <button className="button is-link is-light" style={{marginRight: "5px"}} onClick={loadExample} disabled={isLoading}>Example</button>
                <button className="button is-link is-light" style={{marginRight: "5px"}} onClick={parseMolecule} disabled={isLoading}>Parse</button>
                <button className="button is-link is-light" onClick={clear} disabled={isLoading}>Clear</button>
            </div>
            <div className="columns" style={{marginTop: "5px"}}>
                <div className="column has-text-centered" style={{margin: "5px"}}>
                    <div dangerouslySetInnerHTML={{ __html: svgString }} />
                </div>
                <div className="column has-text-centered" style={{margin: "5px"}}>
                    <pre>{JSON.stringify(monomerGraph, null, 2)}</pre>
                </div>
            </div>
        </div>
    );
};

const ParseProtoCluster = () => {
    return (
        <div className="column is-full">
            Not implemented yet
        </div>
    );
};

const EmbedResults = () => {
    return (
        <div className="column is-full">
            Not implemented yet
        </div>
    );
};

const QueryDatabase = () => {
    return (
        <div className="column is-full">
            Not implemented yet
        </div>
    );
}

const RetroMol = () => {
    const tabs = ["Overview", "Parse molecule", "Parse proto-cluster", "Embed results", "Query database"];
    const [selectedTab, setSelectedTab] = useState("Overview");

    return (
        <div>
            <div className="tabs is-boxed" style={{padding: "20px"}}>
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
            <div className="container">
                {selectedTab === "Overview" && <Overview />}
                {selectedTab === "Parse molecule" && <ParseMolecule />}
                {selectedTab === "Parse proto-cluster" && <ParseProtoCluster />}
                {selectedTab === "Embed results" && <EmbedResults />}
                {selectedTab === "Query database" && <QueryDatabase />}
            </div>
        </div>
    );
};

export default RetroMol;