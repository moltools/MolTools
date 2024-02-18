import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";

const Overview = () => {
    return (
        <div>
            <h1 class="title">PIKAChU</h1>
            <p>
                The Python-based Informatics Kit for the Analysis of Chemical Units (PIKAChU) is cheminformatics toolkit 
                fully written in Python, without any external dependencies. PIKAChU is designed to be a lightweight, 
                easy-to-use, and efficient toolkit for the analysis of chemical units, such as molecules and reactions.
                It is able to draw a molecule using only the SMILES string representation of the molecule.
            </p>
        </div>
    );
};

const DrawMolecule = () => {
    const [smiles, setSmiles] = useState("");
    const [svgString, setSvgString] = useState("");

    const handleDownloadSvgString = () => {
        if (svgString.length > 0) {
            const blob = new Blob([svgString], { type: "image/svg+xml" });
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement("a");
            a.href = url;
            a.download = "pikachu.svg";
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            document.body.removeChild(a);
        } else {
            toast.error("No SVG to download!", { autoClose: true });
        }
    };

    React.useEffect(() => {
        const drawSmiles = async () => {
            try {
                const response = await fetch("/api/draw_smiles_with_pikachu", {
                    method: "POST",
                    headers: {"Content-Type": "application/json"},
                    body: JSON.stringify({ "smiles": smiles })
                });

                if (!response.ok) {
                    throw new Error("Network response was not ok!");
                };

                const json = await response.json();
                
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
                <button class="button is-link is-light" style={{marginRight: "5px"}} onClick={handleDownloadSvgString}>Download</button>
                <button class="button is-link is-light" onClick={() => setSvgString("")}>Clear</button>
            </div>
            <div>
                <div dangerouslySetInnerHTML={{ __html: svgString }} />
            </div>
        </div>
    );
};      

const PIKAChU = () => {
    const tabs = ["Overview", "Draw molecule"];
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
                {selectedTab === "Draw molecule" && <DrawMolecule />}
            </div>
        </div>
    );
};

export default PIKAChU;