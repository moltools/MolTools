import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";

const SmilesInput = ({ 
    smiles, 
    setSmiles, 
    parseMolecule, 
    handleRefresh,
    toggleReviewModal,
}) => {
    const [inputSvg, setInputSvg] = useState("");

    // Load an example SMILES string for erythromycin.
    const loadExample = () => {
        setSmiles("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O");
    };

    // Fetch the SVG representation of the input SMILES string.
    const drawSmiles = async () => {
        if (smiles === "") { return; }

        const data = { 
            "smiles": smiles,
            "height": 240,
            "width": 240
        };

        try {
            const response = await fetch("/api/smiles_to_svg", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ "data": data })
            });

            if (!response.ok) { 
                throw new Error("Network response was not ok!"); 
            };

            const json = await response.json();
            if (json.status === "success") { 
                setInputSvg(json.payload.svg); 
            };
        } catch (error) {
            toast.error(error.message);
        };
    };

    // When the SMILES string changes, fetch the SVG representation.
    useEffect(() => {
        drawSmiles();
    }, [smiles]);

    return (
        <div 
            className="control" 
            style={{
                border: "1px solid #dbdbdb", 
                borderRadius: "5px", 
                marginBottom: "10px"
            }}
        >
            <div className="panel">
                <div className="panel-heading">
                    <div className="title is-5">
                        Input
                    </div>
                </div>
                <div className="panel-block">
                    <div 
                        className="field" 
                        style={{
                            width: "100%", 
                            margin: "10px"
                        }}
                    >
                        <div className="control">
                            <input 
                                className="input" 
                                value={smiles} 
                                type="text" 
                                placeholder="Enter SMILES" 
                                onChange={(e) => setSmiles(e.target.value)} 
                            />
                        </div>
                    </div>
                </div>
                <div 
                    className="panel-block" 
                    style={{ 
                        padding: "10px", 
                        height: "250px",
                        display: "flex",
                        justifyContent: "left",
                    }}
                >
                    <div dangerouslySetInnerHTML={{ __html: inputSvg }} />
                </div>
                <div className="panel-block">
                    <div 
                        className="field has-addons" 
                        style={{ margin: "10px" }}
                    >
                        <div className="control">
                            <button 
                                className="button is-link"
                                style={{
                                    marginRight: "5px", 
                                    marginBottom: "5px"
                                }} 
                                onClick={parseMolecule} 
                            >
                                Parse molecule
                            </button>
                            <button 
                                className="button is-link is-light" 
                                style={{ 
                                    marginRight: "5px", 
                                    marginBottom: "5px" 
                                }} 
                                onClick={loadExample} 
                            >
                                Load example
                            </button>
                            <button 
                                className="button is-link is-light" 
                                style={{
                                    marginRight: "5px", 
                                    marginBottom: "5px"
                                }} 
                                onClick={() => {
                                    setInputSvg("");
                                    handleRefresh();
                                }} 
                            >
                                Refresh
                            </button>
                            <button 
                                className="button is-link is-light"
                                style={{
                                    marginRight: "5px", 
                                    marginBottom: "5px"
                                }} 
                                onClick={() => {
                                    if (smiles === "") { 
                                        toast.warn("No SMILES input to submit.");
                                        return; 
                                    }
                                    toggleReviewModal();
                                }} 
                            >
                                Submit for review
                            </button>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default SmilesInput;