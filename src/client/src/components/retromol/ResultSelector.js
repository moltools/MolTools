import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";

const ResultSelector = ({ 
    results, 
    selectedResult, 
    setSelectedResult,
    linearizedSmiles
}) => {
    const [linearizedSvg, setLinearizedSvg] = useState("");

    // Fetch the SVG representation of the linearized SMILES string.
    const drawSmiles = async () => {
        if (linearizedSmiles === "") { return; }

        const data = { 
            "smiles": linearizedSmiles,
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
                setLinearizedSvg(json.payload.svg); 
            };
        } catch (error) {
            toast.error(error.message);
        };
    };

    // When the linearized SMILES string changes, fetch the SVG representation.
    useEffect(() => {
        drawSmiles();
    }, [linearizedSmiles]);
    
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
                        Results
                    </div>
                </div>
                {linearizedSmiles && (
                    <div 
                    className="panel-block" 
                    style={{ 
                        padding: "10px", 
                        height: "250px",
                        display: "flex",
                        justifyContent: "left",
                        }}
                    >
                        <div dangerouslySetInnerHTML={{ __html: linearizedSvg }} />
                    </div>
                )}
                <div className="panel-block">
                    <div>
                        {results.length === 0 ? (
                            <p style={{ margin: "10px" }}>
                                No results found.
                            </p>
                        ) : (
                            <div 
                                className="field is-grouped is-grouped-multiline" 
                                style={{ 
                                    width: "100%", 
                                    margin: "10px" 
                                }}
                            >
                                {results.map((result, index) => (
                                    <div key={index} className="control">
                                        <div
                                            className={`tags has-addons ${selectedResult === result ? "selected" : ""}`}
                                            style={{ marginRight: "5px", cursor: "pointer" }}
                                            onClick={() => setSelectedResult(result)}
                                        >
                                            <span className={`tag is-link ${selectedResult === result ? "" : "is-light"}`}>
                                                Result {index + 1}
                                            </span>
                                        </div>
                                    </div>
                                ))}
                            </div>
                        )}
                    </div>
                </div>
            </div>
        </div>
    );
};

export default ResultSelector;