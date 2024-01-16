import React, { useState } from "react";
import { toast } from "react-toastify";
import TextInput from "../components/TextInput";

// =====================================================================================================================
// PIKAChU component.
// =====================================================================================================================

/**
 * PICAChU Component
 *
 * This component represents the "PIKAChU" section of the MolTools application. 
 *
 * @returns {JSX.Element} The rendered PIKAChU component.
 */
const PIKAChU = () => {
    const [smiles, setSmiles] = useState("");
    const [svgString, setSvgString] = useState("");

    // Downloads the SVG representation of the molecular model.
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

    // Every time smiles is updated, send smiles to backend and receive svgString of drawn molecule.
    React.useEffect(() => {
        // Draw SMILES.
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

    return (
        <div className="pikachu">
            <div className="biosynfoni-input-container">
                {/* TextInput component */}
                <TextInput 
                    className="pikachu-smiles-input"
                    locked={false} 
                    active={false} 
                    value="" 
                    error="" 
                    label="Enter SMILES" 
                    setValue={setSmiles}
                />
                {/* Button to submit the SMILES */}
                <button 
                    className="pikachu-submit-button"
                    onClick={() => {
                        if (smiles) {
                            handleDownloadSvgString();
                        } else {
                            toast.error("First enter a SMILES string!");
                        }
                    }}  
                >
                    Download
                </button>
            </div>
            <div className="pikachu-container">
                <div className="pikachu-mol-view">
                    {/* Render the SVG string as HTML */}
                    <div dangerouslySetInnerHTML={{ __html: svgString }} />
                </div>
            </div>
            <div className="pikachu-link-container">
                {/* Doi with link to paper */}
                <a 
                    className="pikachu-paper-link"
                    href="https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00616-5"
                    target="_blank"
                    rel="noopener noreferrer"
                >
                    10.1186/s13321-022-00616-5
                </a>
                {/* Link with doi to paper */}
                <a 
                    className="pikachu-github-link"
                    href="https://github.com/BTheDragonMaster/pikachu"
                    target="_blank"
                    rel="noopener noreferrer"
                >
                    github.com/BTheDragonMaster/pikachu
                </a>
            </div>
        </div>
    );
};

// Export the PICAChU component as the default export.
export default PIKAChU;