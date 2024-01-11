import React, { useState } from "react";
import { toast } from "react-toastify";

// =====================================================================================================================
// SmilesInput component.
// =====================================================================================================================

/**
 * SmilesInput component represents an input field with dynamic behavior based on props.
 * 
 * This component represents an input field with dynamic behavior based on props. It is used in the Biosynfoni component
 * to allow users to enter SMILES strings.
 * 
 * @param {Object} props - The props for the SmilesInput component.
 * @param {boolean} props.locked - Determines if the input field is locked.
 * @param {boolean} props.active - Indicates if the input field is active.
 * @param {string} props.value - The current value of the input field.
 * @param {string} props.error - Any error message associated with the input.
 * @param {string} props.label - The label to display as a placeholder.
 * @param {function} props.setValue - A function to handle value changes.
 * @param {Array} props.predicted - An array of predicted values for the input.
 * @returns {JSX.Element} - The rendered SmilesInput component.
 */
const SmilesInput = (props) => {
    // Destructure props into individual variables.
    const { 
        locked,
        active: propsActive, 
        value: propsValue, 
        error: propsError, 
        label: propsLabel, 
        setValue, 
        predicted 
    } = props;

    // Define state variables using the useState hook
    const [active, setActive] = useState(locked ? propsActive : false);
    const [value, setValueState] = useState(propsValue || "");
    const [error, setError] = useState(propsError || "");
    const [label, setLabel] = useState(propsLabel || "Label");

    /**
     * Handles changes in the input value.
     * @param {Object} event - The input change event.
     */
    const changeValue = (event) => {
        const newValue = event.target.value;
        setValue(newValue);
        setValueState(newValue);
        setError("");
    };

    // Determine the CSS class for the input field based on conditions
    const fieldClassName = `field ${(locked ? active : active || value) && "active"} ${locked && !active && "locked"}`;

    return (
        <div className={`${fieldClassName} pikachu-smiles-input`}>
            {active && value && predicted && predicted.includes(value) && <p className="predicted">{predicted}</p>}
            <input
                id={1}
                type="text"
                value={value}
                placeholder={label}
                onChange={changeValue}
                onFocus={() => !locked && setActive(true)}
                onBlur={() => !locked && setActive(false)}
            />
            <label htmlFor={1} className={error && "error"}>
                {error || label}
            </label>
        </div>
    );
};

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
        if (smiles) {
            drawSmiles();
        };
    }, [smiles]);

    return (
        <div className="widget-content pikachu">
            <div className="biosynfoni-input-container">
                {/* SmilesInput component */}
                <SmilesInput 
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
        </div>
    );
};

// Export the PICAChU component as the default export.
export default PIKAChU;