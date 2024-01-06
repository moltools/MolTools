import React, { useState } from "react";
import { toast } from "react-toastify";
import Plot from "react-plotly.js";

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
        <div className={fieldClassName}>
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
// PredictionView component.
// =====================================================================================================================

/**
 * PredictionView component displays a bar chart of predictions.
 * 
 * This component uses the Plotly library to render a bar chart of predictions.
 * 
 * @param {Object} props - The props for the PredictionView component.
 * @param {Object} props.predictions - A dictionary with class names as keys and probabilities as values.
 * @returns {JSX.Element} - The rendered PredictionView component.
 */
const PredictionView = (props) => {
    // Parse predictions into Plotly format. Only keep the top 5 predictions and sort from highest to lowest probability.
    const predictions = props.predictions;
    const sortedPredictions = Object.keys(predictions).sort((a, b) => predictions[b] - predictions[a]);
    const topPredictions = sortedPredictions.slice(0, 5).reverse();
    const topProbabilities = topPredictions.map((prediction) => predictions[prediction]);

    // Format keys to have a maximum of 10 characters. If longer, truncate and add "...".
    for (let i = 0; i < topPredictions.length; i++) {
        if (topPredictions[i].length > 10) {
            topPredictions[i] = topPredictions[i].slice(0, 10) + "...";
        }
    }

    // Plotly data.
    const data = [
        {
            x: topProbabilities,
            y: topPredictions,
            type: "bar",
            orientation: "h",
            marker: { color: "#B9C311" }
        }
    ];

    // Plotly layout.
    const layout = {
        width: 500,
        height: 350,
        xaxis: {
            title: "Probability",
            range: [0, 1.1],
            tickformat: ",.0",
            automargin: true,
            titlefont: { size: 14, family: "HelveticaNeue-Light" },
            tickfont: { size: 12, family: "HelveticaNeue-Light" }
        },
        yaxis: {
            range: [-1, 5],
            automargin: true,
            titlefont: { size: 14, family: "HelveticaNeue-Light" },
            tickfont: { size: 12, family: "HelveticaNeue-Light" }
        },
        margin: { l: 0, r: 0, b: 0, t: 0, pad: 0 }
    };

    // Only return the plot if predictions are available. Check this by verifying if props.predictions has any keys.
    if (Object.keys(props.predictions).length > 0) {
        return (
            <Plot
                data={data}
                layout={layout}
                config={{ responsive: true }}
            />
        );
    } else {
        return <div />;
    };
};

// =====================================================================================================================
// Biosynfoni component.
// =====================================================================================================================

/**
 * Biosynfoni component widget
 * 
 * This component represents a widget for inputting SMILES (Simplified Molecular Input Line Entry System)
 * strings, sending them to a backend API for molecule drawing and biosynthetic class predictions,
 * and displaying the results.
 * 
 * State Variables:
 * - smiles: The SMILES string entered by the user.
 * - svgString: The SVG representation of the drawn molecule.
 * - predictions: Predictions related to the molecule's biosynthetic class.
 * 
 * Side Effects:
 * - Displays a warning toast message when the page loads.
 * - Calls the 'drawSmiles' function to fetch and update 'svgString' when 'smiles' changes.
 * - Calls the 'getPredictions' function to fetch and update 'predictions' when 'smiles' changes.
 * - Displays toast messages for success, error, or warning messages from the backend.
 * 
 * Usage:
 * - Include this component in your React application to create a molecule drawing and prediction widget.
 * - Requires the 'SmilesInput' and 'PredictionView' components.
 * 
 * @returns {JSX.Element} The rendered Biosynfoni component.
 */
const Biosynfoni = () => {
    const [smiles, setSmiles] = useState("");
    const [svgString, setSvgString] = useState("");
    const [predictions, setPredictions] = useState({});

    // Send toast on page load.
    React.useEffect(() => {
        toast.warn("This is a demo page!", { autoClose: false });
    }, []);

    // Send smiles to backend and receive svgString of drawn molecule.
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

    // Send smiles to backend and receive predictions.
    const getPredictions = async () => {
        try {
            const response = await fetch("/api/predict_biosynthetic_class", {
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
                setPredictions(json.payload.predictions);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
       
        } catch (error) {
            const msg = "Could not get predictions!";
            toast.error(msg, { autoClose: true });
            console.error(error);
        };
    };

    // Every time smiles is updated, send smiles to backend and receive svgString of drawn molecule.
    React.useEffect(() => {
        if (smiles) {
            drawSmiles();
        };
    }, [smiles]);

    // Render the Biosynfoni widget.
    return (
        <div className="biosynfoni">
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
                    className="biosynfoni-submit-button"
                    onClick={() => {
                        if (smiles) {
                            getPredictions();
                        } else {
                            toast.error("First enter a SMILES string!");
                        }
                    }}  
                >
                    Submit
                </button>
            </div>
            <div className="biosynfoni-result-container">
                <div className="biosynfoni-mol-view">
                    {/* Render the SVG string as HTML */}
                    <div dangerouslySetInnerHTML={{ __html: svgString }} />
                </div>
                <div className="biosynfoni-pred-view">
                    {/* Render the PredictionView component with 'predictions' */}
                    <PredictionView predictions={predictions} />
                </div>
            </div>
        </div>
    );
};

// Export the Biosynfoni component as the default export.
export default Biosynfoni;