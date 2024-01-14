import React, { useState } from "react";
import { toast } from "react-toastify";
import Plot from "react-plotly.js";
import TextInput from "../components/TextInput";

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
        if (topPredictions[i].length > 20) {
            topPredictions[i] = topPredictions[i].slice(0, 20) + "...";
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
            title: {
                text: "Probability",
                standoff: 20,
                font: { size: 14, family: "HelveticaNeue-Light" }
            },
            range: [0, 1.1],
            tickformat: ",.1",
            hoverformat: ",.3",
            automargin: true,
            titlefont: { size: 14, family: "HelveticaNeue-Light" },
            tickfont: { size: 12, family: "HelveticaNeue-Light" },
            ticklen: 10
        },
        yaxis: {
            range: [-1, 5],
            automargin: true,
            titlefont: { size: 14, family: "HelveticaNeue-Light" },
            tickfont: { size: 12, family: "HelveticaNeue-Light" },
            ticklen: 10,
        },
        margin: { l: 100, r: 0, b: 0, t: 0, pad: 0 }
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
 * - Calls the 'drawSmiles' function to fetch and update 'svgString' when 'smiles' changes.
 * - Calls the 'getPredictions' function to fetch and update 'predictions' when 'smiles' changes.
 * - Displays toast messages for success, error, or warning messages from the backend.
 * 
 * Usage:
 * - Include this component in your React application to create a molecule drawing and prediction widget.
 * - Requires the 'TextInput' and 'PredictionView' components.
 * 
 * @returns {JSX.Element} The rendered Biosynfoni component.
 */
const Biosynfoni = () => {
    const [smiles, setSmiles] = useState("");
    const [svgString, setSvgString] = useState("");
    const [predictions, setPredictions] = useState({});
    const [isLoading, setIsLoading] = useState(false);

    // Send smiles to backend and receive predictions.
    const getPredictions = async () => {
        
        // Set isLoading to true, which disables submit button.
        setIsLoading(true);

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
                setSvgString(json.payload.svg_string);
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

        // Set isLoading to false, which enables submit button.
        setIsLoading(false);
    };

    // Every time smiles is updated, send smiles to backend and receive svgString of drawn molecule.
    React.useEffect(() => {
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

        if (smiles) {
            drawSmiles();
        };
    }, [smiles]);

    // Render the Biosynfoni widget.
    return (
        <div className="biosynfoni">
            <div className="biosynfoni-input-container">
                {/* TextInput component */}
                <TextInput 
                    className="biosynfoni-smiles-input"
                    locked={isLoading} 
                    active={false} 
                    value="" 
                    error="" 
                    label="Enter SMILES" 
                    setValue={setSmiles}
                />
                {/* Button to submit the SMILES */}
                <button 
                    disabled={isLoading}
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