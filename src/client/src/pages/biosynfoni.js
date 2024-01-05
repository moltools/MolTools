import React, { useState } from "react";
import { toast } from "react-toastify";
import Plot from "react-plotly.js";

// =====================================================================================================================
// Component for SMILES input field.
//
// Source: https://medium.com/@whwrd/building-a-beautiful-text-input-component-in-react-f85564cc7e86
// =====================================================================================================================

class SmilesInput extends React.Component {
    constructor(props) {
        super(props);
    
        this.state = {
            active: (props.locked && props.active) || false,
            value: props.value || "",
            error: props.error || "",
            label: props.label || "Label"
        };
    };
  
    changeValue(event) {
        const value = event.target.value;
        this.props.setValue(value);
        this.setState({ value, error: "" });
    };
  
    handleKeyPress(event) {
        if (event.which === 13) {
            this.setState({ value: this.props.predicted });
        };
    };
  
    render() {
        const { active, value, error, label } = this.state;
        const { predicted, locked } = this.props;
        const fieldClassName = `field ${(locked ? active : active || value) && "active"} ${locked && !active && "locked"}`;
    
        return (
            <div className={fieldClassName}>
                {active &&
                    value &&
                    predicted &&
                    predicted.includes(value) && <p className="predicted">{predicted}</p>}
                <input
                    id={1}
                    type="text"
                    value={value}
                    placeholder={label}
                    onChange={this.changeValue.bind(this)}
                    onKeyPress={this.handleKeyPress.bind(this)}
                    onFocus={() => !locked && this.setState({ active: true })}
                    onBlur={() => !locked && this.setState({ active: false })}
                />
                <label htmlFor={1} className={error && "error"}>
                    {error || label}
                </label>
            </div>
        );
    };
};

// =====================================================================================================================
// Component for visualization of predictions.
// =====================================================================================================================

const PredictionView = (props) => {
    
    // Parse predictions into Plotly format. Only keep top 5 predictions and 
    // sort from highest to lowest probability. 
    // NOTE: props.predictions is a dictionary with keys as the class names and values as the probabilities.
    const predictions = props.predictions;
    const sortedPredictions = Object.keys(predictions).sort((b, a) => predictions[b] - predictions[a]);
    const topPredictions = sortedPredictions.slice(0, 5);
    const topProbabilities = topPredictions.map((prediction) => predictions[prediction]);

    // Format keys in such away that they have max 10 characters. If longer, truncate and add "...".
    for (let i = 0; i < topPredictions.length; i++) {
        if (topPredictions[i].length > 10) {
            topPredictions[i] = topPredictions[i].slice(0, 10) + "...";
        };
    };

    const data = [
        {
            x: topProbabilities,
            y: topPredictions,
            type: "bar",
            orientation: "h",
            marker: {color: "#B9C311"}
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
            titlefont: {size: 14, family: "HelveticaNeue-Light"}, 
            tickfont: {size: 12, family: "HelveticaNeue-Light"}
        },
        yaxis: {
            range: [-1, 5],
            automargin: true,
            titlefont: {size: 14, family: "HelveticaNeue-Light"},
            tickfont: {size: 12, family: "HelveticaNeue-Light"}
        },
        margin: {l: 0, r: 0, b: 0, t: 0, pad: 0}
    };

    // Only return plot if predictions are available. Check this seeing if props.predictions has any keys.
    if (Object.keys(props.predictions).length > 0) {
        return (
            <Plot
                data={data}
                layout={layout}
                config={ {responsive: true} }
            />  
        );
    } else {
        return <div/>;
    };
};

// =====================================================================================================================
// Render page.
// =====================================================================================================================

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
        const response = await fetch("/api/draw_smiles", {
            method: "POST",
            headers: {"Content-Type": "application/json"},
            body: JSON.stringify({ "smiles": smiles })
        });
        const body = await response.json();
        if (response.status !== 200) {
            throw Error(body.message);
        };
        setSvgString(body.payload.svg_string);
    };

    // Send smiles to backend and receive predictions.
    const getPredictions = async () => {
        const response = await fetch("/api/predict_biosynthetic_class", {
            method: "POST",
            headers: {"Content-Type": "application/json"},
            body: JSON.stringify({ "smiles": smiles })
        });
        const body = await response.json();
        if (response.status !== 200) {
            throw Error(body.message);
        };
        setPredictions(body.payload.predictions);
    };

    // Every time smiles is updated, send smiles to backend and receive svgString of drawn molecule.
    React.useEffect(() => {
        if (smiles) {
            drawSmiles();
        };
    }, [smiles]);

    // Render widget.
    return (
        <div className="biosynfoni">
            
            <div className="biosynfoni-input-container">
                <SmilesInput 
                    id={1}
                    label="SMILES input"
                    predicted=""
                    locked={false}
                    active={false}
                    setValue={setSmiles}
                />
                <button 
                    className="biosynfoni-submit-button"
                    onClick={() => {
                        if (smiles) {
                            getPredictions();
                            toast.success("SMILES submitted!");
                        } else {
                            toast.error("First enter a SMILES string!");
                        };
                    }}  
                >
                    Submit
                </button>
            </div>

            <div className="biosynfoni-result-container">
                <div className="biosynfoni-mol-view">
                    <div dangerouslySetInnerHTML={{ __html: svgString }} />
                </div>
                <div className="biosynfoni-pred-view">
                    <PredictionView predictions={predictions} />
                </div>
            </div>

        </div>
    );
};

export default Biosynfoni;