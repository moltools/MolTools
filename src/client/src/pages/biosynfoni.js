import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";
import Plot from "react-plotly.js";

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

const Overview = () => {
    return (
        <div>
            <h1 class="title">Biosynfoni</h1>
            <p>
                Biosynfoni is a tool for the prediction of the biosynthetic class of a given molecule. The scope of the tool
                is to provide a user-friendly interface for the prediction of the biosynthetic class of a molecule, based on
                its SMILES string representation. The tool is able to predict the biosynthetic class of a molecule using a
                pre-trained machine learning model.
            </p>
        </div>
    );
};

const PredictMolecule = () => {
    const [smiles, setSmiles] = useState("");
    const [svgString, setSvgString] = useState("");
    const [predictions, setPredictions] = useState({});
    const [isLoading, setIsLoading] = useState(false);

    const getPredictions = async () => {
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

        setIsLoading(false);
    };

    const clearPredictions = () => {
        setSmiles("");
        setPredictions({});
        setSvgString("");
    };

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

    const loadExample = () => {
        setSmiles("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O");
    };

    return (
        <div class="column is-full">
            <div class="field">
                <div 
                    class="control"
                >
                    <input class="input" value={smiles} type="text" placeholder="Enter SMILES" onChange={(e) => setSmiles(e.target.value)} />
                </div>
            </div>
            <div class="control">
                <button class="button is-link is-light" style={{marginRight: "5px"}} onClick={loadExample}>Example</button>
                <button class="button is-link is-light" style={{marginRight: "5px"}} onClick={getPredictions}>Submit</button>
                <button class="button is-link is-light" onClick={clearPredictions}>Clear</button>
            </div>
            <div class="columns" style={{marginTop: "5px"}}>
                <div class="column has-text-centered" style={{margin: "5px"}}>
                    <div dangerouslySetInnerHTML={{ __html: svgString }} />
                </div>
                <div class="column has-text-centered" style={{margin: "5px"}}>
                    <PredictionView predictions={predictions} />
                </div>
            </div>
        </div>
    );
};     

const Biosynfoni = () => {
    const tabs = ["Overview", "Predict biosynthetic class"];
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
                {selectedTab === "Predict biosynthetic class" && <PredictMolecule />}
            </div>
        </div>
    );
};

export default Biosynfoni;