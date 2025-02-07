import React, { useState, } from "react";
import { toast } from "react-toastify";
import Plot from "react-plotly.js";
import SmilesDrawerContainer from "../components/common/SmilesDrawer";
import { Link } from "react-router-dom";
import { 
    AppBar, 
    Box, 
    Button,
    Drawer, 
    IconButton, 
    List, 
    ListItem, 
    ListItemIcon, 
    ListItemText, 
    Toolbar, 
    Typography,
    CssBaseline,
    TextField,
    Tooltip
} from "@mui/material";
import { 
    BugReport as BugReportIcon,
    Home as HomeIcon, 
    Menu as MenuIcon,
    Info as InfoIcon,
    ArrowUpward as ArrowUpwardIcon,
    ArrowDownward as ArrowDownwardIcon,
} from "@mui/icons-material";


const SidebarButtonHome = () => {
    return (
        <ListItem 
            type="button" 
            component={Link} 
            to="/"
            sx={{ textDecoration: "none", color: "#222" }}
        >
            <ListItemIcon>
                <HomeIcon sx={{ color: "#222" }} />
            </ListItemIcon>
            <ListItemText primary="Home" />
        </ListItem>
    );
};


const SidebarButtonBugReport = () => {
    return (
        <ListItem 
            type="button" 
            component="a" 
            href="https://github.com/lucinamay/biosynfoni/issues" 
            target="_blank"
            sx={{ textDecoration: "none", color: "#222" }}
        >
            <ListItemIcon>
                <BugReportIcon sx={{ color: "#222" }} />
            </ListItemIcon>
            <ListItemText primary="Report Bug" />
        </ListItem>
    );
};

const getColor = (name) => {
    const nameToColor = {
        "coenzyme_a" : "#BFBFBF",
        "nadh" : "#BFBFBF",
        "nadph" : "#BFBFBF",
        "standard_amino_acids" : "#FFEAA0",
        "non-standard_amino_acids" : "#FFEAA0",
        "open_pyranose" : "#FFC4CE",  
        "open_furanose" : "#FFC4CE",  
        "pyranose" : "#FFC4CE",  
        "furanose" : "#FFC4CE",  
        "indoleC2N" : "#A783B6",
        "phenylC2N" : "#A783B6",
        "C5N" : "#FF8B61",
        "C4N" : "#FF8B61",
        "phenylC3" : "#A783B6",
        "phenylC2" : "#A783B6",
        "phenylC1" : "#A783B6",
        "isoprene" : "#B9C311",  
        "acetyl" : "#FF8B61",
        "methylmalonyl" : "#FF8B61",
        "ethyl" : "#FF8B61",
        "methyl" : "#FF8B61",
        "phosphate" : "#BFBFBF",
        "sulfonate" : "#BFBFBF",
        "fluorine" : "#BFBFBF",
        "chlorine" : "#BFBFBF",
        "bromine" : "#BFBFBF",
        "iodine" : "#BFBFBF",
        "nitrate" : "#BFBFBF",
        "epoxy" : "#BFBFBF",
        "ether" : "#BFBFBF",
        "hydroxyl" : "#BFBFBF",
        "C3_ring" : "#BFBFBF",
        "C4_ring" : "#BFBFBF",
        "C5_ring" : "#BFBFBF",
        "C6_ring" : "#BFBFBF",
        "C7_ring" : "#BFBFBF",
        "C8_ring" : "#BFBFBF",
        "C9_ring" : "#BFBFBF",
        "C10_ring" : "#BFBFBF",
    };
    const defaultColor = "#ceccca";
    return nameToColor[name] ?? defaultColor;
};

const getLabel = (name) => {
    const nameToLabel = {
        "coenzyme_a" : "Coenzyme A",
        "nadh" : "NADH",
        "nadph" : "NADPH",
        "standard_amino_acids" : "Standard amino acids",
        "non-standard_amino_acids" : "Non-standard amino acids",
        "open_pyranose" : "Open pyranose",
        "open_furanose" : "Open furanose",
        "pyranose" : "Pyranose",
        "furanose" : "Furanose",
        "indoleC2N" : "Indole C2N",
        "phenylC2N" : "Phenyl C2N",
        "C5N" : "C5N",
        "C4N" : "C4N",
        "phenylC3" : "Phenyl C3",
        "phenylC2" : "Phenyl C2",
        "phenylC1" : "Phenyl C1",
        "isoprene" : "Isoprene",
        "acetyl" : "Acetyl",
        "methylmalonyl" : "Methylmalonyl",
        "ethyl" : "Ethyl",
        "methyl" : "Methyl",
        "phosphate" : "Phosphate",
        "sulfonate" : "Sulfonate",
        "fluorine" : "Fluorine",
        "chlorine" : "Chlorine",
        "bromine" : "Bromine",
        "iodine" : "Iodine",
        "nitrate" : "Nitrate",
        "epoxy" : "Epoxy",
        "ether" : "Ether",
        "hydroxyl" : "Hydroxyl",
        "C3_ring" : "C3 ring",
        "C4_ring" : "C4 ring",
        "C5_ring" : "C5 ring",
        "C6_ring" : "C6 ring",
        "C7_ring" : "C7 ring",
        "C8_ring" : "C8 ring",
        "C9_ring" : "C9 ring",
        "C10_ring" : "C10 ring",
    };
    return nameToLabel[name] ?? name;
};

const FingerprintBar = ({ data, selected, setSelected, topPos = [], topNeg = [] }) => {

    return (
        <Box 
            display="flex" 
            sx={{ 
                overflowX: "auto", 
                width: "100%", 
                center: "center", 
                flexWrap: "wrap", 
                p: 5, 
                justifyContent: "center", 
                alignItems: "center",
                backgroundColor: "#ececec",
                borderRadius: "10px",
                boxShadow: "0px 0px 5px 0px rgba(0,0,0.5,0)",
                borderWidth: "5px",
                borderColor: "#ceccca", 
            }}
        >
            {data.map(([name, count], index) => (
                <Tooltip 
                    key={index} 
                    title={
                        <>
                            Feature: {getLabel(name)} <br />
                            Count: {count} <br />
                            {count >= 1 ? "Click to highlight substructure(s) in the molecule." : "No substructures found."}
                        </>
                    }
                    arrow
                >
                    <Box sx={{ 
                        display: "flex", 
                        flexDirection: "column", 
                        alignItems: "center", 
                        justifyContent: "center", 
                        p: "1.5px", 
                        borderWidth: "0px",
                        // backgroundColor: "white",
                    }}>
                        <Box
                            sx={{
                                width: `40px`,
                                height: "40px",
                                backgroundColor: getColor(name),
                                borderRadius: "3px",
                                display: "flex",
                                alignItems: "center",
                                justifyContent: "center",
                                color: selected === name ? "#fff" : "#000",
                                fontSize: "12px",
                                fontFamily: "HelveticaNeue-Light",
                                cursor: count >= 1 ? "pointer" : "default",
                                boxShadow: "rgba(50, 50, 93, 0.25) 0px 50px 100px -20px, rgba(0, 0, 0, 0.3) 0px 30px 60px -30px, rgba(10, 37, 64, 0.35) 0px -2px 6px 0px inset",
                            }}
                            // if double click unselect
                            onClick={() => {
                                // if count < 1 do nothing
                                if (count < 1) {
                                    return;
                                }

                                if (selected === name) {
                                    setSelected("");
                                } else {
                                    setSelected(name);
                                }
                            }}
                        >
                            <Box
                                sx={{
                                    display: "flex",
                                    flexDirection: "column",
                                    alignItems: "center",
                                    justifyContent: "center",
                                    fontFamily: "HelveticaNeue-Light",
                                    fontSize: "12px",
                                    color: selected === name ? "#fff" : "#000",
                                }}
                            >
                                {count}
                                <Box sx={{ display: "flex", alignItems: "center", justifyContent: "center" }}>
                                    {topPos.includes(name) && <ArrowUpwardIcon fontSize="xsmall" />}
                                    {topNeg.includes(name) && <ArrowDownwardIcon fontSize="xsmall" />}
                                </Box>
                            </Box>
                        </Box>
                    </Box>
                </Tooltip>
            ))}
        </Box>
    );
};


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

    // Helper function: assign a color based on the label.
    const getColorForLabel = (label) => {
        if (label.toLowerCase() === "alkaloid") return "#B4CAD8";
        else if (label.toLowerCase() === "amino acid") return "#FFEAA0";
        else if (label.toLowerCase() === "carbohydrate") return "#FFC4CE";
        else if (label.toLowerCase() === "fatty acid") return "#F6B26B";
        else if (label.toLowerCase() === "isoprenoid") return "#B9C311";
        else if (label.toLowerCase() === "phenylpropanoid") return "#A783B6";
        else if (label.toLowerCase() === "polyketide") return "#FF8B61";
        else return "#B9C311"; // Default color.
    };

    // Create an array of colors corresponding to each bar using the helper function.
    const barColors = topPredictions.map(label => getColorForLabel(label));

    // Plotly data.
    const data = [
        {
            x: topProbabilities,
            y: topPredictions,
            type: "bar",
            orientation: "h",
            marker: { color: barColors }
        }
    ];

    // Plotly layout.
    const layout = {
        // title: {text: 'Biosynthetic class predictions', font: { size: 16, family: "HelveticaNeue-Light" }},
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
        margin: { l: 100, r: 0, b: 0, t: 50, pad: 0 }
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

const Biosynfoni = () => {
    const [version, setVersion] = useState("0.0.0");
    const [smiles, setSmiles] = useState("");
    const [isLoading, setIsLoading] = useState(false);
    const [isDrawerOpen, setIsDrawerOpen] = useState(false);

    const [predictions, setPredictions] = useState({});
    const [taggedSmiles, setTaggedSmiles] = useState("");
    const [featureImportance, setFeatureImportance] = useState({});
    const [fingerprint, setFingerprint] = useState([]);
    const [highlights, setHighlights] = useState({});
    const [topPredictedClass, setTopPredictedClass] = useState("");
    const [selectedHighlight, setSelectedHighlight] = useState("");

    // handle for clearing SMILES and predictions.
    const clearData = () => {
        setSmiles("");
        setPredictions({});
        setTaggedSmiles("");
        setFeatureImportance({});
        setFingerprint([]);
        setHighlights({});
        setTopPredictedClass("");
        setSelectedHighlight([]);
    };

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
                setTaggedSmiles(json.payload.tagged_smiles);
                setFeatureImportance(json.payload.feature_importance);
                setFingerprint(json.payload.fingerprint);
                setHighlights(json.payload.highlights);

                // get top predicted class from predictions
                const topClass = Object.keys(json.payload.predictions).reduce((a, b) => json.payload.predictions[a] > json.payload.predictions[b] ? a : b);
                setTopPredictedClass(topClass);

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

    const getVersion = async () => {
        try {
            const response = await fetch("/api/biosynfoni_version");
            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };
            const json = await response.json();
            setVersion(json.payload.version);
        } catch (error) {
            const msg = "Could not get version!";
            toast.error(msg, { autoClose: true });
            console.error(error);
        };
    }

    // Fetch the version number when the component mounts.
    React.useEffect(() => {
        getVersion();
    }, []);

    const toggleDrawer = () => {
        setIsDrawerOpen(!isDrawerOpen);
    };

    // Render the Biosynfoni widget.
    return (
        <Box sx={{ display: "flex", width: "100%", flexDirection: "column" }}>
            <CssBaseline />
            <AppBar position="fixed">
                <Toolbar sx={{ display: "flex", justifyContent: "space-between", backgroundColor: "#7B9204" }}>
                    <Box sx={{ display: "flex", justifyContent: "left", width: "100%", alignItems: "center" }}>
                        <IconButton
                            color="inherit"
                            aria-label="open drawer"
                            edge="start"
                            onClick={toggleDrawer}
                            sx={{ mr: 2 }}
                        >
                            <MenuIcon />
                        </IconButton>
                        <Box sx={{ display: "flex", alignItems: "center" }}>
                            <Box sx={{ display: "flex", alignItems: "center", mr: 2 }}>
                                <img src="/biosynfoni.svg" alt="Biosynfoni" style={{ width: 40, height: 40 }} />
                            </Box>
                            <Typography variant="h6" noWrap>
                                {`Biosynfoni v${version} (web app: v${process.env.REACT_APP_VERSION ? process.env.REACT_APP_VERSION : 'UNKNOWN'})`}
                            </Typography>
                        </Box>
                    </Box>
                </Toolbar>
            <Drawer
                variant="temporary"
                open={isDrawerOpen}
                onClose={toggleDrawer}
                ModalProps={{ keepMounted: true }}
                sx={{ '& .MuiDrawer-paper': { width: 240 } }}
            >
                <List>
                    <SidebarButtonHome />
                    <SidebarButtonBugReport />
                </List>
            </Drawer>
            </AppBar>
            <Box component="main" sx={{ flexGrow: 1, p: 3, pt: 5, height: "100vh", overflow: "auto" }}>
                <Toolbar />
                <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", justifyContent: "center", gap: 1, pb: 5 }}>
                    <TextField 
                        // disabled if loading or if there are predictions
                        disabled={isLoading || Object.keys(predictions).length > 0}
                        label="Enter SMILES"
                        variant="outlined"
                        fullWidth
                        value={smiles}
                        onChange={(event) => setSmiles(event.target.value)} 
                        sx={{ 
                            "& .MuiInputBase-input": { 
                                fontFamily: "HelveticaNeue-Light" 
                            },
                            "& .MuiOutlinedInput-root": {
                                "&:hover .MuiOutlinedInput-notchedOutline": {
                                    borderColor: "#B9C311", // Change outline color on hover
                                },
                                "&.Mui-focused .MuiOutlinedInput-notchedOutline": {
                                    borderColor: "#9FA80E", // Change outline color when focused
                                    borderWidth: "2px", // Optional: make it thicker when focused
                                }
                            },
                            "& .MuiInputLabel-root": {
                                fontFamily: "HelveticaNeue-Light", 
                                color: "#666", // Default label color
                                "&.Mui-focused": {
                                    color: "#7B9204", // Label color when focused
                                }
                            }
                        }}
                    />
                    <Button 
                        sx={{ 
                            height: 56, 
                            width: 100, 
                            fontFamily: "HelveticaNeue-Light", 
                            backgroundColor: "#B9C311", 
                            color: "#222", 
                            "&:hover": {
                                backgroundColor: "#7B9204", // Change to a slightly darker shade on hover
                                color: "#fff" // Change text color on hover
                            }
                        }}
                        disabled={isLoading}
                        onClick={() => {if (smiles) { getPredictions(); } else { toast.error("First enter a SMILES string!");}}}  
                        variant="contained"
                    >
                        Submit
                    </Button>
                    <Button 
                        sx={{ 
                            height: 56, 
                            width: 100, 
                            fontFamily: "HelveticaNeue-Light", 
                            backgroundColor: "#B9C311", 
                            color: "#222", 
                            "&:hover": {
                                backgroundColor: "#7B9204", // Change to a slightly darker shade on hover
                                color: "#fff" // Change text color on hover
                            }
                        }}
                        disabled={isLoading}
                        onClick={clearData}
                        variant="contained"
                    >
                        Clear
                    </Button>
                </Box>
                {fingerprint.length > 0 && (
                    <Box sx={{ display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "center", pb: 3 }}>
                        <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", justifyContent: "center", gap: 0.5 }}>
                            <Typography variant="h6" sx={{ fontFamily: "HelveticaNeue-Light" }}>
                                Fingerprint
                            </Typography>
                            <Tooltip title="Fingerprint that is used to make a prediction. Features with an upward arrow contributed positively to the top predicted class. Features with a downward arrow contributed negatively to the top predicted class. Click a feature to show its corresponding substructure(s) in the molecule." arrow>
                                <InfoIcon fontSize="small" />
                            </Tooltip>
                        </Box>
                        <FingerprintBar
                            data={fingerprint} 
                            selected={selectedHighlight}
                            setSelected={setSelectedHighlight}
                            topPos={featureImportance[topPredictedClass]?.top_pos ?? []}
                            topNeg={featureImportance[topPredictedClass]?.top_neg ?? []}
                        />
                    </Box>
                )}
                <Box
                    sx={{
                        display: "flex",
                        flexWrap: "wrap",
                        alignItems: "center",
                        justifyContent: "center",
                        width: "100%",
                        gap: 3.5,
                    }}
                >
                    <Box sx={{ display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "center", pb: 3 }}>
                        {smiles.length > 0 && (
                            <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", justifyContent: "center", gap: 0.5 }}>
                                <Typography variant="h6" sx={{ fontFamily: "HelveticaNeue-Light" }}>
                                    {fingerprint.length > 0 ? 'Submitted molecule' : 'Molecule for submission'}
                                </Typography>
                            </Box>
                        )}
                        <Box
                            sx={{
                                display: "flex",
                                flexDirection: "column",
                                alignItems: "center",
                                justifyContent: "center",

                                // same background as fingerprint
                                backgroundColor: smiles.length > 0 ? "#ececec" : "#fff",
                                borderRadius: "10px",
                                boxShadow: "0px 0px 5px 0px rgba(0,0,0.5,0)",
                                borderWidth: "5px",
                                borderColor: "#ceccca",
                                p: 3,
                            }}
                        >
                            {smiles.length > 0 && (
                                <SmilesDrawerContainer
                                    identifier={"input-molecule"}
                                    smilesStr={taggedSmiles.length > 0 ? taggedSmiles : smiles}
                                    width={350}
                                    height={350}
                                    highlightAtoms={highlights[selectedHighlight] ?? []}
                                    highlightColor={getColor(selectedHighlight)}
                                />
                            )}
                        </Box>
                    </Box>
                    <Box sx={{ display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "center", pb: 3 }}>
                        {fingerprint.length > 0 && (
                            <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", justifyContent: "center", gap: 0.5 }}>
                                <Typography variant="h6" sx={{ fontFamily: "HelveticaNeue-Light" }}>
                                    Biosynthetic class predictions
                                </Typography>
                            </Box>
                        )}
                        <Box 
                            sx={{ 
                                display: "flex", 
                                flexDirection: "column", 
                                alignItems: "center", 
                                justifyContent: "center",

                                // same background as fingerprint
                                backgroundColor: fingerprint.length > 0 ? "#ececec" : "#fff",
                                borderRadius: "10px",
                                boxShadow: "0px 0px 5px 0px rgba(0,0,0.5,0)",
                                borderWidth: "5px",
                                borderColor: "#ceccca",
                                p: 3,
                            }}
                        >
                            <PredictionView predictions={predictions} />
                        </Box>
                    </Box>
                </Box>

            </Box>
        </Box>
    );
};

export default Biosynfoni;