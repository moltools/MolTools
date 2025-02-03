import React, { useState, } from "react";
import { toast } from "react-toastify";
import Plot from "react-plotly.js";
import SmilesDrawerContainer from "../components/common/SmilesDrawer";
import { Link } from "react-router-dom";
import { 
    useMediaQuery,
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
    TextField
} from "@mui/material";
import { 
    BugReport as BugReportIcon,
    Home as HomeIcon, 
    Menu as MenuIcon,
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

const Biosynfoni = () => {
    const [version, setVersion] = useState("0.0.0");
    const [smiles, setSmiles] = useState("");
    const [predictions, setPredictions] = useState({});
    const [isLoading, setIsLoading] = useState(false);

    const [isDrawerOpen, setIsDrawerOpen] = useState(false);
    const isNarrowScreen = useMediaQuery("(max-width: 600px)");


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
                                {`Biosynfoni (${version})`}
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
                <Box sx={{ display: "flex", flexDirection: "row", alignItems: "center", justifyContent: "center", gap: 2, pb: 10 }}>
                    <TextField 
                        disabled={isLoading}
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
                </Box>
                <Box sx={{ 
                    display: "flex", 
                    flexDirection: isNarrowScreen ? "column" : "row", 
                    alignItems: "center", 
                    justifyContent: "center", 
                    gap: 10 
                }}>
                    <Box>
                        <SmilesDrawerContainer 
                            identifier={"input-molecule"}
                            smilesStr={smiles} 
                            width={300}
                            height={350}
                        />
                    </Box>
                    <Box>
                        <PredictionView predictions={predictions} />
                    </Box>
                </Box>
            </Box>
        </Box>
    );
};

export default Biosynfoni;