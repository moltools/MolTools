import React, { useState, useEffect } from "react";
import ReactDOM from "react-dom/client";
import { BrowserRouter, Routes, Route, useLocation } from "react-router-dom";
import { FiHexagon } from "react-icons/fi";
import { toast } from "react-toastify";

// Import styling.
import "./style/main.css";

// Import widgets.
import Header from "./components/Header";
import Toast from "./components/Toast";

// Import pages.
import Home from "./pages";
import About from "./pages/about";
import Changelog from "./pages/changelog";
import RetroMol from "./pages/retromol";
import CineMol from "./pages/cinemol";
import Biosynfoni from "./pages/biosynfoni";
import MechaMol from "./pages/mechamol";
import Comet from "./pages/comet";

/**
 * Define widgets which will be displayed as tiles on the home page.
 * 
 * @param {string} name - Name of the widget.
 * @param {string} icon - Icon of the widget.
 * @param {string} color - Color of the widget.
 * @param {string} description - Description of the widget.
 * @param {string} path - Path to the widget.
 * 
 * @returns {object} - Object containing the widget data.
 */
const tilesData = [
    { name: "CineMol", icon: "widgets/icon_cinemol.svg", color: "#2f6eb5", description: "Alpha version", path: "/cinemol" },
    { name: "Biosynfoni", icon: "widgets/icon_biosynfoni.svg", color: "#7B9204", description: "Demo version", path: "/biosynfoni" },
    { name: "RetroMol", icon: "widgets/icon_retromol.svg", color: "#F5C900", description: "Under construction", path: "/retromol" },
    { name: "MechaMol", icon: <FiHexagon />, color: "#ccc", description: "Under construction", path: "/mechamol" },
    { name: "Comet", icon: <FiHexagon />, color: "#ccc", description: "Under construction", path: "/comet" }
];

/**
 * Define the display names of the routes.
 * 
 * @param {string} path - Path to the route.
 * @param {string} name - Display name of the route.
 * 
 * @returns {object} - Object containing the route display names.
 */
const routeDisplayNames = {
    "/": "MolTools",
    "/about": "About",
    "/changelog": "Changelog",
    "/retromol": "RetroMol",
    "/cinemol": "CineMol",
    "/biosynfoni": "Biosynfoni",
    "/mechamol": "MechaMol",
    "/comet": "Comet",
};

/**
 * Define the routes of the app.
 */
function AppRoutes () {
    const location = useLocation();
    const [displayLocation, setDisplayLocation] = useState(location);
    const [transitionStage, setTransitionStage] = useState("fade-in");

    useEffect(() => {
        if (location !== displayLocation) setTransitionStage("fade-out");
    }, [location]);

    return (
        <div 
            className={`widget ${transitionStage}`}
            onTransitionEnd={() => {
                if (transitionStage === "fade-out") {
                    setTransitionStage("fade-in");
                    setDisplayLocation(location);
                };
            }}
        >
            <Routes location={displayLocation}>
                <Route path="/" element={<Home tilesData={tilesData}/>} />
                <Route path="/about" element={<About />} />
                <Route path="/changelog" element={<Changelog url={"https://raw.githubusercontent.com/moltools/MolTools/master/CHANGELOG.md"} />} />
                <Route path="/retromol" element={<RetroMol />} />
                <Route path="/cinemol" element={<CineMol />} />
                <Route path="/biosynfoni" element={<Biosynfoni />} />
                <Route path="/mechamol" element={<MechaMol />} />
                <Route path="/comet" element={<Comet />} />
                <Route path="/*" element={<div>404</div>} />
            </Routes>
        </div>
    );
};

function App () {
    const [apiOk, setApiOk] = useState(false);

    /**
     * Check if the API is available.
     */
    useEffect(() => {
        fetch("/api/ping_server")
            .then(response => response.json())
            .then(data => { if (data.status === "success") { setApiOk(true); }; })
            .catch(error => { console.log(error); });
    }, []);

    /**
     * Render the app.
     */
    return (
        <div className="app">
            <BrowserRouter>
                <Header 
                    routeDisplayNames={routeDisplayNames}
                    apiOk={apiOk}
                    navbarLinks={[
                        { name: "Home", path: "/" },
                        { name: "About", path: "/about" },
                        { name: "Changelog", path: "/changelog" }
                    ]}
                />
                <AppRoutes />
                <Toast />
            </BrowserRouter>    
        </div>
    );
};

// Render the app.
const root = ReactDOM.createRoot(document.getElementById("root"));

root.render(<App />);