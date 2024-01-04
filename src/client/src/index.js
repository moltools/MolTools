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
import Footer from "./components/Footer";

// Import pages.
import Home from "./pages";
import About from "./pages/about";
import Changelog from "./pages/changelog";
import RetroMol from "./pages/retromol";
import CineMol from "./pages/cinemol";
import Biosynfoni from "./pages/biosynfoni";

// Define widgets which will be displayed as tiles on the home page.
const tilesData = [
    {
        name: "RetroMol",
        icon: <FiHexagon />,
        description: "",
        path: "/retromol"
    },
    {
        name: "CineMol",
        icon: <FiHexagon />,
        description: "",
        path: "/cinemol"
    },
    {
        name: "Biosynfoni",
        icon: <FiHexagon />,
        description: "",
        path: "/biosynfoni"
    }
];

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
                <Route path="/changelog" element={<Changelog />} />
                <Route path="/retromol" element={<RetroMol />} />
                <Route path="/cinemol" element={<CineMol />} />
                <Route path="/biosynfoni" element={<Biosynfoni />} />
                <Route path="/*" element={<div>404</div>} />
            </Routes>
        </div>
    );
};

function App () {

    // Display toast message on first load.
    useEffect (() => {
        toast.info("Select an app tile to get started!");
    }, []);

    // Ping API to see if it is available.
    useEffect(() => {
        fetch("/api/ping_server")
            .then(response => response.json())
            .then(data => {
                if (data.status === "success") {
                    toast.success("API is available!");
                } else {
                    console.log(data);
                    toast.error("API is unavailable!");
                };
            })
            .catch(error => {
                console.log(error);
                toast.error("API is unavailable!");
            });
    }, []);

    return (
        <div className="app">
            <BrowserRouter>
                <Header navbarLinks={[
                    { name: "Home", path: "/" },
                    { name: "About", path: "/about" },
                    { name: "Change log", path: "/changelog" }
                ]}/>
                <AppRoutes />
                <Footer link={"https://github.com/moltools/MolTools"} />  
                <Toast />
            </BrowserRouter>    
        </div>
    );
};

// Render the app.
const root = ReactDOM.createRoot(document.getElementById("root"));

root.render(<App />);