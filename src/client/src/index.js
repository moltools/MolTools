import React, { useState, useEffect } from "react";
import ReactDOM from "react-dom/client";
import { BrowserRouter, Routes, Route, useLocation } from "react-router-dom";
import { FiHexagon } from "react-icons/fi";
import "./style/main.css";
import Header from "./components/Header";
import Toast from "./components/Toast";
import Home from "./pages";
import About from "./pages/about";
import RetroMol from "./pages/retromol";
import CineMol from "./pages/cinemol";
import Biosynfoni from "./pages/biosynfoni";
import PIKAChU from "./pages/pikachu";

const widgets = [
    { 
        name: "RetroMol", 
        description: "Discover biosynthetically similar molecules", 
        path: "/retromol", 
        component: <RetroMol />
    },
    // { 
    //     name: "CineMol", 
    //     description: "Draw 3D molecule model", 
    //     path: "/cinemol",
    //     component: <CineMol />
    // },
    // {
    //     name: "Biosynfoni", 
    //     description: "Predict biosynthetic class", 
    //     path: "/biosynfoni",
    //     component: <Biosynfoni />
    // },
    // { 
    //     name: "PIKAChU", 
    //     description: "Draw 2D structural formula", 
    //     path: "/pikachu",
    //     component: <PIKAChU /> 
    // },
];

function AppRoutes () {
    const location = useLocation();
    const [displayLocation, setDisplayLocation] = useState(location);
    const [transitionStage, setTransitionStage] = useState("fade-in");

    useEffect(() => {
        if (location !== displayLocation) setTransitionStage("fade-out");
    }, [location, displayLocation]);

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
                <Route path="/" element={<Home widgets={widgets}/>} />
                <Route path="/about" element={<About />} />
                <Route path="/*" element={<div style={{padding: "2rem"}}><h1>404 Not Found</h1></div>} />
                {widgets.map((widget, index) => (
                    <Route key={index} path={widget.path} element={widget.component} />
                ))}
            </Routes>
        </div>
    );
};

function App () {
    return (
        <div>
            <BrowserRouter>
                <Header widgetRoutes={widgets} />
                <AppRoutes />
                <Toast />
            </BrowserRouter>    
        </div>
    );
};

const root = ReactDOM.createRoot(document.getElementById("root"));

root.render(<App />);