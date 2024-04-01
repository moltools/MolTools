import React from "react";
import ReactDOM from "react-dom/client";
import { BrowserRouter, Routes, Route } from "react-router-dom";
import { FiGithub, FiPaperclip } from "react-icons/fi";
import "./style/main.css";

import Header from "./components/general/Header";
import Toast from "./components/general/Toast";

import Home from "./pages";
import About from "./pages/about";
import CineMol from "./pages/cinemol";
import RetroMol from "./pages/retromol";

const widgets = [
    { 
        name: "RetroMol", 
        description: "Discover biosynthetically similar molecules", 
        logo: "/widgets/retromol.svg",
        links: [
            { icon: <FiGithub />, name: "GitHub", url: "https://github.com/moltools/RetroMol"},
        ],
        path: "/retromol", 
        component: <RetroMol />
    },
    { 
        name: "CineMol", 
        description: "Draw 3D molecule model", 
        logo: "/widgets/cinemol.svg",
        links: [
            { icon: <FiGithub />, name: "GitHub", url: "https://github.com/moltools/CineMol"},
            { icon: <FiPaperclip />, name: "Paper", url: "https://chemrxiv.org/engage/chemrxiv/article-details/65bbb3c966c1381729bd6e27"}
        ],
        path: "/cinemol",
        component: <CineMol />
    }
];

function AppRoutes () {
    return (
        <div>
            <Routes>
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