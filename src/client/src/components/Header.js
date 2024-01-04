import React from "react";
import { BrowserRouter, Routes, Route, useLocation } from "react-router-dom";
import Navbar from "./Navbar";

const Header = ({ routeDisplayNames, navbarLinks }) => {

    // Get the display name based on the current route's path using the mapping.
    const location = useLocation();
    const currentDisplayName = routeDisplayNames[location.pathname] || "";

    return (
        <header className="header">
            <div>
                <h2>{currentDisplayName}</h2>
            </div>
            <div className="header-navbar">
                <Navbar links={navbarLinks} />
            </div>
        </header>
    );
};

export default Header;