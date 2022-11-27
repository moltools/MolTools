import React from "react";
import { Link } from "react-router-dom";

function Navbar() {
    return (
        <div className="navbar">
            <Link className="navbar-item" to="/">Home</Link>
        </div>
    );
};

export default Navbar;
