import React from "react";
import Navbar from "./Navbar";

const Header = ({ navbarLinks }) => {
    return (
        <header className="header">
            {/* <div className="header-logo">
                <img src="/logo.svg" alt="logo" />
            </div> */}
            <div className="header-navbar">
                <Navbar links={navbarLinks} />
            </div>
        </header>
    );
};

export default Header;