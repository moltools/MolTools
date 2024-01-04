import React from "react";
import { Link } from "react-router-dom";

const Navbar = ({ links }) => {
    return (
        <div className="navbar">
            {links.map((link, index) => (
                <Link
                    className="link"
                    key={index}
                    to={link.path}
                >
                    {link.name}
                </Link>
            ))}
        </div>
    );
};

export default Navbar;