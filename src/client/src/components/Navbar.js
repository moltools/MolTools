import React from "react";
import { Link } from "react-router-dom";

/**
 * Navbar component that renders a navigation bar with links.
 * @param {Object} props - The props for the Navbar component.
 * @param {Array} props.links - An array of navigation links.
 * @returns {JSX.Element} - The rendered Navbar component.
 */
const Navbar = (props) => {
    // Destructure the props.
    const { links } = props;

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

// Export the Navbar component as the default export.
export default Navbar;