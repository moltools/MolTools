import React from "react";
import { useLocation } from "react-router-dom";
import Navbar from "./Navbar";

/**
 * Header component that displays a header with the current page's display name and a navigation bar.
 * @param {Object} props - The props for the Header component.
 * @param {Object} props.routeDisplayNames - A mapping of route paths to display names.
 * @param {Array} props.navbarLinks - An array of navigation links.
 * @returns {JSX.Element} - The rendered Header component.
 */
const Header = (props) => {
    // Destructure the props.
    const { routeDisplayNames, navbarLinks } = props;

    // Get the current location using the useLocation hook from react-router-dom.
    const location = useLocation();

    // Get the current display name based on the current route's path using the mapping.
    // If the path is not found in the mapping, use an empty string as the default.
    const currentDisplayName = routeDisplayNames[location.pathname] || "";

    return (
        <header className="header">
            <div>
                {/* Display the current display name */}
                <h2>{currentDisplayName}</h2>
            </div>
            <div className="header-navbar">
                {/* Render the Navbar component with the provided links */}
                <Navbar links={navbarLinks} />
            </div>
        </header>
    );
};

// Export the Header component as the default export.
export default Header;