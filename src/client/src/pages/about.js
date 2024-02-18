import React from "react";

// =====================================================================================================================
// About component.
// =====================================================================================================================

/**
 * About Component
 *
 * This component represents the "About" section of the MolTools application, providing
 * information about the purpose and functionality of the application.
 *
 * @returns {JSX.Element} The rendered About component containing information.
 */
const About = () => {
    return (
        <div className="widget-content">
            <p>
                MolTools is a collection of tools for the visualization and analysis of natural product compounds.
            </p>
        </div>
    );
};

// Export the About component as the default export.
export default About;