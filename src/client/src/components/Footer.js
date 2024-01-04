import React from "react";
import { FiGithub } from "react-icons/fi";
import { Link } from "react-router-dom";

const Footer = ({link}) => {
    // Link opens in new tab.
    return (
        <div className="footer">
            <div>
                <Link to={link} target="_blank" rel="noopener noreferrer">
                    <FiGithub size={20} />
                </Link>
            </div>
        </div>
    );
};

export default Footer;