import React from "react";
import { FiGithub } from "react-icons/fi";
import { Link } from "react-router-dom";

const Footer = ({link}) => {
    return (
        <div className="footer">
            <Link to={link}>
                <FiGithub size={20}/>
            </Link>
        </div>
    );
};

export default Footer;