import React from "react";
import { Link } from "react-router-dom";

const CustomLink = ({ to, display_text }) => (
    <Link style={{ textDecoration: "none" }} to={to}>
        <strong>{display_text}</strong>
    </Link>
);

const link_to_antismash = () => {
    return (
        <CustomLink 
            to="https://antismash.secondarymetabolites.org/#!/start"
            display_text="antiSMASH"
        />
    );
};

const link_to_github = () => {
    return (
        <CustomLink 
            to="https://github.com/moltools/RetroMol.GUI"
            display_text="GitHub repository"
        />
    );
};

const About = () => {
    return (
        <div>
            <h1>About</h1>
            <p>
                RetroMol is a Python library and command line tool for mining 
                biosynthetic motifs from natural product structures.
            </p>
            <p>
                This application was created for users to parse individual 
                natural product structures, and to browse a library of already 
                parsed structures through an easy-to-use interface.
            </p>
            <p>
                Additionally, users can upload results from 
                {" "}{link_to_antismash()} in order to find associated 
                compounds.
            </p>
            <p>
                For more information and for reporting issues, please visit the 
                {" "}{link_to_github()}.
            </p>
        </div>
    );
};

export default About;