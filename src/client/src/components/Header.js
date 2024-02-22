import React, { useState, useEffect } from "react";

const Header = ({ widgetRoutes }) => {
    const [apiAvailable, setApiAvailable] = useState(false);
    const [dbAvailable, setDbAvailable] = useState(false);

    const checkApi = async () => {
        try {
            const response = await fetch("/api/ping_server");
            const data = await response.json();
            if (data.status === "success") {
                setApiAvailable(true);
            }
        } catch (error) {
            console.log(error);
        }
    };

    const checkDb = async () => {
        try {
            const response = await fetch("/api/ping_database");
            const data = await response.json();
            if (data.status === "success") {
                console.log(data.message);
                setDbAvailable(true);
            } else {
                console.log(data.message);
            }
        } catch (error) {
            console.log(error);
        }
    };

    useEffect(() => {
        checkApi();
        checkDb();
    });

    return (        
        <nav 
            className="navbar is-fixed-top"
            role="navigation" 
            aria-label="main-navigation"
            style={{top: 0, width: "100%", zIndex: 1000, boxShadow: "0 0 10px rgba(0, 0, 0, 0.5)"}}
        >
            <div className="navbar-brand">
                <a className="navbar-item" href="/" style={{cursor: "pointer"}}>
                    <h1 style={{color: "black",fontWeight: "bold"}}>MolTools</h1>
                </a>
                <a 
                    role="button" className="navbar-burger" aria-label="menu" aria-expanded="tr" data-target="navbarBasicExample"
                    onClick={() => {
                        const burger = document.querySelector(".navbar-burger");
                        const menu = document.querySelector(".navbar-menu");
                        burger.classList.toggle("is-active");
                        menu.classList.toggle("is-active");
                    }}
                >
                    <span aria-hidden="true"></span>
                    <span aria-hidden="true"></span>
                    <span aria-hidden="true"></span>
                </a>
            </div>
            <div id="navbarBasicExample" className="navbar-menu">
                <div className="navbar-start">
                    <a className="navbar-item" href="/">Home</a>
                    <a className="navbar-item" href="/about">About</a>
                    <div className="navbar-item has-dropdown is-hoverable">
                        <a className="navbar-link">More</a>
                        <div className="navbar-dropdown">
                            {widgetRoutes.map((widget, index) => (
                                <a key={index} className="navbar-item" href={widget.path}>{widget.name}</a>
                            ))}
                            <hr className="navbar-divider" />
                            <a className="navbar-item" href="https://github.com/moltools/MolTools/issues" target="_blank" rel="noopener noreferrer">Report an issue</a>
                        </div>
                    </div>
                </div>
                <div className="navbar-end" style={{margin: "1rem"}}>
                    <div 
                        className={`tag ${apiAvailable ? "is-success" : "is-danger"}`} 
                        style={{transition: "background-color 0.5s ease", marginRight: "5px"}}
                    >
                        Server is {apiAvailable ? "online" : "offline"}
                    </div>  
                    <div
                        className={`tag ${dbAvailable ? "is-success" : "is-danger"}`}
                        style={{transition: "background-color 0.5s ease"}}
                    >
                        Database is {dbAvailable ? "online" : "offline"}
                    </div>
                </div>
            </div>
        </nav>
    );
};

export default Header;