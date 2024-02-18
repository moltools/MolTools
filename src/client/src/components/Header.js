import React, { useState, useEffect } from "react";

const Header = ({ widgetRoutes }) => {
    const [apiAvailable, setApiAvailable] = useState(false);

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

    useEffect(() => {
        checkApi();
        }
    );

    return (        
        <nav 
            class="navbar has-navbar-fixed-top"
            role="navigation" 
            aria-label="main-navigation"
            style={{top: 0, width: "100%", zIndex: 100, boxShadow: "0 0 10px rgba(0, 0, 0, 0.5)"}}
        >
            <div class="navbar-brand">
                <a class="navbar-item" href="/" style={{cursor: "pointer"}}>
                    <h1 style={{color: "black",fontWeight: "bold"}}>MolTools</h1>
                </a>
                <a 
                    role="button" class="navbar-burger" aria-label="menu" aria-expanded="tr" data-target="navbarBasicExample"
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
            <div id="navbarBasicExample" class="navbar-menu">
                <div class="navbar-start">
                    <a class="navbar-item" href="/">Home</a>
                    <a class="navbar-item" href="/about">About</a>
                    <div class="navbar-item has-dropdown is-hoverable">
                        <a class="navbar-link">More</a>
                        <div class="navbar-dropdown">
                            {widgetRoutes.map((widget, index) => (
                                <a class="navbar-item" href={widget.path}>{widget.name}</a>
                            ))}
                            <hr class="navbar-divider" />
                            <a class="navbar-item" href="https://github.com/moltools/MolTools/issues" target="_blank" rel="noopener noreferrer">Report an issue</a>
                        </div>
                    </div>
                </div>
                <div class="navbar-end" style={{margin: "1rem"}}>
                    <div 
                        class={`tag ${apiAvailable ? "is-success" : "is-danger"}`} 
                        style={{transition: "background-color 0.5s ease"}}
                    >
                        Server is {apiAvailable ? "online" : "offline"}
                    </div>  
                </div>
            </div>
        </nav>
    );
};

export default Header;