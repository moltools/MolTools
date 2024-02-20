import React from "react";
import { Link } from "react-router-dom";

const Home = ( { widgets }) => {
    return (
        <div style={{padding: "20px"}}>
            {widgets.map((widget, index) => (
                <div 
                    key={index} 
                    class="card"
                    style={{marginBottom: "20px", boxShadow: "0 0 10px rgba(0, 0, 0, 0.25)"}}
                >
                    <div class="card-content">
                        <div class="media">
                            <div class="media-left">
                                <figure class="image is-64x64">
                                    <img src={widget.logo} alt={widget.name} />
                                </figure>
                            </div>
                            <div class="media-content">
                                <p class="title is-4">{widget.name}</p>
                                <p class="subtitle is-6">{widget.description}</p>
                            </div>
                        </div>
                        <div class="content">
                            <div class="buttons" style={{marginBottom: 0}}>
                                {widget.links.map((link, index) => (
                                    <p class="control" key={index} style={{marginTop: 5, marginBottom: 5, marginRight: 5}}>
                                        <a 
                                            href={link.url} 
                                            target="_blank"
                                            rel="noopener noreferrer"
                                            class="button is-link is-light"
                                        >
                                            <span class="icon">{link.icon}</span>
                                            <span>{link.name}</span>
                                        </a>
                                    </p>
                                ))}
                            </div>
                            <p>
                                <Link to={widget.path} class="button is-link is-light">
                                    Open {widget.name} widget
                                </Link>
                            </p>
                        </div>
                    </div>
                </div>
            ))}
        </div>
    );
};

export default Home;