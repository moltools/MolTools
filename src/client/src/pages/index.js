import React from "react";
import { Link } from "react-router-dom";

const Home = ( { widgets }) => {
    return (
        <div style={{padding: "20px"}}>
            {widgets.map((widget, index) => (
                <div 
                    key={index} 
                    className="card"
                    style={{marginBottom: "20px", boxShadow: "0 0 10px rgba(0, 0, 0, 0.25)"}}
                >
                    <div className="card-content">
                        <div className="media">
                            <div className="media-left">
                                <figure className="image is-64x64">
                                    <img src={widget.logo} alt={widget.name} />
                                </figure>
                            </div>
                            <div className="media-content">
                                <p className="title is-4">{widget.name}</p>
                                <p className="subtitle is-6">{widget.description}</p>
                            </div>
                        </div>
                        <div className="content">
                            <div className="buttons" style={{marginBottom: 0}}>
                                {widget.links.map((link, index) => (
                                    <p className="control" key={index} style={{marginTop: 5, marginBottom: 5, marginRight: 5}}>
                                        <a 
                                            href={link.url} 
                                            target="_blank"
                                            rel="noopener noreferrer"
                                            className="button is-link is-light"
                                        >
                                            <span className="icon">{link.icon}</span>
                                            <span>{link.name}</span>
                                        </a>
                                    </p>
                                ))}
                            </div>
                            <p>
                                <Link to={widget.path} className="button is-link is-light">
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