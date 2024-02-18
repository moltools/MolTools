import React from "react";
import { Link } from "react-router-dom";
import { FiImage } from "react-icons/fi";

const Home = ( { widgets }) => {

    return (
        <div style={{padding: "20px"}}>
            {widgets.map((widget, index) => (
                <Link to={widget.path}>
                    <div 
                        key={index} 
                        class="card"
                        style={{marginBottom: "20px", boxShadow: "0 0 10px rgba(0, 0, 0, 0.25)"}}
                    >
                        <div class="card-content">
                            <div class="media">
                                <div class="media-left">
                                    <FiImage size={48} />
                                </div>
                                <div class="media-content">
                                    <p class="title is-4">{widget.name}</p>
                                    <p class="subtitle is-6">{widget.description}</p>
                                </div>
                            </div>
                            <div class="content">
                                No content to display
                            </div>
                        </div>
                    </div>
                </Link>
            ))}
        </div>
    );
};

export default Home;