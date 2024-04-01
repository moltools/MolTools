import React from "react";

const ResultSelector = ({ results, selectedResult, setSelectedResult }) => {
    return (
        <div 
            className="control" 
            style={{
                border: "1px solid #dbdbdb", 
                borderRadius: "5px", 
                marginBottom: "10px"
            }}
        >
            <div className="panel">
                <div className="panel-heading">
                    <div className="title is-5">
                        Results
                    </div>
                </div>
                <div className="panel-block">
                    <div>
                        {results.length === 0 ? (
                            <p style={{ margin: "10px" }}>
                                No results found.
                            </p>
                        ) : (
                            <div 
                                className="field is-grouped is-grouped-multiline" 
                                style={{ 
                                    width: "100%", 
                                    margin: "10px" 
                                }}
                            >
                                {results.map((result, index) => (
                                    <div key={index} className="control">
                                        <div
                                            className={`tags has-addons ${selectedResult === result ? "selected" : ""}`}
                                            style={{ marginRight: "5px", cursor: "pointer" }}
                                            onClick={() => setSelectedResult(result)}
                                        >
                                            <span className={`tag is-link ${selectedResult === result ? "" : "is-light"}`}>
                                                Result {index + 1}
                                            </span>
                                        </div>
                                    </div>
                                ))}
                            </div>
                        )}
                    </div>
                </div>
            </div>
        </div>
    );
};

export default ResultSelector;