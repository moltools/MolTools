import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";
import MultiSelect from "react-widgets/Multiselect";

const MatchSubmission = ({ 
    queryItems, 
    findMatches, 
    matchAgainstMolecules, 
    setMatchAgainstMolecules, 
    matchAgainstProtoClusters, 
    setMatchAgainstProtoClusters,
    selectedBioactivityLabels,
    setSelectedBioactivityLabels
}) => {
    const [allBioactivityLabels, setAllBioactivityLabels] = useState([]);

    const fetchBioactivityLabels = async () => {
        try {
            const response = await fetch("/api/bioactivity_labels", {
                method: "GET",
                headers: { "Content-Type": "application/json" }
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            if (json.status === "success") {
                setAllBioactivityLabels(json.payload.bioactivityLabels);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            toast.error(error.message);
        };
    };

    useEffect(() => {
        fetchBioactivityLabels();
    }, []);

    return (
        <div>
            <div className="panel-block">
                <div 
                    style={{ 
                        width: "100%", 
                        marginTop: "10px", 
                        marginBottom: "10px" 
                    }}
                >
                    <span style={{ margin: "10px" }}>
                        Filter results on known bioactivity labels:
                    </span>
                    <div 
                        className="field has-addons" 
                        style={{ width: "100%" }}
                    >
                        <div 
                            className="control" 
                            style={{ margin: "10px", width: "100%" }}
                        >
                            <MultiSelect
                                defaultValue={selectedBioactivityLabels}
                                value={selectedBioactivityLabels}
                                data={allBioactivityLabels}
                                placeholder="Select bioactivity labels (optional)"
                                onChange={(value) => setSelectedBioactivityLabels(value)}
                            />
                        </div>
                    </div>
                </div>
            </div>
            <div className="panel-block">
                <div className="column is-full">
                    <div className="field has-addons">
                        <div className="control">
                            <button 
                                className="button is-link is-light" 
                                onClick={() => {
                                    if (queryItems.length > 0) {
                                        findMatches();
                                    } else {
                                        toast.warn("Query is empty!");
                                    }
                                }}
                            >
                                Match against database
                            </button>
                        </div>
                    </div>
                    <div className="field has-addons">
                        <div className="control">
                            <label className="checkbox">
                            <input 
                                type="checkbox" 
                                defaultChecked={matchAgainstMolecules}
                                value={matchAgainstMolecules}
                                onChange={() => setMatchAgainstMolecules(!matchAgainstMolecules)}
                            />
                            Match against parsed molecules
                            </label>
                        </div>
                    </div>
                    <div className="field has-addons">
                        <div className="control">
                            <label className="checkbox">
                            <input 
                                type="checkbox" 
                                defaultChecked={matchAgainstProtoClusters}
                                value={matchAgainstProtoClusters}
                                onChange={() => setMatchAgainstProtoClusters(!matchAgainstProtoClusters)}
                            />
                            Match against parsed proto-clusters
                            </label>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    )
};

export default MatchSubmission;