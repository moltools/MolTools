import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";
import MultiSelect from "react-widgets/Multiselect";
import NumberPicker from "react-widgets/NumberPicker";
import InfoPopup from "./InfoPopUp";

const MatchSubmission = ({ 
    queryItems, 
    setMatches,
    setIsLoading,
    toggleModal
}) => {
    const [hasNoLeadingModules, setHasNoLeadingModules] = useState(false);
    const [hasNoTrailingModules, setHasNoTrailingModules] = useState(false);
    const [allBioactivityLabels, setAllBioactivityLabels] = useState([]);
    const [matchAgainstMolecules, setMatchAgainstMolecules] = useState(true);
    const [matchAgainstProtoClusters, setMatchAgainstProtoClusters] = useState(false);
    const [selectedBioactivityLabels, setSelectedBioactivityLabels] = useState([]);
    const [numMatchesToReturn, setNumMatchesToReturn] = useState(50);

    const matchDatabase = async ( ambiguousNotAllowed ) => {
        if (!matchAgainstMolecules && !matchAgainstProtoClusters) {
            toast.warn("No query target selected!");
            return;
        };

        setIsLoading(true);
        setMatches([]);

        const data = {
            "matchItems": queryItems,
            "ambiguousNotAllowed": ambiguousNotAllowed,
            "selectedBioactivityLabels": selectedBioactivityLabels,
            "matchAgainstMolecules": matchAgainstMolecules,
            "matchAgainstProtoClusters": matchAgainstProtoClusters,
            "numMatchesToReturn": numMatchesToReturn,
            "hasNoLeadingModules": hasNoLeadingModules,
            "hasNoTrailingModules": hasNoTrailingModules
        };

        try {
            const response = await fetch("/api/match_database", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ data })
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            if (json.status === "success") {
                toast.success(json.message);
                setMatches(json.payload.matches);
                toggleModal();
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            toast.error(error.message);
        }
        setIsLoading(false);
    };

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
                <div className="column is-full">
                    <div className="field has-addons">
                        <div className="control">
                            <label className="checkbox">
                            <input 
                                type="checkbox" 
                                defaultChecked={hasNoLeadingModules}
                                value={hasNoLeadingModules}
                                onChange={() => setHasNoLeadingModules(!hasNoLeadingModules)}
                            />
                            Has no leading modules (for querying only)
                            </label>
                        </div>
                    </div>
                    <div className="field has-addons">
                        <div className="control">
                            <label className="checkbox">
                            <input 
                                type="checkbox" 
                                defaultChecked={hasNoTrailingModules}
                                value={hasNoTrailingModules}
                                onChange={() => setHasNoTrailingModules(!hasNoTrailingModules)}
                            />
                            Has no trailing modules (for querying only)
                            </label>
                        </div>
                    </div>
                </div>
            </div>
            <div className="panel-block">
                <div 
                    style={{ 
                        width: "100%", 
                        marginTop: "10px", 
                        marginBottom: "10px" 
                    }}
                >
                    <span style={{ margin: "10px" }}>
                        Filter on known bioactivity labels:
                    </span>
                    <div 
                        className="field has-addons" 
                        style={{ width: "100%" }}
                    >
                        <div 
                            className="control" 
                            style={{ 
                                margin: "10px", 
                                width: "100%" 
                            }}
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
                <div 
                    style={{ 
                        marginTop: "10px", 
                        marginBottom: "10px" 
                    }}
                >
                    <span style={{ margin: "10px" }}>
                        Maximum number of matches to report:
                    </span>
                    <div className="field has-addons">
                        <div 
                            className="control" 
                            style={{ 
                                margin: "10px", 
                                width: "100%" 
                            }}
                        >
                            <NumberPicker
                                value={numMatchesToReturn}
                                defaultValue={50}
                                min={1}
                                max={100}
                                step={1}  
                                onChange={(value) => setNumMatchesToReturn(value)}
                            />
                        </div>
                    </div>
                </div>
            </div>
            <div className="panel-block">
                <div className="column is-full">
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
            <div className="panel-block">
                <div className="column is-full">
                    <div className="field is-grouped is-grouped-left">
                        <div className="buttons">
                            <button 
                                className="button is-link is-light" 
                                onClick={() => {
                                    if (queryItems.length > 0) {
                                        matchDatabase(true);
                                    } else {
                                        toast.warn("Query is empty!");
                                    }
                                }}
                            >
                                Match against database
                            </button>
                            <button 
                                className="button is-link is-light" 
                                onClick={() => {
                                    if (queryItems.length > 0) {
                                        matchDatabase(false);
                                    } else {
                                        toast.warn("Query is empty!");
                                    }
                                }}
                            >
                                Query against database
                            </button>
                            <InfoPopup infoText={"1. Match query cannot be ambiguous.\n2. Match query will not show up in results."} />
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default MatchSubmission;