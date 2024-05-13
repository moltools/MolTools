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
    const [selectedMatchType, setSelectedMatchType] = useState("pairwise"); // pairwise or query
    const [selectedPairwiseAlgorithm, setSelectedPairwiseAlgorithm] = useState("global"); // local or global
    const [gapPenalty, setGapPenalty] = useState(2);
    const [endGapPenalty, setEndGapPenalty] = useState(1);
    const [hasNoLeadingModules, setHasNoLeadingModules] = useState(false);
    const [hasNoTrailingModules, setHasNoTrailingModules] = useState(false);
    const [allBioactivityLabels, setAllBioactivityLabels] = useState([]);
    const [matchAgainstMolecules, setMatchAgainstMolecules] = useState(true);
    const [matchAgainstProtoClusters, setMatchAgainstProtoClusters] = useState(false);
    const [selectedBioactivityLabels, setSelectedBioactivityLabels] = useState([]);
    const [minMatchLength, setMinMatchLength] = useState(1);
    const [maxMatchLength, setMaxMatchLength] = useState(100);
    const [numMatchesToReturn, setNumMatchesToReturn] = useState(50);

    const matchDatabase = async () => {
        if (!matchAgainstMolecules && !matchAgainstProtoClusters) {
            toast.warn("No query target selected!");
            return;
        };

        setIsLoading(true);
        setMatches([]);

        const data = {
            "matchItems": queryItems,
            "selectedMatchType": selectedMatchType,
            "selectedPairwiseAlgorithm": selectedPairwiseAlgorithm,
            "gapPenalty": gapPenalty,
            "endGapPenalty": endGapPenalty,
            "hasNoLeadingModules": hasNoLeadingModules,
            "hasNoTrailingModules": hasNoTrailingModules,
            "selectedBioactivityLabels": selectedBioactivityLabels,
            "matchAgainstMolecules": matchAgainstMolecules,
            "matchAgainstProtoClusters": matchAgainstProtoClusters,
            "minMatchLength": minMatchLength, 
            "maxMatchLength": maxMatchLength,
            "numMatchesToReturn": numMatchesToReturn
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
            <div 
                className="panel-heading"
                style={{ borderRadius: "0px" }}
            >
                <div className="title is-6">
                    Options
                </div>
            </div>
            <div className="panel-block">
                <div className="column is-full">
                    <div className="field has-addons">
                        <div className="control">
                            <label className="radio">
                                <input
                                    type="radio"
                                    name="matchType"
                                    value="pairwise"
                                    checked={selectedMatchType === "pairwise"}
                                    onChange={() => setSelectedMatchType("pairwise")}
                                />  
                                <span style={{ marginLeft: "10px" }}>
                                    Match (pairwise alignment)
                                </span>
                            </label>
                        </div>
                    </div>
                    <div className="field has-addons">
                        <div className="control">
                        <label className="radio">
                                <input
                                    type="radio"
                                    name="matchType"
                                    value="query"
                                    checked={selectedMatchType === "query"}
                                    onChange={() => setSelectedMatchType("query")}
                                />  
                                <span style={{ marginLeft: "10px" }}>
                                    Query
                                </span>
                            </label>
                        </div>
                    </div>
                </div>
            </div>
            {selectedMatchType === "pairwise" ? (
                <div className="panel-block">
                    <div className="column is-full">
                        <div className="field has-addons">
                            <div className="control">
                                <label className="radio">
                                    <input
                                        type="radio"
                                        name="selectedPairwiseAlgorithm"
                                        value="local"
                                        checked={selectedPairwiseAlgorithm === "local"}
                                        onChange={() => setSelectedPairwiseAlgorithm("local")}
                                        disabled={true}
                                    />  
                                    <span style={{ marginLeft: "10px" }}>
                                        Local alignment strategy
                                    </span>
                                </label>
                            </div>
                        </div>
                        <div className="field has-addons">
                            <div className="control">
                            <label className="radio">
                                    <input
                                        type="radio"
                                        name="selectedPairwiseAlgorithm"
                                        value="global"
                                        checked={selectedPairwiseAlgorithm === "global"}
                                        onChange={() => setSelectedPairwiseAlgorithm("global")}
                                    />  
                                    <span style={{ marginLeft: "10px" }}>
                                        Global alignment strategy
                                    </span>
                                </label>
                            </div>
                        </div>
                        <div style={{ paddingTop: "10px" }}>
                            <span>
                                Gap penalty:
                            </span>
                            <div className="field has-addons">
                                <div className="control" style={{ maxWidth: "250px", paddingTop: "10px" }}>
                                    <NumberPicker
                                        value={gapPenalty}
                                        defaultValue={gapPenalty}
                                        min={-10}
                                        max={10}
                                        step={1}  
                                        onChange={(value) => setGapPenalty(value)}
                                    />
                                </div>
                            </div>
                        </div>
                        <div style={{ paddingTop: "10px" }}>
                            <span>
                                End gap penalty:
                            </span>
                            <div className="field has-addons">
                                <div className="control" style={{ maxWidth: "250px", paddingTop: "10px" }}>
                                    <NumberPicker
                                        value={endGapPenalty}
                                        defaultValue={endGapPenalty}
                                        min={-10}
                                        max={10}
                                        step={1}  
                                        onChange={(value) => setEndGapPenalty(value)}
                                    />
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            ) : null}
            {selectedMatchType === "query" ? (
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
                                Query has no leading modules
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
                                Query has no trailing modules
                                </label>
                            </div>
                        </div>
                    </div>
                </div>
            ) : null}
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
                                maxWidth: "250px" 
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
                            Include parsed molecules in search space
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
                            Include parsed proto-clusters in search space
                            </label>
                        </div>
                    </div>
                    <div style={{ paddingTop: "10px" }}>
                        <span>
                            Min match length:
                        </span>
                        <div className="field has-addons">
                            <div className="control" style={{ maxWidth: "250px", paddingTop: "10px" }}>
                                <NumberPicker
                                    value={minMatchLength}
                                    defaultValue={minMatchLength}
                                    min={1}
                                    max={100}
                                    step={1}  
                                    onChange={(value) => setMinMatchLength(value)}
                                />
                            </div>
                        </div>
                    </div>
                    <div style={{ paddingTop: "10px" }}>
                        <span>
                            Max match length:
                        </span>
                        <div className="field has-addons">
                            <div className="control" style={{ maxWidth: "250px", paddingTop: "10px" }}>
                                <NumberPicker
                                    value={maxMatchLength}
                                    defaultValue={maxMatchLength}
                                    min={1}
                                    max={100}
                                    step={1}  
                                    onChange={(value) => setMaxMatchLength(value)}
                                />
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            <div className="panel-block">
                <div className="column is-full">
                    <div className="field is-grouped is-grouped-left">
                        <button 
                            className="button is-link is-light" 
                            onClick={() => {
                                if (queryItems.length > 0) {
                                    matchDatabase();
                                } else {
                                    toast.warn("Query is empty!");
                                };
                            }}
                        >
                            Submit
                        </button>
                        {/* <InfoPopup infoText={"1. Match query cannot be ambiguous.\n2. Match query will not show up in results."} /> */}
                    </div>
                </div>
            </div>
        </div>
    );
};

export default MatchSubmission;