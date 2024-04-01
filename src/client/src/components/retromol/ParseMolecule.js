import React, { useState } from "react";
import { toast } from "react-toastify";
import { v4 as uuidv4 } from "uuid";

import AlignmentTable from "./AlignmentTable";
import LoadingOverlay from "./LoadingOverlay";
import Modal from "./Modal";
import QueryBuilder from "./QueryBuilder";
import ResultSelector from "./ResultSelector";
import SmilesInput from "./SmilesInput";

const ParseMolecule = () => {
    // Page state.
    const [isLoading, setIsLoading] = useState(false);
    const [modalActive, setModalActive] = useState(false);

    // Query and result state.
    const [smiles, setSmiles] = useState("");
    const [matchAgainstMolecules, setMatchAgainstMolecules] = useState(true);
    const [matchAgainstProtoClusters, setMatchAgainstProtoClusters] = useState(false);
    const [results, setResults] = useState([]);
    const [selectedResult, setSelectedResult] = useState([]);
    const [matches, setMatches] = useState([]);

    // When setResults is called, give every item in every sequence a unique id.
    const tagResults = (results) => {
        return results.map((seq) => {
            return seq.map((item) => {
                return { ...item, id: uuidv4() };
            });
        });
    };

    // Toggle the results modal.
    const toggleModal = () => {
        setModalActive(!modalActive);
    };

    // Clear the input fields.
    const handleRefresh = () => {
        setSmiles("");
        setMatchAgainstMolecules(true);
        setMatchAgainstProtoClusters(false);
        setResults([]);
        setSelectedResult([]);
        setMatches([]);
    };

    // Parse the input SMILES string.
    const parseMolecule = async () => {
        setIsLoading(true);
        setMatches([]);
        setResults([]);
        setSelectedResult([]);

        const data = { "smiles": smiles };

        try {
            const response = await fetch("/api/parse_retromol", {
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
                console.log(json.payload);
                if (json.payload.sequences) setResults(tagResults(json.payload.sequences));
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

    // Query the database.
    const findMatches = async () => {
        if (!matchAgainstMolecules && !matchAgainstProtoClusters) {
            toast.warn("No query target selected!");
            return;
        };

        setIsLoading(true);
        setMatches([]);

        const data = {
            "matchItems": selectedResult,
            "matchAgainstMolecules": matchAgainstMolecules,
            "matchAgainstProtoClusters": matchAgainstProtoClusters
        };

        try {
            const response = await fetch("/api/find_matches", {
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
                console.log(json.payload)
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

    // Submission element for matching against the database.
    const matchSubmission = (queryItems) => {
        return (
            <div>
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

    return (
        <div>
            {isLoading && <LoadingOverlay />}
            <div className="column is-full">
                <Modal 
                    closeModal={toggleModal} 
                    modalState={modalActive} 
                    title="Matching results"
                >
                    <div>
                        {matches.length === 0 ? (
                            <p>Nothing to see here.</p>
                        ) : (
                            <AlignmentTable data={matches} />
                        )}
                    </div>
                </Modal>
                <SmilesInput 
                    smiles={smiles} 
                    setSmiles={setSmiles} 
                    parseMolecule={parseMolecule} 
                    handleRefresh={handleRefresh}
                />
                <ResultSelector
                    results={results}
                    selectedResult={selectedResult}
                    setSelectedResult={setSelectedResult}
                />
                {selectedResult.length > 0 && 
                    <QueryBuilder 
                        queryItems={selectedResult} 
                        setQueryItems={setSelectedResult} 
                        submissionElement={matchSubmission(selectedResult)}
                    />}
            </div>
        </div>
    );
};

export default ParseMolecule;