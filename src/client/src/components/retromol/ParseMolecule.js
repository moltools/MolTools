import React, { useState } from "react";
import { toast } from "react-toastify";
import { v4 as uuidv4 } from "uuid";

import AlignmentTable from "./AlignmentTable";
import LoadingOverlay from "./LoadingOverlay";
import MatchSubmission from "./MatchSubmission";
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
    const [selectedBioactivityLabels, setSelectedBioactivityLabels] = useState([]);
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

        const data = { 
            "smiles": smiles,
            "selectedBioactivityLabels": selectedBioactivityLabels,
            "matchAgainstMolecules": matchAgainstMolecules,
            "matchAgainstProtoClusters": matchAgainstProtoClusters
        };

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
        };

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
                        submissionElement={
                            <MatchSubmission 
                                queryItems={selectedResult} 
                                findMatches={findMatches} 
                                matchAgainstMolecules={matchAgainstMolecules} 
                                setMatchAgainstMolecules={setMatchAgainstMolecules} 
                                matchAgainstProtoClusters={matchAgainstProtoClusters} 
                                setMatchAgainstProtoClusters={setMatchAgainstProtoClusters}
                                selectedBioactivityLabels={selectedBioactivityLabels}
                                setSelectedBioactivityLabels={setSelectedBioactivityLabels}
                            />
                        }
                    />}
            </div>
        </div>
    );
};

export default ParseMolecule;