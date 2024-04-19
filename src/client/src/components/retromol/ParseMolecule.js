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
    const [results, setResults] = useState([]);
    const [selectedResult, setSelectedResult] = useState([]);
    const [matches, setMatches] = useState([]);

    // When setResults is called, give every item in every sequence a unique id.
    const tagResults = (results) => {
        return results.map((seq) => {
            return seq.map((item) => {
                return [{ ...item, id: uuidv4() }];
            });
        });
    };

    // Toggle the results modal.
    const toggleModal = () => {
        setModalActive(!modalActive);
    };

    // Clear the input fields.
    const clearInputResults = () => {
        setSmiles("");
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
            const response = await fetch("/api/parse_smiles", {
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
                    handleRefresh={clearInputResults}
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
                                setMatches={setMatches}
                                setIsLoading={setIsLoading}
                                toggleModal={toggleModal}
                            />
                        }
                    />}
            </div>
        </div>
    );
};

export default ParseMolecule;