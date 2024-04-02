import React, { useState } from "react";
import { toast } from "react-toastify";

import AlignmentTable from "./AlignmentTable";
import LoadingOverlay from "./LoadingOverlay";
import MatchSubmission from "./MatchSubmission";
import Modal from "./Modal";
import QueryBuilder from "./QueryBuilder";

const QueryDatabase = () => {
    // Page state.
    const [isLoading, setIsLoading] = useState(false);
    const [modalActive, setModalActive] = useState(false);

    // Query and result state.
    const [selectedBioactivityLabels, setSelectedBioactivityLabels] = useState([]);
    const [queryAgainstMolecules, setQueryAgainstMolecules] = useState(true);
    const [queryAgainstProtoClusters, setQueryAgainstProtoClusters] = useState(false);
    const [queryItems, setQueryItems] = useState([]);
    const [matches, setMatches] = useState([]);

    // Toggle the results modal.
    const toggleModal = () => {
        setModalActive(!modalActive);
    };

    // Query the database.
    const queryDatabase = async () => {
        if (!queryAgainstMolecules && !queryAgainstProtoClusters) {
            toast.warn("No query target selected!");
            return;
        };

        setIsLoading(true);
        setMatches([]);

        const data = {
            "queryItems": queryItems,
            "selectedBioactivityLabels": selectedBioactivityLabels,
            "queryAgainstMolecules": queryAgainstMolecules,
            "queryAgainstProtoClusters": queryAgainstProtoClusters
        };

        try {
            const response = await fetch("/api/query_database", {
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
                <QueryBuilder 
                    queryItems={queryItems} 
                    setQueryItems={setQueryItems} 
                    submissionElement={
                        <MatchSubmission 
                            queryItems={queryItems} 
                            findMatches={queryDatabase} 
                            matchAgainstMolecules={queryAgainstMolecules} 
                            setMatchAgainstMolecules={setQueryAgainstMolecules} 
                            matchAgainstProtoClusters={queryAgainstProtoClusters} 
                            setMatchAgainstProtoClusters={setQueryAgainstProtoClusters} 
                            selectedBioactivityLabels={selectedBioactivityLabels}
                            setSelectedBioactivityLabels={setSelectedBioactivityLabels}
                        />}
                />
            </div>
        </div>
    );
};

export default QueryDatabase;