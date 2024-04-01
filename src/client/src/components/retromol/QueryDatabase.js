import React, { useState } from "react";
import { toast } from "react-toastify";

import AlignmentTable from "./AlignmentTable";
import LoadingOverlay from "./LoadingOverlay";
import Modal from "./Modal";
import QueryBuilder from "./QueryBuilder";

const QueryDatabase = () => {
    // Page state.
    const [isLoading, setIsLoading] = useState(false);
    const [modalActive, setModalActive] = useState(false);

    // Query and result state.
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

    // Query submission element.
    const querySubmission = (queryItems) => {
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
                                            queryDatabase();
                                        } else {
                                            toast.warn("Query is empty!");
                                        }
                                    }}
                                >
                                    Query against database
                                </button>
                            </div>
                        </div>
                        <div className="field has-addons">
                            <div className="control">
                                <label className="checkbox">
                                    <input 
                                        type="checkbox" 
                                        defaultChecked={queryAgainstMolecules}
                                        value={queryAgainstMolecules}
                                        onChange={() => setQueryAgainstMolecules(!queryAgainstMolecules)}
                                    />
                                    Query against parsed molecules
                                </label>
                            </div>
                        </div>
                        <div className="field has-addons">
                            <div className="control">
                                <label className="checkbox">
                                    <input 
                                        type="checkbox" 
                                        defaultChecked={queryAgainstProtoClusters}
                                        value={queryAgainstProtoClusters}
                                        onChange={() => setQueryAgainstProtoClusters(!queryAgainstProtoClusters)}
                                    />
                                    Query against parsed proto-clusters
                                </label>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        );
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
                    submissionElement={querySubmission(queryItems)}
                />
            </div>
        </div>
    );
};

export default QueryDatabase;