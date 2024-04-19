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
    const [queryItems, setQueryItems] = useState([]);
    const [matches, setMatches] = useState([]);

    // Toggle the results modal.
    const toggleModal = () => {
        setModalActive(!modalActive);
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
                            setMatches={setMatches}
                            setIsLoading={setIsLoading}
                            toggleModal={toggleModal}
                        />}
                />
            </div>
        </div>
    );
};

export default QueryDatabase;