import React, { useState } from "react";
import { toast } from "react-toastify";
import { v4 as uuidv4 } from "uuid";

import AlignmentTable from "./AlignmentTable";
import LoadingOverlay from "./LoadingOverlay";
import MatchSubmission from "./MatchSubmission";
import Modal from "./Modal";
import QueryBuilder from "./QueryBuilder";
import ResultSelector from "./ResultSelector";
import AntiSmashInput from "./AntiSmashInput";

// Example job ID with url:
// bacteria-f4aedee3-6f4d-4023-ac53-6cf9263e1625
// https://antismash.secondarymetabolites.org/upload/bacteria-f4aedee3-6f4d-4023-ac53-6cf9263e1625/Y16952.json 

const ParseProtoCluster = () => {
    // Page state.
    const [isLoading, setIsLoading] = useState(false);
    const [modalActive, setModalActive] = useState(false);

    // Input state.
    const [selectedInputType, setSelectedInputType] = useState("jobId"); // Either jobId or json.

    // Query and result state.
    const [jobId, setJobId] = useState("");
    const [ncbiAccession, setNcbiAccession] = useState("");
    const [jsonSrc, setJsonSrc] = useState("");
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
        setJobId("");
        setNcbiAccession("");
        setJsonSrc("");
    };

    // Parse input.
    const parseInput = async () => {
        setIsLoading(true);
        setMatches([]);
        setResults([]);
        setSelectedResult([]);

        if (selectedInputType === "jobId") {
            if (jobId === "") {
                toast.warn("Job ID is empty!");
                setIsLoading(false);
                return;
            };
            if (ncbiAccession === "") {
                toast.warn("NCBI accession is empty!");
                setIsLoading(false);
                return;
            };
        };


        if (selectedInputType === "json") {
            if (jsonSrc === "") {
                toast.warn("JSON source is empty!");
                setIsLoading(false);
                return;
            };
        };

        const data = {
            "selectedInputType": selectedInputType,
            "jobId": jobId,
            "ncbiAccession": ncbiAccession,
            "jsonSrc": jsonSrc
        };

        try {
            const response = await fetch("/api/parse_proto_cluster", {
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
                setResults(tagResults(json.payload.results));
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

    // Parse the input data.
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
                        {matches.length == 0 ? (
                            <p>Nothing to see here.</p>
                        ) : (
                            <AlignmentTable data={matches} />
                        )}
                    </div>
                </Modal>
                <AntiSmashInput 
                    selectedInputType={selectedInputType} 
                    setSelectedInputType={setSelectedInputType}
                    jobId={jobId}
                    setJobId={setJobId}
                    ncbiAccession={ncbiAccession}
                    setNcbiAccession={setNcbiAccession}
                    jsonSrc={jsonSrc}
                    setJsonSrc={setJsonSrc}
                    parseInput={parseInput}
                    handleRefresh={handleRefresh}
                />
                <ResultSelector
                    results={results}
                    setSelectedResult={setMatches}
                    tagResults={tagResults}
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

export default ParseProtoCluster;