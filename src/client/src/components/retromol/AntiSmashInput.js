import React from "react";
import { toast } from "react-toastify";

const JobIdInput = ({ jobId, setJobId }) => {
    return (
        <div 
            className="panel-block" 
            style={{ width: "100%" }}
        >
            <div className="column is-full">
                <div className="field is-grouped is-grouped-left">
                    <div className="control" style={{ width: "100%" }}>
                        <div className="field">
                            <div className="control">
                               <input
                                   className="input"
                                   value={jobId}
                                   type="text"
                                   placeholder="Enter job ID"
                                   onChange={(e) => setJobId(e.target.value)}
                               />
                           </div>
                       </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

const JsonInput = ({ jsonSrc, setJsonSrc }) => {

    const handleFileChange = (event) => {
        // const file = event.target.files[0];
        // const reader = new FileReader();

        // reader.onload = (event) => {
        //     setJsonSrc(event.target.result);
        // };

        // reader.readAsText(file);

        toast.error("JSON input is not supported yet!");
    };

    return (
        <div 
            className="panel-block"
            style={{ width: "100%" }}
        >
            <div className="column is-full">
                <div className="field is-grouped is-grouped-left">
                    <div 
                        className="control" 
                        style={{ width: "100%" }}
                    >
                        <input
                            type="file"
                            className="button"
                            style={{width: "100%"}}
                            onChange={handleFileChange}
                        />
                    </div>
                </div>
            </div>
        </div>
    );
};

const AntiSmashInput = ({ 
    selectedInputType, 
    setSelectedInputType,
    jobId,
    setJobId,
    jsonSrc,
    setJsonSrc,
    parseInput,
    handleRefresh
}) => {
    const handleInputTypeChange = (e) => {
        setJobId("");
        setJsonSrc("");
        setSelectedInputType(e.target.value);
    };

    return (
        <div
            className="control"
            style={{
                border: "1px solid #dbdbdb",
                borderRadius: "5px",
                marginBottom: "10px"
            
            }}
        >   
            <div className="panel">
                <div className="panel-heading">
                    <div className="title is-5">
                        Input
                    </div>
                </div>
                <div className="panel-block">
                    <div
                        className="field"
                        style={{
                            width: "100%",
                            margin: "10px"
                        }}
                    >   
                        <div className="column is-full">
                            <div className="field has-addons">
                                <div className="control">
                                    <label className="radio">
                                        <input
                                            type="radio"
                                            name="inputType"
                                            value="jobId"
                                            checked={selectedInputType === "jobId"}
                                            onChange={handleInputTypeChange}
                                        />
                                        <span style={{ marginLeft: "5px" }}>
                                            Job ID
                                        </span>
                                    </label>
                                </div>
                            </div>
                            <div className="field has-addons">
                                <div className="control">
                                    <label className="radio">
                                        <input
                                            type="radio"
                                            name="inputType"
                                            value="json"
                                            checked={selectedInputType === "json"}
                                            onChange={handleInputTypeChange}
                                        />
                                        <span style={{ marginLeft: "5px" }}>
                                            JSON
                                        </span>
                                    </label>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>  
                <div className="panel-block">
                    {selectedInputType === "jobId" ? (
                        <JobIdInput 
                            jobId={jobId} 
                            setJobId={setJobId} 
                        />
                    ) : (
                        <JsonInput 
                            jsonSrc={jsonSrc} 
                            setJsonSrc={setJsonSrc} 
                        />
                    )}
                </div>
                <div className="panel-block">
                    <div
                        className="field has-addons"
                        style={{ margin: "10px" }}
                    >
                        <div className="control">
                            <button
                                className="button is-link is-light"
                                style={{
                                    marginRight: "5px",
                                    marginBottom: "5px"
                                }}
                                onClick={parseInput}
                            >
                                Parse input
                            </button>
                            <button
                                className="button is-link is-light"
                                style={{
                                    marginRight: "5px",
                                    marginBottom: "5px"
                                }}
                                onClick={handleRefresh}
                            >
                                Refresh
                            </button>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default AntiSmashInput;