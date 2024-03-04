import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";
import { v4 as uuidv4 } from "uuid";
import MatchBuilder from "../components/MatchBuilder";
import QueryBuilder from "../components/QueryBuilder";

const Overview = () => {
    return (
        <div> 
            <h1 className="title">RetroMol</h1>
            <p>
                RetroMol is a tool for the generation of retrosynthetic trees for a given molecule, with the purpose
                of identifying potential precursors for the molecule. The scope of the tool is to provide a user-friendly
                interface for the discovery of retrobiosynthesis of modular natural products. Modular natural products
                currently supported by RetroMol include polyketides, non-ribosomal peptides, and polyketide-nonribosomal
                peptide hybrids.
            </p>
            <br />
            <div>
                You can use RetroMol from two starting inputs:
                <ol style={{paddingLeft: "40px"}}>
                    <li>A molecule represented by a SMILES string</li>
                    <li>An output file from <i>antiSMASH</i> containing mined proto-clusters</li>
                </ol>
            </div>
            <br />
            <p>
                Select the appropriate tab to get started. After parsing the molecule or proto-cluster, you can query
                the database to retrieve similar molecules and their retrosynthetic trees.
            </p>
        </div>
    );
};

const Modal = ({ children, closeModal, modalState, title }) => {
    if(!modalState) {
      return null;
    }
    
    return(
      <div className="modal is-active">
        <div className="modal-background" onClick={closeModal} />
        <div className="modal-card" style={{borderRadius: "10px", width: "90%", height: "90%", paddingTop: "50px"}}>
          <header className="modal-card-head">
            <p className="modal-card-title">{title}</p>
            <button className="delete" onClick={closeModal} />
          </header>
          <section className="modal-card-body">
            <div className="content">
              {children}
            </div>
          </section>
        </div>
      </div>
    );
}

const ParseProtoCluster = () => {
    return (
        <div className="column is-full">
            Not implemented yet
        </div>
    );
};

const querySubmission = (queryItems, setMatches) => {
    return (
        <div>
            <div className="panel-block">
                <div className="field has-addons">
                    <div className="control">
                        <button 
                            className="button is-link is-light" 
                            onClick={() => {
                                if (queryItems.length > 0) {
                                    console.log(queryItems)
                                    setMatches([]);
                                } else {
                                    toast.warn("Query is empty!");
                                }
                            }}
                            // disabled={isLoading || monomerGraphData.nodes.length === 0}
                            // disabled={isLoading || results.length === 0}
                        >
                            Match against database
                        </button>
                    </div>
                </div>
            </div>
            <div className="panel-block">
                <label className="checkbox">
                <input 
                    type="checkbox" 
                    defaultChecked={true}
                    // checked={matchOptions.matchMolecules} 
                    // onChange={() => handleOptionChange('matchMolecules')}
                />
                Match against parsed molecules
                </label>
            </div>
            <div className="panel-block">
                <label className="checkbox">
                <input 
                    type="checkbox" 
                    defaultChecked={false}
                    // checked={matchOptions.matchProtoClusters} 
                    // onChange={() => handleOptionChange('matchProtoClusters')}
                    disabled={true}
                />
                Match against parsed proto-clusters
                </label>
            </div>
        </div>
    )
};

const AlignmentTable = ({ data }) => {
    // Get the list of identifiers
    const identifiers = data.map(item => item.identifier);

    // Get the length of the sequences
    const sequenceLength = data[0].sequence.length;

    // Define colors based on starting characters
    const getColor = (char) => {
        if (char.startsWith('AA')) {
            return '#ed9ea8';
        } else if (char.startsWith('A')) {
            return '#209bef';
        } else if (char.startsWith('B')) {
            return '#fede57';
        } else if (char.startsWith('C')) {
            return '#47c774';
        } else if (char.startsWith('D')) {
            return '#ff3960';
        } else {
            return 'transparent';
        }
    };

    return (
        <div style={{ overflowX: 'auto', overflowY: 'auto' }}>
            <table style={{ borderCollapse: 'collapse' }}>
                <thead>
                    <tr>
                        <th>Identifier</th>
                        {/* Create table headers for each sequence */}
                        {Array.from(Array(sequenceLength).keys()).map(index => (
                            // <th key={index}>Module {index + 1}</th>
                            <th key={index} style={{ textAlign: 'center' }}>
                                {index + 1}
                            </th>
                        ))}
                        <th style={{ textAlign: 'left' }}>Bioactivity</th>
                    </tr>
                </thead>
                <tbody>
                    {/* Iterate over each item in data */}
                    {data.map(item => (
                        <tr key={item.identifier}>
                            <td style={{whiteSpace: "nowrap"}}>
                                {/* Render link if URL exists, else render plain text */}
                                {item.url ? (
                                    <a href={item.url} target="_blank" rel="noopener noreferrer">{item.identifier}</a>
                                ) : (
                                    item.identifier
                                )}
                            </td>
                            {/* Render each character in the sequence */}
                            {item.sequence.map((char, index) => (
                                <td key={index} style={{ backgroundColor: char === 'GAP' ? 'transparent' : getColor(char), textAlign: 'center' }}>
                                    {char === 'GAP' || char === '???' ? '' : char}
                                </td>
                            ))}
                            <td style={{ textAlign: "left",whiteSpace: 'nowrap' }}>
                                {item.bioactivity.length ? item.bioactivity.join(', ') : "N/A"}
                            </td>
                        </tr>
                    ))}
                </tbody>
            </table>
        </div>
    );
};

const ParseMolecule = () => {
    const [smiles, setSmiles] = useState("");
    const [isLoading, setIsLoading] = useState(false);

    const [inputSmiles, setInputSmiles] = useState("");
    const [inputSmilesSVG, setInputSmilesSVG] = useState("");
    const [linearizedSmiles, setLinearizedSmiles] = useState("");
    const [linearizedSmilesSVG, setLinearizedSmilesSVG] = useState("");

    const [results, setResults] = useState([]);
    const [selectedResult, setSelectedResult] = useState([]);

    const [modalActive, setModalActive] = useState(false);
    const [matches, setMatches] = useState([]);

    const [matchOptions, setMatchOptions] = useState({
        matchMolecules: true,
        matchProtoClusters: false
    });

    // when setResults is called, give every item in every seq a unique id
    const tagResults = (results) => {
        return results.map((seq) => {
            return seq.map((item) => {
                return { ...item, id: uuidv4() };
            });
        });
    };

    const handleOptionChange = (option) => {
        setMatchOptions({ ...matchOptions, [option]: !matchOptions[option] });
    };

    const toggleModal = () => {
        setModalActive(!modalActive);
    };

    const drawSmiles = async (smiles, setSvgString) => {
        try {
            const response = await fetch("/api/draw_smiles", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ "smiles": smiles })
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();
            
            // Unpack response.
            if (json.status === "success") {
                setSvgString(json.payload.svg_string);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
    
        } catch (error) {
            const msg = "Could not draw molecule!";
            toast.error(msg, { autoClose: true });
            console.error(error);
        };
    };

    const parseMolecule = async () => {
        setIsLoading(true);
        setMatches([]);
        setResults([]);
        setSelectedResult([]);

        try {
            const response = await fetch("/api/parse_retromol", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ "smiles": smiles })
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();
            
            if (json.status === "success") {
                toast.success(json.message);
                console.log(json.payload);
                if (json.payload.sequences) setResults(tagResults(json.payload.sequences));
                if (json.payload.input_smiles) setInputSmiles(json.payload.input_smiles);
                if (json.payload.linearized_smiles) setLinearizedSmiles(json.payload.linearized_smiles);
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

    // Load example -- sets the SMILES input to a default value
    const loadExample = () => {
        setSmiles("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O");
    };

    const clear = () => {
        setSmiles("");
        setMatches([]);
        setResults([]);
        setSelectedResult([]);
        setLinearizedSmiles("");
        setInputSmiles("");
        setInputSmilesSVG("");
        setLinearizedSmilesSVG("");
    };

    const findMatches = async () => {
        setIsLoading(true);
        setMatches([]);
        try {
            const response = await fetch("/api/find_matches", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ "queryItems": selectedResult })
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

    const matchSubmission = (queryItems, setMatches) => {

        return (
            <div>
                <div className="panel-block">
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
                </div>
                <div className="panel-block">
                    <label className="checkbox">
                    <input 
                        type="checkbox" 
                        defaultChecked={true}
                        disabled={true}
                        // checked={matchOptions.matchMolecules} 
                        // onChange={() => handleOptionChange('matchMolecules')}
                    />
                    Match against parsed molecules
                    </label>
                </div>
                <div className="panel-block">
                    <label className="checkbox">
                    <input 
                        type="checkbox" 
                        defaultChecked={false}
                        // checked={matchOptions.matchProtoClusters} 
                        // onChange={() => handleOptionChange('matchProtoClusters')}
                        disabled={true}
                    />
                    Match against parsed proto-clusters
                    </label>
                </div>
            </div>
        )
    };

    // on effect when input smiles and linearized smiles change they need to be drawn 
    useEffect(() => {
        if (inputSmiles) {
            drawSmiles(inputSmiles, setInputSmilesSVG);
            drawSmiles(linearizedSmiles, setLinearizedSmilesSVG);
        }
    }, [inputSmiles]);

    // useEffect(() => {
    //     if (linearizedSmiles) {
    //         drawSmiles(linearizedSmiles, setLinearizedSmilesSVG);
    //     }
    // }, [linearizedSmiles]);

    return (
        <div>
            {/* {isLoading &&  <div className="loader" style={{position: "fixed", top: "50%", left: "50%", transform: "translate(-50%, -50%)", backgroundColor: "rgba(1, 0, 0, 0.9)", zIndex: "1000"}} />} */}
            {isLoading && (
                <div className="loader-overlay" style={{position: "fixed", top: 0, left: 0, width: "100%", height: "100%", backgroundColor: "rgba(0, 0, 0, 0.5)", display: "flex", justifyContent: "center", alignItems: "center", zIndex: 1000}}>
                    <div className="loader" style={{position: "absolute", top: "50%", left: "50%", transform: "translate(-50%, -50%)"}} />
                </div>
            )}
            <div id="content-space" className="column is-full">
                <Modal closeModal={toggleModal} modalState={modalActive} title="Matching results">
                    <div>
                        {matches.length === 0 ? (
                            <p>Nothing to see here.</p>
                        ) : (
                            <AlignmentTable data={matches} />
                        )}
                    </div>
                </Modal>
                <div className="control" style={{border: "1px solid #dbdbdb", borderRadius: "5px", marginBottom: "10px"}}>
                    <div className="panel">
                        <div className="panel-heading">
                            <div className="title is-5">Input</div>
                        </div>
                        <div className="panel-block">
                            <div className="field" style={{width: "100%"}}>
                                <div 
                                    className="control"
                                >
                                    <input className="input" value={smiles} type="text" placeholder="Enter SMILES" onChange={(e) => setSmiles(e.target.value)} disabled={isLoading} />
                                </div>
                            </div>
                        </div>
                        <div className="panel-block">
                            <div className="field has-addons">
                                <div className="control">
                                    <button className="button is-link is-light" style={{marginRight: "5px", marginBottom: "5px"}} onClick={loadExample} disabled={isLoading}>Load example</button>
                                    <button className="button is-link is-light" style={{marginRight: "5px", marginBottom: "5px"}} onClick={parseMolecule} disabled={isLoading}>Parse molecule</button>
                                    <button className="button is-link is-light" style={{marginRight: "5px", marginBottom: "5px"}} onClick={clear} disabled={isLoading}>Clear</button>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                <div className="control" style={{border: "1px solid #dbdbdb", borderRadius: "5px", marginBottom: "10px"}}>
                    <div className="panel">
                        <div className="panel-heading">
                            <div className="title is-5">Results</div>
                        </div>
                        <div className="panel-block">
                            <div>
                                {results.length === 0 ? (
                                    <p>No results found.</p>
                                ) : (
                                    <div className="field is-grouped is-grouped-multiline" style={{ width: "100%" }}>
                                        {results.map((result, index) => (
                                            <div key={index} className="control">
                                                <div
                                                    className={`tags has-addons ${selectedResult === result ? 'selected' : ''}`}
                                                    style={{ marginRight: "5px", cursor: "pointer" }}
                                                    onClick={() => setSelectedResult(result)}
                                                >
                                                    <span className={`tag is-link ${selectedResult === result ? '' : 'is-light'}`}>
                                                        Result {index + 1}
                                                    </span>
                                                </div>
                                            </div>
                                        ))}
                                    </div>
                                )}
                            </div>
                        </div>
                        <div className="panel-block">
                            <div>
                                {results.length === 0 ? (
                                    <p>No results to draw.</p>
                                ) : (
                                    <div className="columns is-centered" style={{margin: "5px", verticalAlign: "middle", textAlign: "center", width: "100%"}}>
                                        <div className="column has-text-centered is-half">
                                            {/* <h3 style={{marginBottom: "10px"}}>
                                                Input molecule
                                            </h3> */}
                                            <div style={{width: "50%", margin: "0 auto"}}>
                                                <div dangerouslySetInnerHTML={{ __html: inputSmilesSVG }} />
                                            </div>
                                        </div>
                                        <div className="column has-text-centered is-half">  
                                            {/* <h3 style={{marginBottom: "10px"}}>
                                                Linearized backbone
                                            </h3> */}
                                            <div style={{width: "50%", margin: "0 auto"}}>
                                                <div dangerouslySetInnerHTML={{ __html: linearizedSmilesSVG }} />
                                            </div>
                                        </div>
                                    </div>
                                )}
                            </div>
                        </div>
                    </div>
                </div>
                {selectedResult.length > 0 && 
                    <MatchBuilder 
                        queryItems={selectedResult} 
                        setQueryItems={setSelectedResult} 
                        submissionElement={matchSubmission(selectedResult, setMatches)}
                    />}
            </div>
        </div>
    );
};

const QueryDatabase = () => {
    const [queryItems, setQueryItems] = useState([]);
    const [isLoading, setIsLoading] = useState(false);
    const [matches, setMatches] = useState([]);

    const [modalActive, setModalActive] = useState(false);

    const toggleModal = () => {
        setModalActive(!modalActive);
    };


    const queryDatabase = async () => {
        setIsLoading(true);
        setMatches([]);
        try {
            const response = await fetch("/api/query_database", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ "queryItems": queryItems })
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

    const querySubmission = (queryItems, setQueryMatches) => {
        return (
            <div>
                <div className="panel-block">
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
                </div>
                <div className="panel-block">
                    <label className="checkbox">
                    <input 
                        type="checkbox" 
                        defaultChecked={true}
                        disabled={true}
                        // checked={matchOptions.matchMolecules} 
                        // onChange={() => handleOptionChange('matchMolecules')}
                    />
                    Query against parsed molecules
                    </label>
                </div>
                <div className="panel-block">
                    <label className="checkbox">
                    <input 
                        type="checkbox" 
                        defaultChecked={false}
                        // checked={matchOptions.matchProtoClusters} 
                        // onChange={() => handleOptionChange('matchProtoClusters')}
                        disabled={true}
                    />
                    Query against parsed proto-clusters
                    </label>
                </div>
            </div>
        )
    };

    return (
        <div>
            {/* {isLoading &&  <div className="loader" style={{position: "fixed", top: "50%", left: "50%", transform: "translate(-50%, -50%)", backgroundColor: "rgba(1, 0, 0, 0.9)", zIndex: "1000"}} />} */}
            {isLoading && (
                <div className="loader-overlay" style={{position: "fixed", top: 0, left: 0, width: "100%", height: "100%", backgroundColor: "rgba(0, 0, 0, 0.5)", display: "flex", justifyContent: "center", alignItems: "center", zIndex: 1000}}>
                    <div className="loader" style={{position: "absolute", top: "50%", left: "50%", transform: "translate(-50%, -50%)"}} />
                </div>
            )}
            <div id="content-space" className="column is-full">
                <Modal closeModal={toggleModal} modalState={modalActive} title="Matching results">
                    <div>
                        {matches.length === 0 ? (
                            <p>Nothing to see here.</p>
                        ) : (
                            <AlignmentTable data={matches} />
                            // <p>Results</p>
                        )}
                    </div>
                </Modal>
                <div className="control" style={{border: "1px solid #dbdbdb", borderRadius: "5px", marginBottom: "10px"}}></div>
                <QueryBuilder 
                    queryItems={queryItems} 
                    setQueryItems={setQueryItems} 
                    submissionElement={querySubmission(queryItems, setMatches)}
                />
            </div>
        </div>
    );
}

const RetroMol = () => {
    const tabs = ["Overview", "Parse molecule", "Parse proto-cluster", "Query database"];
    const [selectedTab, setSelectedTab] = useState("Overview");

    return (
        <div>
            <div className="tabs is-boxed" style={{padding: "20px"}}>
                <ul>
                    {tabs.map((tab, index) => (
                        <li 
                            key={index}
                            className={selectedTab === tab ? "is-active" : ""}
                            onClick={() => setSelectedTab(tab)}

                        >
                            <a>{tab}</a>
                        </li>
                    ))}
                </ul>
            </div>
            <div className="container" style={{padding: "20px"}}>
                {selectedTab === "Overview" && <Overview />}
                {selectedTab === "Parse molecule" && <ParseMolecule />}
                {selectedTab === "Parse proto-cluster" && <ParseProtoCluster />}
                {selectedTab === "Query database" && <QueryDatabase />}
            </div>
        </div>
    );
};

export default RetroMol;