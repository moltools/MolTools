import React, { useState, useRef, useEffect } from "react";
import { toast } from "react-toastify";
import { ForceGraph2D } from "react-force-graph";
import QueryDesigner from "../components/QueryDesigner";

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

function isListInListOfLists(list, listOfLists) {
    return listOfLists.some(subList => {
        return subList[0] === list[0] && subList[1] === list[1];
    });
};

const Modal = ({ children, closeModal, modalState, title }) => {
    if(!modalState) {
      return null;
    }
    
    return(
      <div className="modal is-active">
        <div className="modal-background" onClick={closeModal} />
        <div className="modal-card" style={{borderRadius: "10px"}}>
          <header className="modal-card-head">
            <p className="modal-card-title">{title}</p>
            <button className="delete" onClick={closeModal} />
          </header>
          <section className="modal-card-body">
            <div className="content">
              {children}
            </div>
          </section>
          {/* <footer className="modal-card-foot">
            <a className="button" onClick={closeModal}>Cancel</a>
          </footer> */}
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

const QueryBuilder = () => {
    const [queryItems, setQueryItems] = useState([]);
  
    const handleAddItem = () => {
      setQueryItems([...queryItems, 'New Item']);
    };
  
    const handleRemoveItem = (index) => {
      const updatedQueryItems = [...queryItems];
      updatedQueryItems.splice(index, 1);
      setQueryItems(updatedQueryItems);
    };
  
    const handleDragStart = (event, index) => {
      event.dataTransfer.setData('index', index.toString());
    };
  
    const handleDrop = (event) => {
      const indexFrom = parseInt(event.dataTransfer.getData('index'));
      const indexTo = parseInt(event.target.dataset.index);
  
      const updatedQueryItems = [...queryItems];
      const [removed] = updatedQueryItems.splice(indexFrom, 1);
      updatedQueryItems.splice(indexTo, 0, removed);
      setQueryItems(updatedQueryItems);
    };
  
    return (
      <div>
        <div className="panel">
          <div className="panel-heading">
            <button className="button is-primary" onClick={handleAddItem}>Add Item</button>
          </div>
        </div>
        <div className="field">
          {queryItems.map((item, index) => (
            <div
              key={index}
              data-index={index}
              className="query-item"
              draggable
              onDragStart={(event) => handleDragStart(event, index)}
              onDrop={handleDrop}
              onDragOver={(event) => event.preventDefault()}
            >
              <div className="panel">
                <div className="panel-block">
                  <div>{item}</div>
                  <button className="button is-danger" onClick={() => handleRemoveItem(index)}>Remove</button>
                </div>
              </div>
            </div>
          ))}
        </div>
      </div>
    );
};

const ParseMolecule = () => {
    const [smiles, setSmiles] = useState("");
    const [isLoading, setIsLoading] = useState(false);
    // const [moleculeGraphData, setMoleculeGraphData] = useState({nodes: [], links: []});
    // const [monomerGraphData, setMonomerGraphData] = useState({nodes: [], links: []});
    // const [monomerGraphEdges, setMonomerGraphEdges] = useState([]);
    // const [primarySequence, setPrimarySequence] = useState([]);
    // const [selectedAtomIds, setSelectedAtomIds] = useState([]);

    const [results, setResults] = useState([]);
    const [selectedResult, setSelectedResult] = useState(null);

    const [modalActive, setModalActive] = useState(false);
    const [matches, setMatches] = useState([]);

    const [matchOptions, setMatchOptions] = useState({
        matchMolecules: true,
        matchProtoClusters: false
    });

    const handleOptionChange = (option) => {
        setMatchOptions({ ...matchOptions, [option]: !matchOptions[option] });
    };

    const toggleModal = () => {
        setModalActive(!modalActive);
    };

    const parseMolecule = async () => {
        setIsLoading(true);
        // setMoleculeGraphData({nodes: [], links: []});
        // setMonomerGraphData({nodes: [], links: []});
        // setMonomerGraphEdges([]);
        // setPrimarySequence([]);
        // setSelectedAtomIds([]);
        setMatches([]);
        setResults([]);
        setSelectedResult(null);

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
                // if (json.payload.molecule_graph_data) setMoleculeGraphData(json.payload.molecule_graph_data);
                // if (json.payload.monomer_graph_data) setMonomerGraphData(json.payload.monomer_graph_data);
                // if (json.payload.primary_seq_monomer_ids) setMonomerGraphEdges(json.payload.primary_seq_monomer_ids);
                // if (json.payload.primary_seq) setPrimarySequence(json.payload.primary_seq);
                if (json.payload.sequences) setResults(json.payload.sequences);
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

    // const embedSeq = async () => {
    //     setIsLoading(true);
    //     try {
    //         const response = await fetch("/api/embed_retromol", {
    //             method: "POST",
    //             headers: {"Content-Type": "application/json"},
    //             body: JSON.stringify({ "primary_seq": primarySequence })
    //         });

    //         if (!response.ok) {
    //             throw new Error("Network response was not ok!");
    //         };

    //         const json = await response.json();

    //         if (json.status === "success") {
    //             toast.success(json.message);
    //             if (json.payload.matches) setMatches(json.payload.matches);
    //             // console.log(json.payload.matches);
    //             toggleModal();
    //         } else if (json.status === "warning") {
    //             toast.warn(json.message);
    //         } else if (json.status === "failure") {
    //             toast.error(json.message);
    //         };
    //     } catch (error) {
    //         toast.error(error.message);
    //     }
    //     setIsLoading(false);
    // };

    // Load example -- sets the SMILES input to a default value
    const loadExample = () => {
        setSmiles("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O");
    };

    const clear = () => {
        setSmiles("");
        // setMoleculeGraphData({nodes: [], links: []});
        // setMonomerGraphData({nodes: [], links: []});
        // setMonomerGraphEdges([]);
        // setPrimarySequence([]);
        // setSelectedAtomIds([]);
        setMatches([]);
        setResults([]);
        setSelectedResult(null);
    };

    // const columnsContainer = document.getElementById("content-space");
    // const width = columnsContainer ? (columnsContainer.clientWidth / 2) - 100 : 0;
    // const height = 400;

    const graphRef = useRef(null);
    const monomerGraphRef = useRef(null);

    return (
        <div id="content-space" className="column is-full">
            <Modal closeModal={toggleModal} modalState={modalActive} title="Embedding results">
                <div>
                    Nothing to see here...
                </div>
                {/* {primarySequence.length > 0 && (
                    <div>
                        <div className="columns">
                            <div className="column">
                                {matches.length > 0 ? (
                                    <table className="table is-fullwidth">
                                        <thead>
                                            <tr>
                                                <th>Match</th>
                                                <th>Name</th>
                                                <th>Distance</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            {matches.map((match, index) => (
                                                // console.log(match),
                                                <tr key={index}>
                                                    <td>{match.index + 1}</td>
                                                    {match.link == null ? (
                                                        <td>{match.label}</td>
                                                    ) : (
                                                        <td><a href={match.link} target="_blank" rel="noreferrer">{match.label}</a></td>
                                                    )}
                                                    <td>{match.distance}</td>
                                                </tr>
                                            ))}
                                        </tbody>
                                    </table>
                                ) : (
                                    <p>No matches found</p>
                                )}
                            </div>
                        </div>
                    </div>
                )} */}
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
                        <div class={`dropdown ${results.length > 0 ? "is-hoverable" : ""}`}>
                            <div class="dropdown-trigger">
                                <button class="button" aria-haspopup="true" aria-controls="dropdown-menu">
                                <span>{results.length > 0 ? "Select result" : "No results found"}</span>
                                <span class="icon is-small">
                                    {results.length > 0 ? <i class="fas fa-angle-down" aria-hidden="true"></i> : <i class="fas fa-times" aria-hidden="true"></i>}
                                </span>
                                </button>
                            </div>
                            <div class="dropdown-menu" id="dropdown-menu" role="menu">
                                <div class="dropdown-content">
                                {results.map((result, index) => (
                                    <a class="dropdown-item" key={index} onClick={() => setSelectedResult(result)}>
                                        {"Sequence " + (index + 1)}
                                    </a>
                                ))}
                                </div>
                            </div>
                        </div>
                    </div>
                    {selectedResult && (
                        <div className="panel-block">
                            <div className="content">
                                {Object.keys(selectedResult).length > 0 && (
                                    <div>
                                        <pre>{JSON.stringify(selectedResult, null, 2)}</pre>
                                    </div>
                                )}
                            </div>
                        </div>
                    )}
                    <div className="panel-block">
                        <div className="field has-addons">
                        <div className="control">
                            <button 
                                className="button is-link is-light" 
                                onClick={() => {
                                    if (selectedResult) {
                                        toast.warn("Not implemented yet!");
                                    } else {
                                        toast.warn("Select a result first!");
                                    }
                                }}
                                // disabled={isLoading || monomerGraphData.nodes.length === 0}
                                disabled={isLoading || results.length === 0}
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
                            checked={matchOptions.matchMolecules} 
                            onChange={() => handleOptionChange('matchMolecules')}
                        />
                        Match against parsed molecules
                        </label>
                    </div>
                    <div className="panel-block">
                        <label className="checkbox">
                        <input 
                            type="checkbox" 
                            checked={matchOptions.matchProtoClusters} 
                            onChange={() => handleOptionChange('matchProtoClusters')}
                            disabled={true}
                        />
                        Match against parsed proto-clusters
                        </label>
                    </div>
                </div>
            {/* </div>
                <div id="columns" className="columns">
                    <div id="left-column" className="column has-text-centered">
                        <div className="control">
                            <div className="panel" style={{marginBottom: "10px"}}>
                                <div id="left-panel-block" className="panel-block">
                                    <div className="field has-addons">
                                        <div className="control">
                                        <ForceGraph2D
                                            ref={graphRef}
                                            style={{position: "relative"}}
                                            graphData={moleculeGraphData}
                                            nodeRelSize={0.3}
                                            nodeLabel={(node) => `${node.name}`}
                                            nodeColor={(node) => {
                                                if (selectedAtomIds.includes(node.id)) {
                                                    return "orange";
                                                } else {
                                                    return node.color;
                                                }
                                            }}
                                            linkLabel={(link) => {
                                                if (link.bondtype === 1.0) {
                                                    return "Single";
                                                } else if (link.bondtype === 2.0) {
                                                    return "Double";
                                                } else if (link.bondtype === 3.0) {
                                                    return "Triple";
                                                } else {
                                                    return "Unknown";
                                                }
                                            }}
                                            linkWidth={6}
                                            width={
                                                document.getElementById("left-panel-block") ? 
                                                document.getElementById("left-panel-block").clientWidth - 25 : 
                                                0
                                            }
                                            // width={width}
                                            height = {
                                                document.getElementById("left-panel-block") ? 
                                                document.getElementById("left-panel-block").clientWidth - 25 : 
                                                0
                                            }
                                            // backgroundColor="#f5f5f5"
                                            backgroundColor="#fff"
                                            cooldownTicks={0}
                                            onEngineStop={() => {
                                                graphRef.current.zoomToFit();
                                            }}
                                            enableNodeDrag={false}
                                            enablePanInteraction={true}
                                            enableZoomInteraction={false}
                                        />
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                    <div id="right-column" className="column has-text-centered">
                        <div className="control">
                            <div className="panel" style={{marginBottom: "10px"}}>
                                <div id="right-panel-block" className="panel-block">
                                    <div className="field has-addons">
                                        <div className="control">
                                            <ForceGraph2D
                                                ref={monomerGraphRef}
                                                style={{position: "relative"}}
                                                graphData={monomerGraphData}
                                                nodeRelSize={0.5}
                                                nodeLabel={(node) => "Select"}
                                                onNodeClick={(node) => {
                                                    setSelectedAtomIds(node.atom_ids);
                                                }}
                                                nodeColor={(node) => node.color}
                                                linkColor={(link) => {
                                                    if (isListInListOfLists([link.source.id, link.target.id], monomerGraphEdges) || isListInListOfLists([link.target.id, link.source.id], monomerGraphEdges)) {
                                                        return "rgba(255, 0, 0, 0.5)";
                                                    };
                                                }}
                                                linkWidth={(link) => {
                                                    if (isListInListOfLists([link.source.id, link.target.id], monomerGraphEdges) || isListInListOfLists([link.target.id, link.source.id], monomerGraphEdges)) {
                                                        return 5;
                                                    } else {
                                                        return 2;
                                                    }
                                                }}
                                                // width={width}
                                                // height={height}
                                                width={
                                                    document.getElementById("right-panel-block") ? 
                                                    document.getElementById("right-panel-block").clientWidth - 25 : 
                                                    0
                                                }
                                                height = {
                                                    document.getElementById("right-panel-block") ? 
                                                    document.getElementById("right-panel-block").clientWidth - 25 : 
                                                    0
                                                }
                                                // backgroundColor="#f5f5f5"
                                                backgroundColor="#fff"
                                                cooldownTicks={0}
                                                onEngineStop={() => {
                                                    monomerGraphRef.current.zoomToFit();
                                                }}
                                                enableNodeDrag={false}
                                                enablePanInteraction={true}
                                                enableZoomInteraction={false}
                                                nodeCanvasObject={(node, ctx, globalScale) => {
                                                    const label = node.name;
                                                    const fontSize = 12/globalScale;
                                                    ctx.font = `${fontSize}px Sans-Serif`;
                                                    const textWidth = ctx.measureText(label).width;
                                                    const bckgDimensions = [textWidth, fontSize].map(n => n + fontSize * 0.2); // some padding

                                                    // ctx.fillStyle = "rgba(255, 255, 255, 0.5)";
                                                    ctx.fillStyle = "rgba(255, 140, 0, 0.5)";
                                                    ctx.fillRect(node.x - bckgDimensions[0] / 2, node.y - bckgDimensions[1] / 2, ...bckgDimensions);

                                                    ctx.textAlign = "center";
                                                    ctx.textBaseline = "middle";
                                                    ctx.fillStyle = "black";
                                                    ctx.fillText(label, node.x, node.y);
                                                }}
                                            />
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div> */}
                </div>
        </div>
    );
};

const QueryDatabase = () => {
    return (
        <div className="column is-full">
            {/* <QueryBuilder /> */}
            <QueryDesigner />
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