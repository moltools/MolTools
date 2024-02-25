import React, { useState, useRef, useEffect } from "react";
import { toast } from "react-toastify";
import { ForceGraph2D } from "react-force-graph";

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

const ParseMolecule = () => {
    const [smiles, setSmiles] = useState("");
    const [isLoading, setIsLoading] = useState(false);
    const [moleculeGraphData, setMoleculeGraphData] = useState({nodes: [], links: []});
    const [monomerGraphData, setMonomerGraphData] = useState({nodes: [], links: []});
    const [monomerGraphEdges, setMonomerGraphEdges] = useState([]);
    const [primarySequence, setPrimarySequence] = useState([]);
    const [selectedAtomIds, setSelectedAtomIds] = useState([]);

    const [modalActive, setModalActive] = useState(false);
    const [matches, setMatches] = useState([]);

    const toggleModal = () => {
        setModalActive(!modalActive);
    };

    const parseMolecule = async () => {
        setIsLoading(true);
        setMoleculeGraphData({nodes: [], links: []});
        setMonomerGraphData({nodes: [], links: []});
        setMonomerGraphEdges([]);
        setPrimarySequence([]);
        setSelectedAtomIds([]);
        setMatches([]);
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
                if (json.payload.molecule_graph_data) setMoleculeGraphData(json.payload.molecule_graph_data);
                if (json.payload.monomer_graph_data) setMonomerGraphData(json.payload.monomer_graph_data);
                if (json.payload.primary_seq_monomer_ids) setMonomerGraphEdges(json.payload.primary_seq_monomer_ids);
                if (json.payload.primary_seq) setPrimarySequence(json.payload.primary_seq);
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

    const embedSeq = async () => {
        setIsLoading(true);
        try {
            const response = await fetch("/api/embed_retromol", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ "primary_seq": primarySequence })
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            if (json.status === "success") {
                toast.success(json.message);
                if (json.payload.matches) setMatches(json.payload.matches);
                // console.log(json.payload.matches);
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

    // Load example -- sets the SMILES input to a default value
    const loadExample = () => {
        setSmiles("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O");
    };

    const clear = () => {
        setSmiles("");
        setMoleculeGraphData({nodes: [], links: []});
        setMonomerGraphData({nodes: [], links: []});
        setMonomerGraphEdges([]);
        setPrimarySequence([]);
        setSelectedAtomIds([]);
        setMatches([]);
    };

    const margin = 30;
    const rightContainer = document.getElementById("right-column");
    const widthRightContainer = rightContainer ? rightContainer.offsetWidth - margin : 0;
    const leftContainer = document.getElementById("left-column");
    const widthLeftContainer = leftContainer ? leftContainer.offsetWidth - margin: 0;
    const height = 400;

    const graphRef = useRef(null);
    const monomerGraphRef = useRef(null);

    return (
        <div className="column is-full">
            <Modal closeModal={toggleModal} modalState={modalActive} title="Embedding results">
                {primarySequence.length > 0 && (
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
                )}
            </Modal>
            <div className="field">
                <div 
                    className="control"
                >
                    <input className="input" value={smiles} type="text" placeholder="Enter SMILES" onChange={(e) => setSmiles(e.target.value)} disabled={isLoading} />
                </div>
            </div>
            <div className="control">
                <button className="button is-link is-light" style={{marginRight: "5px", marginBottom: "5px"}} onClick={loadExample} disabled={isLoading}>Load example</button>
                <button className="button is-link is-light" style={{marginRight: "5px", marginBottom: "5px"}} onClick={parseMolecule} disabled={isLoading}>Parse</button>
                <button className="button is-link is-light" style={{marginRight: "5px", marginBottom: "5px"}} onClick={clear} disabled={isLoading}>Clear</button>
                <button 
                    className="button is-link is-light" 
                    onClick={() => {
                        embedSeq();
                    }}
                    // enable when there is a monomergraph
                    disabled={isLoading || monomerGraphData.nodes.length === 0}
                >
                    Query database
                </button>
            </div>
            <div className="columns" style={{marginTop: "5px"}}>
                <div id="left-column" className="column has-text-centered">
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
                        width={widthLeftContainer}
                        height={height}
                        backgroundColor="#f5f5f5"
                        cooldownTicks={0}
                        onEngineStop={() => {
                            graphRef.current.zoomToFit();
                        }}
                        enableNodeDrag={false}
                        // enablePanInteraction={false}
                        // enableZoomInteraction={false}
                    />
                </div>
                <div id="right-column" className="column has-text-centered">
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
                        width={widthRightContainer}
                        height={height}
                        backgroundColor="#f5f5f5"
                        cooldownTicks={0}
                        onEngineStop={() => {
                            monomerGraphRef.current.zoomToFit();
                        }}
                        enableNodeDrag={false}
                        // enablePanInteraction={false}
                        // enableZoomInteraction={false}
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
    );
};

const ParseProtoCluster = () => {
    return (
        <div className="column is-full">
            Not implemented yet
        </div>
    );
};

const QueryDatabase = () => {
    return (
        <div className="column is-full">
            Not implemented yet
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