import React from "react";

const SmilesInput = ({ 
    smiles, 
    setSmiles, 
    parseMolecule, 
    handleRefresh 
}) => {
    // Load an example SMILES string for erythromycin.
    const loadExample = () => {
        setSmiles("CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O");
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
                        <div className="control">
                            <input 
                                className="input" 
                                value={smiles} 
                                type="text" 
                                placeholder="Enter SMILES" 
                                onChange={(e) => setSmiles(e.target.value)} 
                            />
                        </div>
                    </div>
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
                                onClick={loadExample} 
                            >
                                Load example
                            </button>
                            <button 
                                className="button is-link is-light" 
                                style={{
                                    marginRight: "5px", 
                                    marginBottom: "5px"
                                }} 
                                onClick={parseMolecule} 
                            >
                                Parse molecule
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

export default SmilesInput;