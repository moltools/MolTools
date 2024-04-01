import React from "react";
import { toast } from "react-toastify";

const IdentityPicker = ({ item, index, queryItems, setQueryItems }) => {
    return (
        <div 
            className="dropdown is-hoverable" 
            style={{ marginRight: "5px" }}
        >
            <div className="dropdown-trigger">
                <button 
                    className="button" 
                    aria-haspopup="true" 
                    aria-controls="dropdown-menu"
                >
                    <span style={{ textTransform: 'capitalize' }}>
                        {item.identifier}
                    </span>
                    <span className="icon is-small">
                        <i 
                            className="fas fa-angle-down" 
                            aria-hidden="true"
                        />
                    </span>
                </button>
            </div>
            <div 
                className="dropdown-menu" 
                id="dropdown-menu" 
                role="menu"
            >
                <div className="dropdown-content">
                <a 
                    className="dropdown-item"
                    onClick={() => {
                        const updatedItems = queryItems.map((queryItem, i) => {
                            if (i === index) { 
                                queryItem.identifier = "any"; 
                                queryItem.properties = {
                                    size: null
                                };
                            };
                            return queryItem;
                        });
                        setQueryItems(updatedItems);
                    }}
                >
                    Any
                </a>
                <a 
                    className="dropdown-item"
                    onClick={() => {
                        const updatedItems = queryItems.map((queryItem, i) => {
                            if (i === index) { 
                                queryItem.identifier = "polyketide"; 
                                queryItem.properties = { 
                                    accessory_domains: null,
                                    decoration_type: null
                                };
                            };
                            return queryItem;
                        });
                        setQueryItems(updatedItems);
                    }}
                >
                    Polyketide
                </a>
                <a 
                    className="dropdown-item"
                    onClick={() => {
                        const updatedItems = queryItems.map((queryItem, i) => {
                            if (i === index) { 
                                queryItem.identifier = "peptide"; 
                                queryItem.properties = {
                                    classification: null,
                                    pubchem_cid: null
                                };
                            };
                            return queryItem;
                        });
                        setQueryItems(updatedItems);
                    }}
                >
                    Peptide
                </a>
                </div>
            </div>
        </div>
    );
};

const PolyketideTypePicker = ({ item, index, queryItems, setQueryItems }) => {
    // General function to handle the accessory domain and decoration type. 
    const itemType = ({ name, accessory_domains, decoration_type }) => {
        return (
            <a
                className="dropdown-item"
                onClick={() => {
                    const updatedItems = queryItems.map((queryItem, i) => {
                        if (i === index) { 
                            queryItem.properties.accessory_domains = accessory_domains; 

                            if (decoration_type !== null) {
                                queryItem.properties.decoration_type = decoration_type;
                            };
                        }
                        return queryItem;
                    });
                    setQueryItems(updatedItems);
                }}
            >
                {name}
            </a>
        );
    };

    return (
        <div className="dropdown is-hoverable" style={{marginRight: "5px"}}>
            <div className="dropdown-trigger">
                <button className="button" aria-haspopup="true" aria-controls="dropdown-menu">
                <span style={{ textTransform: 'capitalize' }}>
                    {
                        item.properties.accessory_domains === null ? "Any" :
                        item.properties.accessory_domains.length === 0 ? "A" :
                        item.properties.accessory_domains.length === 1 ? "B" :
                        item.properties.accessory_domains.length === 2 ? "C" :
                        item.properties.accessory_domains.length === 3 ? "D" : "Any"
                    }
                </span>
                <span className="icon is-small">
                    <i className="fas fa-angle-down" aria-hidden="true"></i>
                </span>
                </button>
            </div>
            <div className="dropdown-menu" id="dropdown-menu" role="menu">
                <div className="dropdown-content">
                    {itemType({ name: "Any", accessory_domains: null, decoration_type: null })}
                    {itemType({ name: "A", accessory_domains: [], decoration_type: null })}
                    {itemType({ name: "B", accessory_domains: ["KR"], decoration_type: null })}
                    {itemType({ name: "C", accessory_domains: ["KR", "DH"], decoration_type: null })}
                    {itemType({ name: "D", accessory_domains: ["KR", "DH", "ER"], decoration_type: null })}
                </div>
            </div>
        </div>
    );
};

const PolyketideDecorationTypePicker = ({ item, index, queryItems, setQueryItems }) => {
    // General function to handle the accessory domain and decoration type.
    const itemType = ({ name, decoration_type }) => {
        return (
            <a 
                className="dropdown-item" 
                onClick={() => {
                        const updatedItems = queryItems.map((queryItem, i) => { 
                            if (i === index) { 
                                queryItem.properties.decoration_type = decoration_type; 
                            }; 
                            
                            return queryItem;
                        }); 
                        
                        setQueryItems(updatedItems);
                    }
                }
            >
                {name}
            </a>
        );
    };

    return (
        <div className="dropdown is-hoverable" style={{marginRight: "5px"}}>
            <div className="dropdown-trigger">
                <button className="button" aria-haspopup="true" aria-controls="dropdown-menu">
                <span>
                    {
                        item.properties.decoration_type === null ? "Any" :
                        item.properties.decoration_type
                    }
                </span>
                <span className="icon is-small">
                    <i className="fas fa-angle-down" aria-hidden="true"></i>
                </span>
                </button>
            </div>
            <div className="dropdown-menu" id="dropdown-menu" role="menu">
                <div className="dropdown-content" style={{height: "auto", maxHeight: "200px", overflowY: "auto"}}>
                    <a
                        className="dropdown-item"
                        onClick={() => {
                            const updatedItems = queryItems.map((queryItem, i) => {
                                if (i === index) {
                                    queryItem.properties.decoration_type = null;
                                }
                                return queryItem;
                            });
                            setQueryItems(updatedItems);
                        }}
                    >
                        Any
                    </a>  
                    {itemType({ name: 1, decoration_type: 1 })}
                    {itemType({ name: 2, decoration_type: 2 })}
                    {itemType({ name: 3, decoration_type: 3 })}
                    {itemType({ name: 4, decoration_type: 4 })}
                    {itemType({ name: 5, decoration_type: 5 })}
                    {itemType({ name: 6, decoration_type: 6 })}
                    {itemType({ name: 7, decoration_type: 7 })}
                    {itemType({ name: 8, decoration_type: 8 })}
                    {itemType({ name: 9, decoration_type: 9 })}
                    {itemType({ name: 10, decoration_type: 10 })}
                    {itemType({ name: 11, decoration_type: 11 })}
                    {itemType({ name: 12, decoration_type: 12 })}
                </div>
            </div>
        </div>
    );
};

const PolyketideAny = ({ item, index, queryItems, setQueryItems }) => {
    return (
        <div>
            <div 
                className="field has-addons" 
                style={{ marginRight: "5px" }}
            >
                <div className="control">
                    <input
                        className="input"
                        type="text"
                        placeholder="Any size"
                        value={item.properties.size || ""}
                        onChange={e => {
                            // Trim any leading/trailing spaces.
                            const value = e.target.value.trim();

                            const updatedItems = queryItems.map((queryItem, i) => {
                                if (i === index) {
                                    // Check if the value is not an empty string and is a valid integer
                                    const parsedValue = parseInt(value);

                                    if (!isNaN(parsedValue)) { 
                                        // Check if the parsed value is a valid number.
                                        queryItem.properties.size = parsedValue;
                                    } else {
                                        // Set to null if the value is not a valid integer.
                                        queryItem.properties.size = null; 
                                        toast.error("Please enter a valid integer as length!");
                                    };
                                };

                                return queryItem;
                            });

                            setQueryItems(updatedItems);
                        }}
                    />
                </div>
            </div>
            {item.properties.size !== null && (
                <button 
                    onClick={() => {
                        const updatedItems = queryItems.map((queryItem, i) => {
                            if (i === index) { 
                                // Set to null.
                                queryItem.properties.size = null; 
                            };

                            return queryItem;
                        });
                        setQueryItems(updatedItems);
                    }}
                >
                    Set to undefined
                </button>
            )}
        </div>
    );
};

const PeptideTypePicker = ({ item, index, queryItems, setQueryItems }) => {
    // General function to handle the accessory domain and decoration type.
    const itemType = ({ name, classification }) => {
        return (
            <a 
                className="dropdown-item" 
                onClick={() => { 
                    const updatedItems = queryItems.map((queryItem, i) => { 
                        if (i === index) { 
                            queryItem.properties.classification = "polar and charged"; 
                        };
                        
                        return queryItem; 
                    }); 
                    
                    setQueryItems(updatedItems); 
                }}
            >
                {name}
            </a>
        );
    };

    return (
        <div 
            className="dropdown is-hoverable" 
            style={{ marginRight: "5px" }}
        >
            <div className="dropdown-trigger">
                <button 
                    className="button" 
                    aria-haspopup="true" 
                    aria-controls="dropdown-menu"
                >
                <span>
                    {
                        item.properties.classification === null ? "Any" :
                        item.properties.classification === "polar and charged" ? "Polar and charged" :
                        item.properties.classification === "small hydrophobic" ? "Small hydrophobic" :
                        item.properties.classification === "small non-hydrophobic" ? "Small non-hydrophobic" :
                        item.properties.classification === "tiny" ? "Tiny" :
                        item.properties.classification === "bulky" ? "Bulky" :
                        item.properties.classification === "cyclic aliphatic" ? "Cyclic aliphatic" :
                        item.properties.classification === "cysteine-like" ? "Cysteine-like"
                        : "Any"
                    }
                </span>
                    <span className="icon is-small">
                        <i 
                            className="fas fa-angle-down" 
                            aria-hidden="true"
                        />
                    </span>
                </button>
            </div>
            <div 
                className="dropdown-menu" 
                id="dropdown-menu" 
                role="menu"
            >
                <div className="dropdown-content">
                    <a
                        className="dropdown-item"
                        onClick={() => {
                            const updatedItems = queryItems.map((queryItem, i) => {
                                if (i === index) { 
                                    queryItem.properties.classification = null; 
                                    queryItem.properties.pubchem_cid = null;
                                };

                                return queryItem;
                            });

                            setQueryItems(updatedItems);
                        }}  
                    >
                        Any
                    </a>
                    {itemType({ name: "Polar and charged", classification: "polar and charged" })}
                    {itemType({ name: "Small hydrophobic", classification: "small hydrophobic" })}
                    {itemType({ name: "Small non-hydrophobic", classification: "small non-hydrophobic" })}
                    {itemType({ name: "Tiny", classification: "tiny" })}
                    {itemType({ name: "Bulky", classification: "bulky" })}
                    {itemType({ name: "Cyclic aliphatic", classification: "cyclic aliphatic" })}
                    {itemType({ name: "Cysteine-like", classification: "cysteine-like" })}
                </div>
            </div>
        </div>
    );
};

const QueryItem = ({ item, index, queryItems, setQueryItems, handleDelete }) => {
    return (
        <div 
            className="panel" 
            style={{
                backgroundColor: "white", 
                border: "1px solid #dbdbdb", 
                borderRadius: "5px", 
                margin: "10px"
            }}
        >
            <div className="panel-block">
                <div className="field has-addons">
                    <div 
                        className="control" 
                        style={{ verticalAlign: "center" }}
                    >
                        <span className="icon is-small">
                            <i 
                                className="fas fa-grip-vertical" 
                                aria-hidden="true"
                            />
                        </span>
                        {" Module " + (index + 1)}
                    </div>
                </div>
            </div>
            <div className="panel-block">
                <div className="field has-addons">
                    <IdentityPicker 
                        item={item} 
                        index={index} 
                        queryItems={queryItems} 
                        setQueryItems={setQueryItems}
                    />
                    {item.identifier === "polyketide" ? (
                        <PolyketideTypePicker
                            item={item}
                            index={index}
                            queryItems={queryItems}
                            setQueryItems={setQueryItems}
                        />
                    ) : <div></div>}
                    {item.identifier === "polyketide" ? (
                        <PolyketideDecorationTypePicker
                            item={item}
                            index={index}
                            queryItems={queryItems}
                            setQueryItems={setQueryItems}
                        />
                    ) : <div></div>}
                    {item.identifier === "any" ? (
                        <PolyketideAny
                            item={item}
                            index={index}
                            queryItems={queryItems}
                            setQueryItems={setQueryItems}
                        />
                    ) : <div></div>}
                    {item.identifier === "peptide" ? (
                        <PeptideTypePicker
                            item={item}
                            index={index}
                            queryItems={queryItems}
                            setQueryItems={setQueryItems}
                        />
                    ) : <div></div>}
                    <button
                        className="button is-danger is-light"
                        onClick={() => {handleDelete(index);}}
                    >
                        Delete
                    </button>
                </div>
            </div>
        </div>
    );
};

export default QueryItem;