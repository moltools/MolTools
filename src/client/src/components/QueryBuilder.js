import React, { useState } from "react";
import { v4 as uuidv4 } from "uuid";
import { toast } from "react-toastify";
import { DragDropContext, Droppable, Draggable } from "react-beautiful-dnd";

const reorder = (list, startIndex, endIndex) => {
    const result = Array.from(list);
    const [removed] = result.splice(startIndex, 1);
    result.splice(endIndex, 0, removed);
    return result;
};

const QueryBuilder = ({ queryItems, setQueryItems, submissionElement }) => {  
    const handleAddWildcardItem = () => {
        const newItem = {
            id: uuidv4(),
            identifier: "any",
            properties: {
                size: null
            }
        }
        setQueryItems([newItem, ...queryItems]);
    };

    const handleDelete = index => {
        try {
            const updatedItems = queryItems.filter((item, i) => i !== index);
            setQueryItems(updatedItems);
        } catch (error) {
            console.error("Error while deleting:", error);
        }
    };

    const onDragEnd = result => {
        if (!result.destination) {
            return;
        }
        try {
            const updatedItems = reorder(
                queryItems,
                result.source.index,
                result.destination.index
            );
            setQueryItems(updatedItems);
        } catch (error) {
            console.error("Error while reordering:", error);
        }
    };
  
    return (
        <div className="control" style={{border: "1px solid #dbdbdb", borderRadius: "5px", marginBottom: "10px"}}>
            <div className="panel">
                <div className="panel-heading">
                    <div className="title is-5">Query builder</div>
                </div>
                <div className="panel-block">
                    <div className="field has-addons">
                        <div className="control">
                            <button className="button is-link is-light" onClick={handleAddWildcardItem} style={{marginRight: "5px"}}>
                                Add module
                            </button>
                        </div>
                    </div>
                </div>
                    {queryItems.length != 0 ? (
                    <div className="panel-block">
                        <div className="field" style={{width: "100%"}}>
                            <DragDropContext onDragEnd={onDragEnd}>
                                <Droppable droppableId="droppable">
                                    {(provided, snapshot) => (
                                        <div ref={provided.innerRef} {...provided.droppableProps}>
                                            {queryItems.map((item, index) => (
                                                <Draggable key={item.id} draggableId={item.id} index={index}>
                                                    {(provided, snapshot) => (
                                                        <div ref={provided.innerRef} {...provided.draggableProps} {...provided.dragHandleProps}>
                                                            <div className="panel" style={{backgroundColor: "white", border: "1px solid #dbdbdb", borderRadius: "5px", marginBottom: "10px"}}>
                                                                <div className="panel-block">
                                                                    <div className="field has-addons">
                                                                        <div className="control" style={{verticalAlign: "center"}}>
                                                                            <span className="icon is-small">
                                                                                <i className="fas fa-grip-vertical" aria-hidden="true"></i>
                                                                            </span>
                                                                            {" Module " + (index + 1)}
                                                                        </div>
                                                                    </div>
                                                                </div>
                                                                <div className="panel-block">
                                                                    <div className="field has-addons">
                                                                        <div className="dropdown is-hoverable" style={{marginRight: "5px"}}>
                                                                            <div className="dropdown-trigger">
                                                                                <button className="button" aria-haspopup="true" aria-controls="dropdown-menu">
                                                                                <span style={{ textTransform: 'capitalize' }}>{item.identifier}</span>
                                                                                <span className="icon is-small">
                                                                                    <i className="fas fa-angle-down" aria-hidden="true"></i>
                                                                                </span>
                                                                                </button>
                                                                            </div>
                                                                            <div className="dropdown-menu" id="dropdown-menu" role="menu">
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
                                                                                            }
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
                                                                                            }
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
                                                                                            }
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
                                                                        {item.identifier === "polyketide" ? (
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
                                                                                        <a
                                                                                            className="dropdown-item"
                                                                                            onClick={() => {
                                                                                                const updatedItems = queryItems.map((queryItem, i) => {
                                                                                                    if (i === index) { 
                                                                                                        queryItem.properties.accessory_domains = null; 
                                                                                                        queryItem.properties.decoration_type = null;
                                                                                                    }
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
                                                                                                    if (i === index) { queryItem.properties.accessory_domains = []; }
                                                                                                    return queryItem;
                                                                                                });
                                                                                                setQueryItems(updatedItems);
                                                                                            }}
                                                                                        >
                                                                                            A
                                                                                        </a>
                                                                                        <a
                                                                                            className="dropdown-item"
                                                                                            onClick={() => {
                                                                                                const updatedItems = queryItems.map((queryItem, i) => {
                                                                                                    if (i === index) { queryItem.properties.accessory_domains = ["KR"]; }
                                                                                                    return queryItem;
                                                                                                });
                                                                                                setQueryItems(updatedItems);
                                                                                            }}
                                                                                        >
                                                                                            B
                                                                                        </a>
                                                                                        <a
                                                                                            className="dropdown-item"
                                                                                            onClick={() => {
                                                                                                const updatedItems = queryItems.map((queryItem, i) => {
                                                                                                    if (i === index) { queryItem.properties.accessory_domains = ["KR", "DH"]; }
                                                                                                    return queryItem;
                                                                                                });
                                                                                                setQueryItems(updatedItems);
                                                                                            }}
                                                                                        >
                                                                                            C
                                                                                        </a>
                                                                                        <a
                                                                                            className="dropdown-item"
                                                                                            onClick={() => {
                                                                                                const updatedItems = queryItems.map((queryItem, i) => {
                                                                                                    if (i === index) { queryItem.properties.accessory_domains = ["KR", "DH", "ER"]; }
                                                                                                    return queryItem;
                                                                                                });
                                                                                                setQueryItems(updatedItems);
                                                                                            }}
                                                                                        >
                                                                                            D
                                                                                        </a>
                                                                                    </div>
                                                                                </div>
                                                                            </div>
                                                                            
                                                                        ) : (
                                                                            <div></div>
                                                                        )}
                                                                        {item.identifier === "polyketide" ? (
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
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 1; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            1
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 2; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            2
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 3; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            3
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 4; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            4
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 5; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            5
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 6; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            6
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 7; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            7
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 8; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            8
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 9; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            9
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 10; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            10
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 11; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            11
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => {const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.decoration_type = 12; } return queryItem;}); setQueryItems(updatedItems);}}>
                                                                                            12
                                                                                        </a>
                                                                                    </div>
                                                                                </div>
                                                                            </div>
                                                                            
                                                                        ) : (
                                                                            <div></div>
                                                                        )}
                                                                        {item.identifier === "any" ? (
                                                                            <div>
                                                                                <div className="field has-addons" style={{marginRight: "5px"}}>
                                                                                    <div className="control">
                                                                                        <input
                                                                                            className="input"
                                                                                            type="text"
                                                                                            placeholder="Any size"
                                                                                            value={item.properties.size || ''}
                                                                                            onChange={e => {
                                                                                                const value = e.target.value.trim(); // trim any leading/trailing spaces
                                                                                                const updatedItems = queryItems.map((queryItem, i) => {
                                                                                                    if (i === index) {
                                                                                                        // Check if the value is not an empty string and is a valid integer
                                                                                                        const parsedValue = parseInt(value);
                                                                                                        if (!isNaN(parsedValue)) { // Check if the parsed value is a valid number
                                                                                                            queryItem.properties.size = parsedValue;
                                                                                                        } else {
                                                                                                            queryItem.properties.size = null; // Set to null if the value is not a valid integer
                                                                                                            toast.error("Please enter a valid integer as length!");
                                                                                                        }
                                                                                                    }
                                                                                                    return queryItem;
                                                                                                });
                                                                                                setQueryItems(updatedItems);
                                                                                            }}
                                                                                        />
                                                                                    </div>
                                                                                </div>
                                                                                {item.properties.size !== null && (
                                                                                    <button onClick={() => {
                                                                                        const updatedItems = queryItems.map((queryItem, i) => {
                                                                                            if (i === index) { 
                                                                                                queryItem.properties.size = null; // set to null
                                                                                            }
                                                                                            return queryItem;
                                                                                        });
                                                                                        setQueryItems(updatedItems);
                                                                                    }}>Set to undefined</button>
                                                                                )}
                                                                            </div>
                                                                        ) : (
                                                                            <div></div>
                                                                        )}
                                                                        {item.identifier === "polyketide" ? (
                                                                            <div></div>
                                                                        ) : (
                                                                            <div></div>
                                                                        )}
                                                                        {item.identifier === "peptide" ? (
                                                                            <div className="dropdown is-hoverable" style={{marginRight: "5px"}}>
                                                                                <div className="dropdown-trigger">
                                                                                    <button className="button" aria-haspopup="true" aria-controls="dropdown-menu">
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
                                                                                        <i className="fas fa-angle-down" aria-hidden="true"></i>
                                                                                    </span>
                                                                                    </button>
                                                                                </div>
                                                                                <div className="dropdown-menu" id="dropdown-menu" role="menu">
                                                                                    <div className="dropdown-content">
                                                                                        <a
                                                                                            className="dropdown-item"
                                                                                            onClick={() => {
                                                                                                const updatedItems = queryItems.map((queryItem, i) => {
                                                                                                    if (i === index) { 
                                                                                                        queryItem.properties.classification = null; 
                                                                                                        queryItem.properties.pubchem_cid = null;
                                                                                                    }
                                                                                                    return queryItem;
                                                                                                });
                                                                                                setQueryItems(updatedItems);
                                                                                            }}  
                                                                                        >
                                                                                            Any
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => { const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.classification = "polar and charged"; } return queryItem; }); setQueryItems(updatedItems); }}>
                                                                                            Polar and charged
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => { const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.classification = "small hydrophobic"; } return queryItem; }); setQueryItems(updatedItems); }}>
                                                                                            Small hydrophobic
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => { const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.classification = "small non-hydrophobic"; } return queryItem; }); setQueryItems(updatedItems); }}>
                                                                                            Small non-hydrophobic
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => { const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.classification = "tiny"; } return queryItem; }); setQueryItems(updatedItems); }}>
                                                                                            Tiny
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => { const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.classification = "bulky"; } return queryItem; }); setQueryItems(updatedItems); }}>
                                                                                            Bulky
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => { const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.classification = "cyclic aliphatic"; } return queryItem; }); setQueryItems(updatedItems); }}>
                                                                                            Cyclic aliphatic
                                                                                        </a>
                                                                                        <a className="dropdown-item" onClick={() => { const updatedItems = queryItems.map((queryItem, i) => { if (i === index) { queryItem.properties.classification = "cysteine-like"; } return queryItem; }); setQueryItems(updatedItems); }}>
                                                                                            Cysteine-like
                                                                                        </a>
                                                                                    </div>
                                                                                </div>
                                                                            </div>
                                                                        ) : (
                                                                            <div></div>
                                                                        )}
                                                                        <button
                                                                            className="button is-danger is-light"
                                                                            onClick={() => {handleDelete(index);}}
                                                                        >
                                                                            Delete
                                                                        </button>
                                                                    </div>
                                                                </div>
                                                            </div>
                                                        </div>
                                                    )}
                                                </Draggable>
                                            ))}
                                            {provided.placeholder}
                                        </div>
                                    )}
                                </Droppable>
                            </DragDropContext>  
                        </div>
                    </div>
                ) : (
                    <div className="panel-block">
                        No modules added yet.
                    </div>
                )}
                { submissionElement && (
                    submissionElement
                )}
            </div>
        </div>
    );
};

export default QueryBuilder;