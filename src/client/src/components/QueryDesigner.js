
import React, { useState, useEffect } from "react";
import { v4 as uuidv4 } from "uuid";
import { DragDropContext, Droppable, Draggable } from "react-beautiful-dnd";
import { toast } from "react-toastify";

const reorder = (list, startIndex, endIndex) => {
    const result = Array.from(list);
    const [removed] = result.splice(startIndex, 1);
    result.splice(endIndex, 0, removed);
    return result;
};

const grid = 8;

// const getItemStyle = (isDragging, draggableStyle, itemName, itemType) => ({
//     userSelect: "none",
//     display: "flex", // Use flexbox for centering
//     justifyContent: "center", // Center horizontally
//     alignItems: "center", // Center vertically
//     minWidth: 100,
//     maxWidth: 100,
//     minHeight: 50,
//     maxHeight: 50,
//     margin: `0 ${grid / 2}px ${grid}px 0`, // Add margin on all sides
//     background:
//         // isDragging ? "lightgreen"
//         itemName === "PKUnit" ? "#D3D3D3"
//         : itemName === "NRPUnit" ? "#708090"
//         : "gray",
//     borderRadius: "4px", // Add border radius for a rounded appearance
//     boxShadow: isDragging
//         ? "0 4px 8px rgba(0, 0, 0, 0.2)" // Add shadow when dragging
//         : "0 2px 4px rgba(0, 0, 0, 0.1)", // Slight shadow when not dragging
//     ...draggableStyle,
// });


// const getListStyle = (isDraggingOver, itemCount) => ({
//     // background: isDraggingOver ? "#F5C900" : "#928E95",
//     minHeight: "60px",
//     padding: grid,
//     display: "flex",
//     justifyContent: "center", // Center items horizontally
//     alignItems: "center", // Center items vertically
//     overflowX: "auto",
//     width: "100%",
//     transition: "background-color 0.3s ease-in-out",
//     overflowY: "hidden",
//     borderRadius: "3px", // Add border radius for a rounded appearance
//     height: itemCount * (50 + grid) + 2 * grid, // Height based on number of items
// });

const QueryDesigner = () => {
    const [items, setItems] = useState([]);

    const [pksConfigValid, setPksConfigValid] = useState(true);

    const isPksConfigValid = () => {
        if (pksOptions.KRAccessory === false && pksOptions.DHAccessory === false && pksOptions.ERAccessory === false) {
            return true;
        } else if (pksOptions.KRAccessory === true && pksOptions.DHAccessory === false && pksOptions.ERAccessory === false) {
            return true;
        } else if (pksOptions.KRAccessory === true && pksOptions.DHAccessory === true && pksOptions.ERAccessory === false) {
            return true;
        } else if (pksOptions.KRAccessory === true && pksOptions.DHAccessory === true && pksOptions.ERAccessory === true) {
            return true;
        } else {
            return false;
        }
    };

    const [pksOptions, setPksOptions] = useState({
        KRAccessory: false,
        DHAccessory: false,
        ERAccessory: false
    });

    const handlePkOptionsChange = (option) => {
        setPksOptions({ ...pksOptions, [option]: !pksOptions[option] });
    };

    const onDragEnd = result => {
        if (!result.destination) {
            return;
        }
        try {
            const updatedItems = reorder(
                items,
                result.source.index,
                result.destination.index
            );
            setItems(updatedItems);
        } catch (error) {
            console.error("Error while reordering:", error);
        }
    };

    const addNewItem = (name, type, subtype) => {
        try {
            const newItem = {
                id: uuidv4(), // Generate a unique ID
                // content: name === `PKUnit` ? `PK ${type}${subtype}` : `NRP ${type === "AlphaAminoAcid" ? "α" : "β"}`,
                content: name === `PKUnit` ? `${type}` : "AA",
                name: name,
                type: type,
                subtype: subtype
            };
            setItems([...items, newItem]);
        } catch (error) {
            console.error("Error while adding new item:", error);
        }
    };

    const handleDelete = index => {
        try {
            const updatedItems = items.filter((item, i) => i !== index);
            setItems(updatedItems);
        } catch (error) {
            console.error("Error while deleting:", error);
        }
    };

    const handleClear = () => {
        try {
            setItems([]);
        } catch (error) {
            console.error("Error while clearing:", error);
        }
    };

    // every time tags PKS are toggle, check if the configuration is valid
    useEffect(() => {
        setPksConfigValid(isPksConfigValid());
    }, [pksOptions]);

    const Toolbar = () => {
        return (
            <div className="column is-full">
                <div className="control">
                    <div className="panel" style={{marginBottom: "10px"}}>
                        {/* <h3>Toolbar</h3> */}
                        <div className="panel-block">
                            <div className="field has-addons">
                                <div className="control">
                                    <button 
                                        className="button is-link is-light" 
                                        onClick = {() => {
                                            handleClear();
                                        }}
                                    >
                                        Start new query
                                    </button>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                {/* <div className="control">
                    <div className="panel" style={{marginBottom: "10px"}}>
                        <div className="panel-block">
                            <div className="field has-addons">
                                <div className="control">
                                    <button 
                                        className="button is-link is-light" 
                                        style={{marginRight: "5px"}}
                                        onClick = {() => {
                                            toast.error("Not implemented yet!");
                                        }}
                                    >
                                        Add start token
                                    </button>
                                    <button
                                        className="button is-link is-light"
                                        onClick = {() => {
                                            toast.error("Not implemented yet!");
                                        }}
                                    >
                                        Add end token
                                    </button>
                                </div>
                            </div>
                        </div>
                    </div>
                </div> */}
                {/* <div className="control">
                    <div className="panel" style={{marginBottom: "10px"}}>
                        <div className="panel-block">
                            <div className="field has-addons">
                                <div className="control">
                                    <button 
                                        className="button is-link is-light" 
                                        style={{marginRight: "5px"}}
                                        onClick = {() => {
                                            toast.error("Not implemented yet!");
                                        }}
                                    >
                                        Add wild card
                                    </button>
                                    <div className="control">
                                        <input 
                                            className="input" 
                                            type="text" 
                                            placeholder="Min length"
                                            style={{width: "110px"}}
                                        />
                                    </div>
                                    <div className="control">
                                        <input 
                                            className="input" 
                                            type="text" 
                                            placeholder="Max length"
                                            style={{width: "110px"}}
                                        />
                                    </div>
                                </div>
                            </div>
                        </div>
                        <div className="panel-block">
                            <label className="checkbox">
                            <input 
                                type="checkbox" 
                                checked={true}
                            />
                            Could be any amino acid
                            </label>
                        </div>
                        <div className="panel-block">
                            <label className="checkbox">
                            <input 
                                type="checkbox" 
                                checked={true}
                            />
                            Could be any polyketide synthase domain
                            </label>
                        </div>
                    </div>
                </div> */}
                <div className="control">
                    <div className="panel" style={{marginBottom: "10px"}}>
                        <div className="panel-block">
                            <div className="field has-addons">
                                <div className="control">
                                    <button 
                                        className="button is-link is-light" 
                                        onClick = {() => {
                                            addNewItem("NRPUnit", "AlphaAminoAcid", null);
                                        }}
                                    >
                                        Add adenylation domain
                                    </button>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                <div className="control">
                    <div className="panel" style={{marginBottom: "10px"}}>
                        <div className="panel-block">
                            <div className="field has-addons">
                                <div className="control">
                                    <button 
                                        className="button is-link is-light" 
                                        disabled = {!pksConfigValid}
                                        onClick = {() => {
                                            if (pksOptions.KRAccessory === false && pksOptions.DHAccessory === false && pksOptions.ERAccessory === false) {
                                                addNewItem("PKUnit", "A", 1);
                                            } else if (pksOptions.KRAccessory === true && pksOptions.DHAccessory === false && pksOptions.ERAccessory === false) {
                                                addNewItem("PKUnit", "B", 1);
                                            } else if (pksOptions.KRAccessory === true && pksOptions.DHAccessory === true && pksOptions.ERAccessory === false) {
                                                addNewItem("PKUnit", "C", 1);
                                            } else if (pksOptions.KRAccessory === true && pksOptions.DHAccessory === true && pksOptions.ERAccessory === true) {
                                                addNewItem("PKUnit", "D", 1);
                                            } else {
                                                toast.error("Invalid polyketide synthase domain configuration!");
                                            }
                                        }}
                                    >
                                        Add polyketide synthase domain
                                    </button>
                                </div>
                            </div>
                        </div>
                        <div className="panel-block">
                            <label className="checkbox">
                            <input 
                                type="checkbox" 
                                checked={pksOptions.KRAccessory}
                                onChange={() => handlePkOptionsChange('KRAccessory')}
                            />
                            Has ketoreductase accessory domain
                            </label>
                        </div>
                        <div className="panel-block">
                            <label className="checkbox">
                            <input 
                                type="checkbox" 
                                checked={pksOptions.DHAccessory}
                                onChange={() => handlePkOptionsChange('DHAccessory')}
                            />
                            Has dehydratase accessory domain
                            </label>
                        </div>
                        <div className="panel-block">
                            <label className="checkbox">
                            <input
                                type="checkbox" 
                                checked={pksOptions.ERAccessory}
                                onChange={() => handlePkOptionsChange('ERAccessory')}
                            />
                            Has enoylreductase accessory domain
                            </label>
                        </div>
                    </div>
                </div>
            </div>
        );
    };

    return (
        <div className="container">
            <Toolbar />
            <div className="column is-full">
                <div className="control">
                    <div className="panel" style={{marginBottom: "10px"}}>
                        <div className="panel-block">
                            <div className="field has-addons">
                                <div className="column is-full">
                                <div className="control">
                                    <DragDropContext onDragEnd={onDragEnd}>
                                        <Droppable droppableId="droppable">
                                            {(provided, snapshot) => (
                                                <div ref={provided.innerRef} {...provided.droppableProps}>
                                                    {items.map((item, index) => (
                                                        <Draggable key={item.id} draggableId={item.id} index={index}>
                                                            {(provided, snapshot) => (
                                                                <div ref={provided.innerRef} {...provided.draggableProps} {...provided.dragHandleProps}>
                                                                    <div className="panel" style={{marginBottom: "10px"}}>
                                                                        <div className="panel-block">
                                                                            <div className="field has-addons">
                                                                                <div className="control">
                                                                                    <button
                                                                                        className="button is-danger is-light"
                                                                                        onClick={() => {
                                                                                            handleDelete(index);
                                                                                        }}
                                                                                    >
                                                                                        Delete
                                                                                    </button>
                                                                                </div>
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
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
};

export default QueryDesigner;