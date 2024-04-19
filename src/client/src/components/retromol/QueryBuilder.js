import React from "react";
import { v4 as uuidv4 } from "uuid";
import { toast } from "react-toastify";
import { DragDropContext, Droppable, Draggable } from "react-beautiful-dnd";

import QueryItem from "./QueryItem";

// Reorders the list of items after drag and drop.
const reorder = (list, startIndex, endIndex) => {
    const result = Array.from(list);
    const [removed] = result.splice(startIndex, 1);
    result.splice(endIndex, 0, removed);
    return result;
};

const QueryBuilder = ({ queryItems, setQueryItems, submissionElement }) => {  
    // Add a new wildcard item to the list of query items.
    const handleAddWildcardItem = () => {
        const newItem = [{
            id: uuidv4(),
            identifier: "any",
            properties: {
                size: null
            }
        }];
        setQueryItems([newItem, ...queryItems]);
    };

    // Add a new wildcard item to the query item list.
    const handleAddWildcardItemToList = (listIndex) => {
        const newItem = {
            id: uuidv4(),
            identifier: "any",
            properties: {
                size: null
            }
        };
        setQueryItems(queryItems.map((itemList, index) => index === listIndex ? [...itemList, newItem] : itemList));
    };

    // Delete an item from the list of query items.
    const handleDelete = index => {
        try {
            const updatedItems = queryItems.filter((_, i) => i !== index);
            setQueryItems(updatedItems);
        } catch (error) {
            console.error("Error while deleting:", error);
        };
    };

    // Delete an item from the query item list.
    const handleDeleteItemFromList = (listIndex, itemIndex) => {
        setQueryItems(queryItems.map((itemList, index) => index === listIndex ? itemList.filter((_, i) => i !== itemIndex) : itemList));
    };

    // Reorder the list of query items after drag and drop.
    const onDragEnd = result => {
        if (!result.destination) {
            return;
        };

        try {
            const updatedItems = reorder(
                queryItems,
                result.source.index,
                result.destination.index
            );
            setQueryItems(updatedItems);
        } catch (error) {
            toast.error("Error while reordering!");
            console.error("Error while reordering:", error);
        };
    };
  
    return (
        <div 
            className="control" 
            style={{
                border: "1px solid #dbdbdb", 
                borderRadius: "5px"
            }}
        >
            <div className="panel">
                <div className="panel-heading">
                    <div className="title is-5">
                        Query builder
                    </div>
                </div>

                <div className="panel-block">
                    <div className="field has-addons">
                        <div 
                            className="control" 
                            style={{ margin: "10px" }}
                        >
                            <button 
                                className="button is-link is-light" 
                                onClick={handleAddWildcardItem} 
                                style={{ marginRight: "5px" }}
                            >
                                Add module
                            </button>
                        </div>
                    </div>
                </div>

                {queryItems.length != 0 ? (
                    <div className="panel-block">
                        <div 
                            className="field" 
                            style={{ width: "100%" }}
                        >
                            <DragDropContext onDragEnd={onDragEnd}>
                                <Droppable droppableId="droppable">
                                    {(provided, _) => (
                                        <div 
                                            ref={provided.innerRef} 
                                            {...provided.droppableProps}
                                        >
                                            {queryItems.map((itemList, index) => (
                                                <Draggable 
                                                    key={itemList[0]?.id} 
                                                    draggableId={itemList[0]?.id} 
                                                    index={index}
                                                >
                                                    {(provided, _) => (
                                                        <div 
                                                            ref={provided.innerRef} 
                                                            {...provided.draggableProps} 
                                                            {...provided.dragHandleProps}
                                                        >
                                                            <QueryItem
                                                                itemList={itemList}
                                                                index={index}
                                                                queryItems={queryItems}
                                                                setQueryItems={setQueryItems}
                                                                handleDelete={handleDelete}
                                                                handleDeleteItemFromList={handleDeleteItemFromList}
                                                                handleAddItemToList={handleAddWildcardItemToList}
                                                            />
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
                        <span style={{ margin: "10px" }}>
                            No modules added yet.
                        </span>
                    </div>
                )}
                {submissionElement && submissionElement}
            </div>
        </div>
    );
};

export default QueryBuilder;