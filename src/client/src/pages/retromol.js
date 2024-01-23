import React, { useState } from "react";
import { toast } from "react-toastify";

// =====================================================================================================================
// RetroMol component.
// =====================================================================================================================

/**
 * RetroMol Component
 *
 * This component represents the "RetroMol" section of the MolTools application, indicating that it is under construction.
 *
 * @returns {JSX.Element} The rendered RetroMol component displaying "Under construction."
 */
const RetroMol = () => {
    const [selectedWindow, setSelectedWindow] = useState(null);
    const [windows, setWindows] = useState({
        "Enter compound SMILES": <div>Option 1 Value</div>,
        "Upload antiSMASH JSON": <div>Option 2 Value</div>,
    });
    const handleDropdownChangeWindow = (event) => {
        setSelectedWindow(event.target.value);
    };

    const [selectedItem, setSelectedItem] = useState(null);
    const [items, setItems] = useState({});
    const handleDropdownChangeItem = (event) => {
        setSelectedItem(event.target.value);
    };

    const buttons = [
        <button 
            key="1"
            onClick={() => {toast.warning("Button not yet implemented!");}}
        >
            Parse
        </button>,
        <button 
            key="2"
            onClick={() => {toast.warning("Button not yet implemented!");}}
        >
            GitHub
        </button>,
    ];

    return (
        <div className="retromol">
            <div className="retromol sidebar">
                {buttons}
            </div>
            <div className="retromol workspace">
                <div className="retromol workspace vertical-container">
                    <div className="header">Parse input</div>
                    <select value={selectedWindow} onChange={handleDropdownChangeWindow}>
                        <option value={null}>No input type selected</option>
                        {Object.keys(windows).map((key) => (
                            <option key={key} value={key}>
                                {key}
                            </option>
                        ))}
                    </select>
                    {selectedWindow ? windows[selectedWindow] : null}
                </div>
                <div className="retromol workspace vertical-container">
                    <div className="header">Query result</div>
                    <select 
                        value={selectedItem} onChange={handleDropdownChangeItem}
                        disabled={Object.keys(items).length === 0}
                    >
                        <option value={null}>No item selected</option>
                        {Object.keys(items).map((key) => (
                            <option key={key} value={key}>
                                {key}
                            </option>
                        ))}
                    </select>
                    <div className="content">
                        {selectedItem ? items[selectedItem] : null}
                    </div>
                </div>
            </div>
        </div>
    );
};

// Export the RetroMol component as the default export.
export default RetroMol;