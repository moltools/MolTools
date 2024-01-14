import React, { useState } from "react";

/**
 * TextInput component represents an input field with dynamic behavior based on props.
 * 
 * This component represents an input field with dynamic behavior based on props.
 * 
 * @param {Object} props - The props for the TextInput component.
 * @param {string} props.className - The CSS class for the input field.
 * @param {boolean} props.locked - Determines if the input field is locked.
 * @param {boolean} props.active - Indicates if the input field is active.
 * @param {string} props.value - The current value of the input field.
 * @param {string} props.error - Any error message associated with the input.
 * @param {string} props.label - The label to display as a placeholder.
 * @param {function} props.setValue - A function to handle value changes.
 * @param {Array} props.predicted - An array of predicted values for the input.
 * @returns {JSX.Element} - The rendered TextInput component.
 */
const TextInput = (props) => {
    // Destructure props into individual variables.
    const { 
        className,
        locked,
        active: propsActive, 
        value: propsValue, 
        error: propsError, 
        label: propsLabel, 
        setValue, 
        predicted 
    } = props;

    // Define state variables using the useState hook
    const [active, setActive] = useState(locked ? propsActive : false);
    const [value, setValueState] = useState(propsValue || "");
    const [error, setError] = useState(propsError || "");
    const label = propsLabel || "Label";

    /**
     * Handles changes in the input value.
     * @param {Object} event - The input change event.
     */
    const changeValue = (event) => {
        const newValue = event.target.value;
        setValue(newValue);
        setValueState(newValue);
        setError("");
    };

    // Determine the CSS class for the input field based on conditions
    const fieldClassName = `field ${(locked ? active : active || value) && "active"} ${locked && !active && "locked"}`;

    return (
        <div className={`${fieldClassName} ${className}`}>
            {active && value && predicted && predicted.includes(value) && <p className="predicted">{predicted}</p>}
            <input
                id={1}
                type="text"
                value={value}
                placeholder={label}
                onChange={changeValue}
                onFocus={() => !locked && setActive(true)}
                onBlur={() => !locked && setActive(false)}
            />
            <label htmlFor={1} className={error && "error"}>
                {error || label}
            </label>
        </div>
    );
};

// Export the TextInput component as the default export.
export default TextInput;