import React, { useState, useEffect } from "react";
import { FiCheckCircle, FiCircle } from "react-icons/fi";
import { toast } from "react-toastify";
import FilePickerButton from "../components/FilePicker";

// Define toggle component to set file content type.
const Toggle = ({ options, value, setValue }) => (
    <div className="toggle">
        {options.map((option) => (
            <label key={option}>
                <input
                    type="radio"
                    name="toggle"
                    value={option}
                    checked={value === option}
                    onChange={(event) => setValue(event.target.value)}
                />
                {option}
            </label>
        ))}
    </div>
);

// Define container that shows if a file was uploaded with checkmark or cross.
const FileUploaded = ({ fileContent }) => {
    if (fileContent === null) {
        return (
            <div>File uploaded: <FiCircle /></div>
        );
    } else {
        return (
            <div>File uploaded: <FiCheckCircle /></div>
        );
    };  
};

// SMILES input container, should trim whitespaces upon setting value.
const SMILESInput = ({ value, setValue }) => (
    <div>
        <label>SMILES (optional):</label>
        <input
            type="text"
            value={value}
            onChange={(event) => setValue(event.target.value.trim())}
        />
    </div>
);

// Checkbox with label.
const Checkbox = ({ label, value, setValue }) => (<label>
    <input
        type="checkbox"
        checked={value}
        onChange={(event) => setValue(event.target.checked)}
    />
        {label}
    </label>
    
);

// Sider.
const SliderWithDynamicMax = (props) => {
    // Unpack props.
    const { value, setValue, maxValue } = props;

    // const handleInputChange = (event) => {
    //     const newValue = parseInt(event.target.value, 10);
    //     setValue(newValue);
    //     if (!isNaN(newValue)) {
    //         setMaxValue(newValue);
    //     } else {
    //         setMaxValue(100);
    //     }
    // };

    return (
        <div>
            {/* <input
                type="number"
                value={value}
                onChange={handleInputChange}
                placeholder="Enter a value"
            /> */}
            <input
                type="range"
                min="0"
                max={maxValue}
                value={value}
                onChange={(event) => setValue(parseInt(event.target.value, 10))}
            />
        </div>
    );
};

// =====================================================================================================================
// Comet component.
// =====================================================================================================================

/**
 * PARASECT Component
 *
 * This component represents the "PARASECT" section of the MolTools application.
 *
 * @returns {JSX.Element} The rendered PARASECT component.
 */
const PARASECT = () => {
    // Define state.
    const [isLoading, setIsLoading] = useState(false);
    const [fileContent, setFileContent] = useState(null);
    const [fileContentType, setFileContentType] = useState(null);
    const [smiles, setSmiles] = useState(null);
    const [prediction, setPrediction] = useState(null);
    const [excludeStandardSubstrates, setExcludeStandardSubstrates] = useState(false);

    const [sliderValue, setSliderValue] = useState(3);
    const [sliderMaxValue, setSliderMaxValue] = useState(34);

    // Set the file content.
    const handleFileSelected = (content) => {
        setFileContent(content);
    };

    // Set the file content type.
    const handleFileContentType = (contentType) => {
        setFileContentType(contentType);
    };

    // If smiles has length, set max to + 1
    useEffect(() => { 
        if (smiles !== null && smiles.length > 0) {
            setSliderMaxValue(35);
        } else {
            setSliderMaxValue(34);
        };
    }, [smiles]);

    // API call with file content and content type.
    const handlePredictPARASECT = async () => {
        // Set is loading to true, which grays out model and deactivates buttons.
        setIsLoading(true);

        if (fileContent === null || fileContentType === null) {
            // Set is loading to false again.
            setIsLoading(false);
            toast.warn("Please select a file and its file type first!");
            return;
        };
    
        try {
            const response = await fetch("/api/predict_parasect", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({}),
            });
    
            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };
    
            const json = await response.json();
            
            // Unpack response.
            if (json.status === "success") {
                toast.success(json.message);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
    
            // Set is loading to false again.
            setIsLoading(false);

        } catch (error) {
            const msg = "Could not make prediction!";
            toast.error(msg, { autoClose: true });
            console.error(error);
    
            // Set is loading to false again.
            setIsLoading(false);
        };
    };

    return (
        <div className="widget-content">
            <div>
                <Toggle 
                    options={["GBK", "Fasta"]}
                    value={fileContentType}
                    setValue={handleFileContentType}
                />
                <FilePickerButton 
                    className="btn btn-primary" 
                    onFileSelected={handleFileSelected} 
                    disabled={isLoading}
                />
                <FileUploaded fileContent={fileContent} />
                <SMILESInput value={smiles} setValue={setSmiles} />
                <Checkbox 
                    label="Exclude standard substrates"
                    value={excludeStandardSubstrates}
                    setValue={setExcludeStandardSubstrates}
                />
                { excludeStandardSubstrates ? (
                    <div></div>
                ) : (
                    <SliderWithDynamicMax 
                        value={sliderValue}
                        setValue={setSliderValue}
                        maxValue={sliderMaxValue}
                    />
                )}
            </div>

            <button
                className="btn btn-primary"
                onClick={handlePredictPARASECT}
                disabled={isLoading}
            >
                Predict
            </button>
            <div>
                {prediction !== null && (
                    <div>{prediction}</div>
                ) || (
                    <div>No prediction made yet.</div>
                )}
            </div>
        </div>
    );
};

export default PARASECT;