import React, { useState } from "react";
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

// =====================================================================================================================
// Comet component.
// =====================================================================================================================

/**
 * PARAS Component
 *
 * This component represents the "PARAS" section of the MolTools application.
 *
 * @returns {JSX.Element} The rendered PARAS component.
 */
const PARAS = () => {
    // Define state.
    const [isLoading, setIsLoading] = useState(false);
    const [fileContent, setFileContent] = useState(null);
    const [fileContentType, setFileContentType] = useState(null);
    const [prediction, setPrediction] = useState(null);

    // Set the file content.
    const handleFileSelected = (content) => {
        setFileContent(content);
    };

    // Set the file content type.
    const handleFileContentType = (contentType) => {
        setFileContentType(contentType);
    };

    // API call with file content and content type.
    const handlePredictPARAS = async () => {
        // Set is loading to true, which grays out model and deactivates buttons.
        setIsLoading(true);

        if (fileContent === null || fileContentType === null) {
            // Set is loading to false again.
            setIsLoading(false);
            toast.warn("Please select a file and its file type first!");
            return;
        };
    
        try {
            const response = await fetch("/api/predict_paras", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({
                    "file_content": fileContent,
                    "file_content_type": fileContentType,
                }),
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

    // Render the file content.
    return (
        <div className="widget-content">
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
            <button
                className="btn btn-primary"
                onClick={handlePredictPARAS}
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

export default PARAS;