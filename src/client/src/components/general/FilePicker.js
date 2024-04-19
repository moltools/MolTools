import React, { useRef } from "react";

function FilePickerButton({ onFileSelected, className }) {
    const fileInputRef = useRef(null);

    const handleFileSelect = () => {
        fileInputRef.current.click();
    };

    const handleInputChange = async (event) => {
        const selectedFile = event.target.files[0];
        if (selectedFile) {
            try {
                const contents = await readFileContents(selectedFile);
                onFileSelected(contents);
            } catch (error) {
                console.error("Error reading file:", error);
            }
        }
    };

    const readFileContents = (file) => {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = (event) => {
                resolve(event.target.result);
            };
            reader.onerror = (error) => {
                reject(error);
            };
            reader.readAsText(file);
        });
    };

    return (
        <div>
            <input
                type="file"
                ref={fileInputRef}
                style={{ display: 'none' }}
                onChange={handleInputChange}
            />
            <button 
                className={className} 
                onClick={handleFileSelect}
            >
                Select File
            </button>
        </div>
    );
}

export default FilePickerButton;