import React, { useRef } from 'react';

/**
 * FilePickerButton component.
 * 
 * @param {Object} props - The props for the FilePickerButton component.
 * @param {Function} props.onFileSelected - The function to call when a file is selected.
 * @param {string} props.className - The class name to apply to the button.
 */
function FilePickerButton(props) {
    // Destructure the props.
    const { onFileSelected, className } = props;

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
                console.error('Error reading file:', error);
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
            <button className={className} onClick={handleFileSelect}>Select File</button>
        </div>
    );
}

// Export the FilePickerButton component as the default export.
export default FilePickerButton;
