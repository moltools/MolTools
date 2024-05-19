import React from "react";

const AlignmentTable = ({ data }) => {

    // Define a function to handle downloading the results as a JSON file.
    const handleDownload = () => {
        const filename = "results.json";
        const dataJson = JSON.stringify(data, null, 4);
        const blob = new Blob([dataJson], { type: "application/json" });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = filename;
        link.click();
    };

    // Get the length of the sequences by checking the first item in the data array.
    const sequenceLength = data[0].sequence.length;

    // Define colors based on starting characters.
    const getColor = (char) => {
        if (char.startsWith("AA")) {
            return "#ed9ea8";
        } else if (char.startsWith("A")) {
            return "#209bef";
        } else if (char.startsWith("B")) {
            return "#fede57";
        } else if (char.startsWith("C")) {
            return "#47c774";
        } else if (char.startsWith("D")) {
            return "#ff3960";
        } else if (char === "???") {
            return "#ccc";
        } else {
            return "transparent";
        }
    };

    return (
        <div style={{ overflowX: "auto", overflowY: "auto" }}>
            <button
                className="button is-link is-light"
                style={{ marginBottom: "10px" }}
                onClick={handleDownload}
            >
                Download results
            </button>
            <table style={{ borderCollapse: "collapse" }}>
                <thead>
                    <tr>

                        {/* Create table header for identifier. */}
                        <th>
                            Identifier
                        </th>

                        {/* Create table headers for each sequence. */}
                        {Array.from(Array(sequenceLength).keys()).map(index => (
                            <th 
                                key={index} 
                                style={{ textAlign: "center" }}
                            >
                                {index + 1}
                            </th>
                        ))}

                        {/* Create table header for bioactivity. */}
                        <th style={{ textAlign: "left" }}>Bioactivity</th>

                    </tr>
                </thead>
                <tbody>

                    {/* Iterate over each item in data. */}
                    {data.map(item => (
                        <tr key={item.identifier}>

                            <td style={{ whiteSpace: "nowrap" }}>
                                {/* Render link if URL exists, else render plain text. */}
                                {item.url ? (
                                    <a 
                                        href={item.url} 
                                        target="_blank" 
                                        rel="noopener noreferrer"
                                    >
                                        {item.identifier}
                                    </a>
                                ) : (
                                    <a
                                        href={`https://pubchem.ncbi.nlm.nih.gov/#query=${item.identifier}`} 
                                        target="_blank" 
                                        rel="noopener noreferrer"
                                    >
                                        {item.identifier}
                                    </a>
                                )}
                            </td>

                            {/* Display each character in the sequence. */}
                            {item.sequence.map((char, index) => (
                                <td 
                                    key={index} 
                                    style={{ 
                                        backgroundColor: char === "GAP" ? "transparent" : getColor(char), 
                                        textAlign: "center"
                                    }}
                                >
                                    {char === "GAP" ? "" : char}
                                </td>
                            ))}

                            {/* Display bioactivity labels. */}
                            <td 
                                style={{ 
                                    textAlign: "left",
                                    whiteSpace: "nowrap" 
                                }}
                            >
                                {item.bioactivity.length ? item.bioactivity.join(", ") : "N/A"}
                            </td>
                        </tr>
                    ))}
                </tbody>
            </table>
        </div>
    );
};

export default AlignmentTable;