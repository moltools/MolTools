import React, { useEffect, useState } from "react";
import ReactMarkdown from "react-markdown";
import { FiRefreshCw } from "react-icons/fi";
import { toast } from "react-toastify";

const Changelog = () => {
    const [changelogContent, setChangelogContent] = useState("");
    const [isFetching, setIsFetching] = useState(true);

    const fetchChangelog = () => {
        setIsFetching(true);

        fetch("https://raw.githubusercontent.com/moltools/RetroMol.GUI/main/CHANGELOG.md")
            .then((response) => {
                if (!response.ok) {
                    throw new Error("Network response was not ok");
                }
                return response.text();
            })
            .then((data) => {
                setChangelogContent(data);
            })
            .catch((error) => {
                console.error("Error fetching CHANGELOG.md:", error);
                toast.error("Error fetching CHANGELOG.md");
            })
            .finally(() => {
                setIsFetching(false);
            });
    };

    useEffect(() => {
        fetchChangelog();
    }, []);

    return (
        <div>
            <h1>
                Changelog
                <FiRefreshCw 
                    style={{
                        marginLeft: "0.5rem",
                        cursor: isFetching ? "not-allowed" : "pointer",
                    }}
                    onClick={isFetching ? null : fetchChangelog}
                    disabled={isFetching}
                    size={15}
                    className={isFetching ? "spin" : ""}
                />  
            </h1>
            {changelogContent.length ? (
                <ReactMarkdown>{changelogContent}</ReactMarkdown>
            ) : (
                <p>No changelog available.</p>
            )}
        </div>
    );
};

export default Changelog;