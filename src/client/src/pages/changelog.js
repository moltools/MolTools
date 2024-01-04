import React, { useEffect, useState } from "react";
import ReactMarkdown from "react-markdown";
import { FiRefreshCw } from "react-icons/fi";
import { toast } from "react-toastify";

const Changelog = ({url}) => {
    const [changelogContent, setChangelogContent] = useState("");
    const [isFetching, setIsFetching] = useState(true);

    const fetchChangelog = () => {
        setIsFetching(true);
        
        fetch(url)
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
                console.error("Error fetching changelog:", error);
                toast.error("Error fetching changelog!");
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
            <h2>
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
            </h2>
            {changelogContent.length ? (
                <ReactMarkdown>{changelogContent}</ReactMarkdown>
            ) : (
                <p>No changelog available.</p>
            )}
        </div>
    );
};

export default Changelog;