import React, { useEffect, useState } from "react";
import ReactMarkdown from "react-markdown";
import { toast } from "react-toastify";

const Changelog = ({url}) => {
    const [changelogContent, setChangelogContent] = useState("");

    const fetchChangelog = () => {
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
            .finally(() => {});
    };

    useEffect(() => {
        fetchChangelog();
    }, []);

    return (
        <div className="changelog-container">
            {changelogContent.length ? (
                <ReactMarkdown>{changelogContent}</ReactMarkdown>
            ) : (
                <p>No changelog available.</p>
            )}
        </div>
    );
};

export default Changelog;