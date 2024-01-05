import React, { useState, useEffect } from "react";
import { toast } from 'react-toastify';
import ReactMarkdown from 'react-markdown';

// =====================================================================================================================
// Changelog component.
// =====================================================================================================================

/**
 * Changelog Component
 *
 * This component fetches and displays the changelog content from a given URL.
 *
 * @param {Object} props - The props for the Changelog component.
 * @param {string} props.url - The URL from which to fetch the changelog content.
 * @returns {JSX.Element} The rendered Changelog component displaying the changelog content.
 */
const Changelog = ({ url }) => {
    // State variable to store the fetched changelog content
    const [changelogContent, setChangelogContent] = useState("");

    /**
     * Fetch the changelog content from the provided URL.
     */
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
                console.error(error);
                toast.error("Error fetching changelog!");
            })
            .finally(() => {});
    };

    // Use the useEffect hook to fetch the changelog content on component mount
    useEffect(() => {
        fetchChangelog();
    }, []);

    return (
        <div className="widget-content">
            {/* Display the changelog content as Markdown */}
            {changelogContent.length ? (
                <ReactMarkdown>{changelogContent}</ReactMarkdown>
            ) : (
                <p>No changelog available!</p>
            )}
        </div>
    );
};

// Export the Changelog component as the default export.
export default Changelog;