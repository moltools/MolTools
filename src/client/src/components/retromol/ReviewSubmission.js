import React, { useState } from "react";
import { toast } from "react-toastify";

const ReviewSubmission = ({ smiles, toggleModal }) => {
    const [msg, setMsg] = useState("");
    const [email, setEmail] = useState("");

    // Submit the input SMILES string for review.
    const submitForReview = async () => {
        if (smiles === "") { return; }

        const data = { 
            "smiles": smiles,
            "msg": msg,
            "email": email,
            "timestamp": new Date().toLocaleString()
        };

        try {
            const response = await fetch("/api/submit_for_review", {
                method: "POST",
                headers: {"Content-Type": "application/json"},
                body: JSON.stringify({ "data": data })
            });

            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();
            if (json.status === "success") {
                toast.success(json.message);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
        } catch (error) {
            toast.error(error.message);
        };
    };

    return (
        <div>
            <div className="field">
                <label className="label">SMILES</label>
                <div className="control">
                    <input
                        className="input"
                        type="text"
                        value={smiles}
                        disabled
                    />
                </div>
            </div>
            <div className="field">
                <label className="label">Message (required)</label>
                <div className="control">
                    <textarea 
                        className="textarea" 
                        placeholder="Enter a message for the reviewer..." 
                        onChange={(e) => setMsg(e.target.value)}
                    >
                    </textarea>
                </div>
            </div>
            <div className="field">
                <label className="label">Email (optional)</label>
                <div className="control">
                    <input 
                        className="input" 
                        type="email" 
                        placeholder="Enter your email address..." 
                        onChange={(e) => setEmail(e.target.value)}
                    />
                </div>
            </div>
            <div className="field">
                <div className="control">
                    <button 
                        className="button is-primary" 
                        onClick={() => {
                            if (msg === "") {
                                toast.error("Please enter a message for the reviewer!");
                                return;
                            };
                            submitForReview();
                            toggleModal();
                        }}
                    >
                        Submit for Review
                    </button>
                </div>
            </div>
        </div>
    );
};

export default ReviewSubmission;
