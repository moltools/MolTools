import React from "react";

const Overview = () => {
    return (
        <div> 
            <h1 className="title">
                RetroMol
            </h1>
            <p>
                RetroMol is a tool for the generation of retrosynthetic trees for a given molecule, with the purpose
                of identifying potential precursors for the molecule. The scope of the tool is to provide a user-friendly
                interface for the discovery of retrobiosynthesis of modular natural products. Modular natural products
                currently supported by RetroMol include polyketides, non-ribosomal peptides, and polyketide-nonribosomal
                peptide hybrids.
            </p>
            <br />
            <div>
                You can use RetroMol from two starting inputs:
                <ol style={{ paddingLeft: "40px" }}>
                    <li>A molecule represented by a SMILES string</li>
                    <li>An output file from <i>antiSMASH</i> containing mined proto-clusters</li>
                </ol>
            </div>
            <br />
            <p>
                Select the appropriate tab to get started. After parsing the molecule or proto-cluster, you can query
                the database to retrieve molecules and clusters with similar architectures.
            </p>
        </div>
    );
};

export default Overview;