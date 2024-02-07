import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";
import { 
    BsArrowRepeat,
    BsBrushFill, 
    BsChevronDoubleLeft, 
    BsChevronDoubleRight, 
    BsCircleHalf, 
    BsDropletFill, 
    BsEyeFill, 
    BsFillCloudDownloadFill,
    BsFillCloudUploadFill, 
    BsFillDatabaseFill,
    BsFillLightningFill, 
    BsGithub, 
    BsGlobe2
} from "react-icons/bs";

// =====================================================================================================================
// Example SDF string.
// =====================================================================================================================

/**
 * Example SDF string for Penicillin G.
 * @type {string}
 */
const exampleSdfString = `5904
  -OEChem-12182318453D

 41 43  0     1  0  0  0  0  0999 V2000
   -0.8019    1.2308    0.5170 S   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2842   -2.5451   -1.2026 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.3517    1.0760   -0.8170 O   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1157   -0.6970    0.5961 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1598   -2.0405    1.2167 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4781   -0.7369    0.3018 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5677   -1.3807   -0.3375 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3100   -0.4177    1.1279 C   0  0  1  0  0  0  0  0  0  0  0  0
   -2.5670    1.6679    0.1134 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1671    0.3360   -0.3862 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.6142   -1.6325    0.4990 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.9181   -1.8261   -0.3007 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2193    2.2224    1.3871 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5771    2.7313   -0.9863 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6296    0.1576   -0.1297 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8730   -1.6107    0.1024 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9020   -1.2655   -0.9563 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8679   -0.1949   -0.5225 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5670    1.1373   -0.7687 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0457   -0.5556    0.1171 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4619    2.1290   -0.3670 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9405    0.4362    0.5190 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6486    1.7786    0.2769 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5167   -0.4852    2.1999 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9961    0.1926   -1.4616 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4509   -2.4585    1.2025 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2655    1.4761    2.1882 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6724    3.0921    1.7709 H   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2446    2.5572    1.1944 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0666    3.6467   -0.6644 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0842    2.3779   -1.8995 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6008    3.0119   -1.2561 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4359   -1.0113   -1.2762 H   0  0  0  0  0  0  0  0  0  0  0  0
   -6.3150    0.9762   -0.6611 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.4453   -2.1891   -1.1950 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.4229   -0.9646   -1.8960 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.6443    1.4206   -1.2675 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.2818   -1.5979    0.3120 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.2343    3.1743   -0.5548 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.8642    0.1634    1.0209 H   0  0  0  0  0  0  0  0  0  0  0  0
    6.3451    2.5508    0.5901 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  8  1  0  0  0  0
  1  9  1  0  0  0  0
  2 12  2  0  0  0  0
  3 15  1  0  0  0  0
  3 34  1  0  0  0  0
  4 15  2  0  0  0  0
  5 16  2  0  0  0  0
  6  8  1  0  0  0  0
  6 10  1  0  0  0  0
  6 12  1  0  0  0  0
  7 11  1  0  0  0  0
  7 16  1  0  0  0  0
  7 33  1  0  0  0  0
  8 11  1  0  0  0  0
  8 24  1  0  0  0  0
  9 10  1  0  0  0  0
  9 13  1  0  0  0  0
  9 14  1  0  0  0  0
 10 15  1  0  0  0  0
 10 25  1  0  0  0  0
 11 12  1  0  0  0  0
 11 26  1  0  0  0  0
 13 27  1  0  0  0  0
 13 28  1  0  0  0  0
 13 29  1  0  0  0  0
 14 30  1  0  0  0  0
 14 31  1  0  0  0  0
 14 32  1  0  0  0  0
 16 17  1  0  0  0  0
 17 18  1  0  0  0  0
 17 35  1  0  0  0  0
 17 36  1  0  0  0  0
 18 19  2  0  0  0  0
 18 20  1  0  0  0  0
 19 21  1  0  0  0  0
 19 37  1  0  0  0  0
 20 22  2  0  0  0  0
 20 38  1  0  0  0  0
 21 23  2  0  0  0  0
 21 39  1  0  0  0  0
 22 23  1  0  0  0  0
 22 40  1  0  0  0  0
 23 41  1  0  0  0  0
M  END
$$$$`;

// =====================================================================================================================
// Sidebar buttons.
// =====================================================================================================================

/**
 * SidebarButton Component
 *
 * This component represents a button element typically used within the sidebar of the MolTools application. It provides
 * options for user interaction and displays an icon along with a title.
 *
 * @param {Object} props - The props for the SidebarButton component.
 * @param {boolean} props.disabled - Indicates whether the button is disabled.
 * @param {JSX.Element} props.icon - The JSX element representing an icon to display within the button.
 * @param {string} props.title - The title or label to display alongside the icon.
 * @param {function} props.onClick - A function to be called when the button is clicked.
 * @returns {JSX.Element} The rendered SidebarButton component with icon, title, and click functionality.
 */
const SidebarButton = (props) => {
    // Destructure props into individual variables.
    const { disabled, icon, title, onClick } = props;

    return (
        <button className="cinemol-sidebar-button" onClick={onClick} disabled={disabled}>
            {icon}
            <span>{title}</span>
        </button>
    );
};

/**
 * SidebarCounter Component
 *
 * This component represents a counter with increment and decrement buttons, typically used within the sidebar of the
 * MolTools application. It provides options for user interaction, displays an icon, title, current value, and allows
 * users to increment or decrement the value.
 *
 * @param {Object} props - The props for the SidebarCounter component.
 * @param {boolean} props.disabled - Indicates whether the counter buttons are disabled.
 * @param {JSX.Element} props.icon - The JSX element representing an icon to display.
 * @param {string} props.title - The title or label to display alongside the icon.
 * @param {number} props.value - The current value of the counter.
 * @param {function} props.onIncrement - A function to be called when the increment button is clicked.
 * @param {function} props.onDecrement - A function to be called when the decrement button is clicked.
 * @returns {JSX.Element} The rendered SidebarCounter component with icon, title, value, and increment/decrement buttons.
 */
const SidebarCounter = (props) => {
    // Destructure props into individual variables.
    const { disabled, icon, title, value, onIncrement, onDecrement } = props;

    return (
        <div className="cinemol-sidebar-counter-container">
            <div>
                {icon}
                <span className="cinemol-sidebar-counter-title">{title}: {value}</span>
            </div>
            <div>
                <button className="cinemol-sidebar-counter-button left" onClick={onDecrement} disabled={disabled}>{"<"}</button>
                <button className="cinemol-sidebar-counter-button right" onClick={onIncrement} disabled={disabled}>{">"}</button>
            </div>
        </div>
    );
};

/**
 * RotationCounter Component
 * 
 * This component represents a counter with increment and decrement buttons, typically used within the sidebar of the
 * MolTools application. It provides options for user interaction, displays an icon, title, current value, and allows
 * users to increment or decrement the value.
 * 
 * @param {Object} props - The props for the RotationCounter component.
 * @param {string} props.title - The title or label to display alongside the icon.
 * @param {boolean} props.isLoading - Indicates whether the counter buttons are disabled.
 * @param {number} props.rotation - The current value of the counter.
 * @param {function} props.setRotation - A function to be called when the increment button is clicked.
 * @returns {JSX.Element} The rendered RotationCounter component with icon, title, value, and increment/decrement buttons.
 */
const RotationCounter = (props) => {
    // Destructure props into individual variables.
    const { title, isLoading, rotation, setRotation } = props;

    // Maximum rotation is 2 * pi - pi / 12.
    const maxRotation = 2 * Math.PI - Math.PI / 12;
    
    /**
     * Increment the rotation.
     */
    const handleIncrement = () => {
        const incrementedValue = rotation + Math.PI / 12;
        setRotation(() =>
            incrementedValue <= maxRotation
            ? Math.round(incrementedValue * 10) / 10
            : Math.round((incrementedValue % maxRotation) * 10) / 10
        );
    };
    
    /**
     * Decrement the rotation.
     */
    const handleDecrement = () => {
        const decrementedValue = rotation - Math.PI / 12;
        setRotation(() =>
            decrementedValue >= 0
            ? Math.round(decrementedValue * 10) / 10
            : Math.round((maxRotation + decrementedValue) * 10) / 10
        );
    };
  
    return (
        <SidebarCounter
            disabled={isLoading}
            icon={<BsGlobe2 />}
            title={title}
            value={rotation}
            onIncrement={handleIncrement}
            onDecrement={handleDecrement}
        />
    );
};  

// =====================================================================================================================
// Component for Sidebar.
// =====================================================================================================================

/**
 * Sidebar Component
 *
 * This component represents the sidebar in the "CineMol" section of the MolTools application. It provides options to
 * toggle the sidebar, display version information, and render buttons for user interaction.
 *
 * @param {Object} props - The props for the Sidebar component.
 * @param {boolean} props.isOpen - Indicates whether the sidebar is open or closed.
 * @param {function} props.toggleSidebar - A function to toggle the visibility of the sidebar.
 * @param {string} props.version - The version information to display.
 * @param {JSX.Element} props.buttons - Buttons or components for user interaction within the sidebar.
 * @returns {JSX.Element} The rendered Sidebar component with toggle, version, and buttons.
 */
const Sidebar = (props) => {
    // Destructure props into individual variables.
    const { isOpen, toggleSidebar, version, buttons } = props;

    return (
        <div className={"cinemol-sidebar" + (isOpen ? " open" : "")}>
            {/* Button to toggle the sidebar */}
            <button onClick={toggleSidebar} className="cinemol-sidebar-toggle">
                {isOpen ? 
                    <BsChevronDoubleLeft className="cinemol-sidebar-toggle-icon" />  
                    : 
                    <BsChevronDoubleRight className="cinemol-sidebar-toggle-icon" /> 
                }
            </button>
            <div className="cinemol-sidebar-content">
                {/* Display the version information */}
                <div className="cinemol-sidebar-version">Version {version}</div>
                {/* Render buttons or components for user interaction */}
                {buttons}
            </div>
        </div>
    );
};

// =====================================================================================================================
// CineMol component.
// =====================================================================================================================

/**
 * CineMol component widget
 * 
 * This component represents the "CineMol" section of the MolTools application. It provides options for user interaction
 * and displays a 3D molecular model.
 * 
 * State Variables:
 * - version: The version of the CineMol component.
 * - mode: Dark or light background of the molecular model.
 * - svgString: SVG representation of the molecular model.
 * - sidebarOpen: Determines if the sidebar is collapsed or not.
 * - isLoading: App grays out when loading.
 * - initialRender: Initial render of the molecular model.
 * - sdfString: SDF string of the molecular model.
 * - style: Style of the molecular model.
 * - look: Look of the molecular model.
 * - includeHydrogens: Include hydrogens in the molecular model.
 * - resolution: Resolution of the molecular model.
 * - rotationX: Rotation over the x-axis of the molecular model in radians.
 * - rotationY: Rotation over the y-axis of the molecular model in radians.
 * - rotationZ: Rotation over the z-axis of the molecular model in radians.
 * - viewBox: View box of the molecular model.
 * 
 * Side Effects:
 * - When the component is rendered, the molecular model is drawn.
 * - When the component is rendered, the version of the CineMol component is retrieved.
 * - When the component is rendered, the sidebar is collapsed.
 * 
 * Usage:
 * - Include this component in your React application to create a molecule drawing widget.
 * - Requires the Sidebar, SidebarButton, and SidebarCounter components.
 * 
 * @returns {JSX.Element} The rendered CineMol component.
 */
const CineMol = () => {
    // General state variables.
    const [version, setVersion] = useState("0.0.0");                    // Version of the CineMol component.
    const [mode, setMode] = useState("dark");                           // Dark or light background of the molecular model.
    const [svgString, setSvgString] = useState("");                     // SVG representation of the molecular model.
    const [sidebarOpen, setSideBarOpen] = useState(true);               // Determines if the sidebar is collapsed or not.
    const [isLoading, setIsLoading] = useState(false);                  // App grays out when loading.
    const [initialRender, setInitialRender] = useState(true);           // Initial render of the molecular model.

    // Button state variables.
    const [sdfString, setSdfString] = useState("");                     // SDF string of the molecular model.
    const [style, setStyle] = useState("space-filling");                // Style of the molecular model.
    const [look, setLook] = useState("glossy");                         // Look of the molecular model.
    const [includeHydrogens, setIncludeHydrogens] = useState(false);    // Include hydrogens in the molecular model.
    const [resolution, setResolution] = useState(50);                   // Resolution of the molecular model.
    const [rotationX, setRotationX] = useState(0.0);                    // Rotation over the x-axis of the molecular model in radians.
    const [rotationY, setRotationY] = useState(0.0);                    // Rotation over the y-axis of the molecular model in radians.
    const [rotationZ, setRotationZ] = useState(0.0);                    // Rotation over the z-axis of the molecular model in radians.
    const [viewBox, setViewBox] = useState(undefined);                  // View box of the molecular model.

    // Downloads the SVG representation of the molecular model.
    const handleDownloadSvgString = () => {
        if (svgString.length > 0) {
            const blob = new Blob([svgString], { type: "image/svg+xml" });
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement("a");
            a.href = url;
            a.download = "model.svg";
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            document.body.removeChild(a);
        } else {
            toast.error("No SVG model to download!", { autoClose: true });
        }
    };

    // Uploads an SDF file and sets the SDF string.
    const handleUploadSdfFile = () => {
        const fileInput = document.createElement("input");
        fileInput.type = "file";
        fileInput.accept = ".sdf";
        fileInput.onchange = event => {
            const file = event.target.files[0];

            // Check if the file size is greater than 1 MB (1,000,000 bytes).
            if (file.size > 1000000) {
                toast.error("File size exceeds 1 MB limit!", { autoClose: true });
                return;
            };

            const reader = new FileReader();
            reader.onload = event => {
                const sdfString = event.target.result;
                toast.success("SDF file uploaded!", { autoClose: true });
                setSdfString(sdfString);
            };
            reader.readAsText(file);
        };
        fileInput.click();
    };

    // Get the version of the CineMol component.
    const handleGetVersion = async () => {
        try {
            const response = await fetch("/api/fetch_cinemol_version");
            
            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };

            const json = await response.json();

            // Unpack response.
            if (json.status === "success") {
                setVersion(json.payload.version);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };

        } catch (error) {
            const msg = "Could not retrieve version!";
            toast.error(msg, { autoClose: true });
            console.error(error);
        };
    };

    // Draw the molecular model. 
    const handleDrawModel = async () => {
        // Set is loading to true, which grays out model and deactivates buttons.
        setIsLoading(true);

        // Get dimensions parent container.
        const container = document.querySelector(".cinemol-viewer-container");
        var sidebarWidth = 0;
        if (sidebarOpen) {
            const sidebar = document.querySelector(".cinemol-sidebar");
            sidebarWidth = sidebar.clientWidth;
        };
        const width = container.clientWidth - sidebarWidth - 5; // Add small margin to prevent horizontal scrollbar.
        const height = container.clientHeight - 5; // Add small margin to prevent vertical scrollbar.

        if (sdfString.length === 0) {
            // Set is loading to false again.
            setIsLoading(false);
            return;
        };
    
        try {
            const response = await fetch("/api/draw_model", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({
                    "sdf_string": sdfString,
                    "style": style,
                    "look": look,
                    "include_hydrogens": includeHydrogens,
                    "resolution": resolution,
                    "rotation_x": rotationX,
                    "rotation_y": rotationY,
                    "rotation_z": rotationZ,
                    "view_box": viewBox,
                    "width": width,
                    "height": height,
                }),
            });
    
            if (!response.ok) {
                throw new Error("Network response was not ok!");
            };
    
            const json = await response.json();
            
            // Unpack response.
            if (json.status === "success") {
                setSvgString(json.payload.svg_string);
                setViewBox(json.payload.view_box);
            } else if (json.status === "warning") {
                toast.warn(json.message);
            } else if (json.status === "failure") {
                toast.error(json.message);
            };
    
            // Set is loading to false again.
            setIsLoading(false);

        } catch (error) {
            const msg = "Could not parse SDF string!";
            toast.error(msg, { autoClose: true });
            console.error(error);
    
            // Set is loading to false again.
            setIsLoading(false);
        };
    };

    /**
     * Toggle style.
     */
    const handleToggleStyle = () => {
        if        (style === "space-filling")   { setStyle("ball-and-stick");
        } else if (style === "ball-and-stick")  { setStyle("tube");
        } else if (style === "tube")            { setStyle("wireframe");
        } else if (style === "wireframe")       { setStyle("space-filling");
        } else                                  { setStyle("ball-and-stick"); }
    };

    /**
     * Toggle look.
     */
    const handleToggleLook = () => {
        if        (look === "glossy")   { setLook("cartoon");
        } else if (look === "cartoon")  { setLook("glossy");
        } else                          { setLook("cartoon"); }
    };
    
    /**
     * Toggle background style.
     */
    const handleToggleMode = () => {
        if (mode === "dark") { setMode("light"); } else { setMode("dark"); };
    };

    /**
     * Reset the state variables.
     */
    const handleReset = () => {
        setMode("dark");
        setSvgString("");
        setSdfString("");
        setStyle("space-filling");
        setLook("glossy");
        setIncludeHydrogens(false);
        setResolution(50);
        setRotationX(0.0);
        setRotationY(0.0);
        setRotationZ(0.0);
        setViewBox(undefined);
    };

    /**
     * Get the version of the CineMol component and skip the first render. When any of the following state variables 
     * change, the component will re-render and the useEffect hook will run again.
     */ 
    useEffect(() => {
        if (initialRender) {
            handleGetVersion();
            setInitialRender(false);
        } else {
            handleDrawModel();
        }
    }, [sdfString, style, look, includeHydrogens, resolution, rotationX, rotationY, rotationZ]);

    // Sidebar buttons.
    const buttons = [
        <SidebarButton key={1} disabled={isLoading} icon={<BsFillCloudUploadFill />} title="Upload SDF file" onClick={handleUploadSdfFile} />,
        <SidebarButton key={2} disabled={isLoading} icon={<BsFillCloudDownloadFill />} title="Download model as SVG" onClick={handleDownloadSvgString} />,
        <SidebarButton key={3} disabled={isLoading} icon={<BsFillDatabaseFill />} title="Load example: penicillin G" onClick={ () => setSdfString(exampleSdfString) } />,
        <SidebarButton key={4} disabled={isLoading} icon={<BsBrushFill />} title={`Toggle style: ${style}`}  onClick={handleToggleStyle} />,
        <SidebarButton key={5} disabled={isLoading} icon={<BsEyeFill />} title={`Toggle look: ${look}`}  onClick={handleToggleLook} />,
        <SidebarButton key={6} disabled={isLoading} icon={<BsDropletFill />} title={`Toggle hydrogens: ${includeHydrogens ? "true" : "false"}`} onClick={ () => setIncludeHydrogens(!includeHydrogens) } />,
        <SidebarButton key={7} disabled={isLoading} icon={<BsCircleHalf />} title={`Toggle background: ${mode}`} onClick={handleToggleMode} />,
        <SidebarCounter key={9} disabled={isLoading} icon={<BsFillLightningFill />} title="Resolution" value={resolution} 
            onIncrement={ () => { if (resolution < 100) { setResolution(resolution + 5) } } } 
            onDecrement={ () => { if (resolution >= 30) { setResolution(resolution - 5) } } }
        />,
        <RotationCounter key={10} title="Rotation X (rad)" isLoading={isLoading} rotation={rotationX} setRotation={setRotationX}/>,
        <RotationCounter key={11} title="Rotation Y (rad)" isLoading={isLoading} rotation={rotationY} setRotation={setRotationY}/>,
        <RotationCounter key={12} title="Rotation Z (rad)" isLoading={isLoading} rotation={rotationZ} setRotation={setRotationZ}/>,
        <SidebarButton key={13} disabled={isLoading} icon={<BsArrowRepeat />} title="Reset" onClick={handleReset} />,
        <SidebarButton key={8} disabled={isLoading} icon={<BsGithub />} title="Open GitHub issues" onClick={ () => window.open("https://github.com/moltools/CineMol/issues", "_blank") } />,
    ];

    return (
        <div className="cinemol">
            {/* Sidebar component for user interaction */}
            <Sidebar isOpen={sidebarOpen} toggleSidebar={() => setSideBarOpen(!sidebarOpen)} version={version} buttons={buttons} />
            <div className={`cinemol-viewer-container ${mode}`}>
                {/* Display the molecular model */}
                <div className={`cinemol-viewer ${isLoading ? 'grayed-out' : ''}`} dangerouslySetInnerHTML={{ __html: svgString }} />
            </div>
        </div>
    );
};

// Export the CineMol component as the default export.
export default CineMol;