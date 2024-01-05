import React, { useState, useEffect } from "react";
import { toast } from "react-toastify";
import { 
    BsChevronDoubleLeft, BsChevronDoubleRight, 
    BsFillCloudUploadFill, BsFillCloudDownloadFill,
    BsCircleHalf, BsFillDatabaseFill,
    BsBrushFill, BsEyeFill, BsDropletFill, 
    BsGithub, BsFillLightningFill, BsGlobe2
} from "react-icons/bs";

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

class UploadButton extends React.Component {
    handleUploadSdfFile = () => {
        const fileInput = document.createElement("input");
        fileInput.type = "file";
        fileInput.accept = ".sdf";
        fileInput.onchange = event => {
            const file = event.target.files[0];
            const reader = new FileReader();
            reader.onload = event => {
                const sdfString = event.target.result;
                toast.success("SDF file uploaded!", { autoClose: true });
                this.props.setSdfString(sdfString);
            };
            reader.readAsText(file);
        };
        fileInput.click();
    };

    render() {
        return (
            <button className="square-button" onClick={this.handleUploadSdfFile} disabled={this.props.isLoading}>
                <BsFillCloudUploadFill />
                <span>Upload SDF file</span>
            </button>
        );
    };
}

class LoadExampleButton extends React.Component {
    handleLoadExample = () => {
        this.props.setSdfString(exampleSdfString);
    };

    render() {
        return (
            <button className="square-button" onClick={this.handleLoadExample} disabled={this.props.isLoading}>
                <BsFillDatabaseFill />
                <span>Load Penicillin G example</span>
            </button>
        );
    };
};

const DownloadSvgButton = (props) => {
    return (
        <button className="square-button" onClick={props.downloadSvgString} disabled={props.isLoading}>
            <BsFillCloudDownloadFill />
            <span>Download model as SVG</span>
        </button>
    );
};

class ToggleStyleButton extends React.Component {
    toggleStyle = () => {
        if (this.props.style === "Space-filling") {
            this.props.setStyle("Ball-and-stick");
        } else if (this.props.style === "Ball-and-stick") {
            this.props.setStyle("Tube");
        } else if (this.props.style === "Tube") {
            this.props.setStyle("Wireframe");
        } else if (this.props.style === "Wireframe") {
            this.props.setStyle("Space-filling");
        } else {
            toast.error("Unknown style!", { autoClose: true });
            this.props.setStyle("Ball-and-stick");
        }
    };

    render() {
        return (
            <button className="square-button" onClick={this.toggleStyle} disabled={this.props.isLoading}>
                <BsBrushFill />
                <span>Toggle style: {this.props.style}</span>
            </button>
        );
    };
};

class ToggleLookButton extends React.Component {
    toggleLook = () => {
        if (this.props.look === "Cartoon") {
            this.props.setLook("Glossy");
        } else if (this.props.look === "Glossy") {
            this.props.setLook("Cartoon");
        } else {
            toast.error("Unknown look!", { autoClose: true });
            this.props.setLook("Cartoon");
        }
    };

    render() {
        return (
            <button className="square-button" onClick={this.toggleLook} disabled={this.props.isLoading}>
                <BsEyeFill />
                <span>Toggle look: {this.props.look}</span>
            </button>
        );
    };
};

class IncludeHydrogensButton extends React.Component {
    toggleIncludeHydrogens = () => {
        this.props.setIncludeHydrogens(!this.props.includeHydrogens);
    };

    render() {
        return (
            <button className="square-button" onClick={this.toggleIncludeHydrogens} disabled={this.props.isLoading}>
                <BsDropletFill />
                <span>Toggle include hydrogens: {this.props.includeHydrogens ? "True" : "False"}</span>
            </button>
        );
    };
};

class ToggleBackGroundButton extends React.Component {
    toggleBackground = () => {
        if (this.props.mode === "light") {
            this.props.setMode("dark");
        }  else {
            this.props.setMode("light");
        }
    };

    render() {
        return (
            <button className="square-button" onClick={this.toggleBackground} disabled={this.props.isLoading}>
                <BsCircleHalf />
                <span>Toggle background: {this.props.mode.charAt(0).toUpperCase() + this.props.mode.slice(1)}</span>
            </button>
        );
    };
};


class OpgenGithubIssuesButton extends React.Component {
    openNewTab = () => { window.open("https://github.com/moltools/CineMol/issues", "_blank"); };
  
    render() {
        return (
            <button className="square-button" onClick={this.openNewTab} disabled={this.props.isLoading}>
                <BsGithub />
                <span>Open GitHub issues</span>
            </button>
        );
    };
};

class CounterButton extends React.Component {  
    handleIncrement = () => {
        if (this.props.resolution < 100) {
            this.props.setResolution(this.props.resolution + 5);
        }
    };
  
    handleDecrement = () => {
        if (this.props.resolution >= 30) {
            this.props.setResolution(this.props.resolution - 5);
        }
    };
  
    render() {
        return (
            <div className="increment-container">
                <div className="increment-container-part-left">
                    <BsFillLightningFill />
                    <span className="increment-title">Resolution: {this.props.resolution}</span>
                </div>
                <div className="increment-container-part-right">
                    <button className="increment-button left" onClick={this.handleDecrement} disabled={this.props.isLoading}>{"<"}</button>
                    <button className="increment-button right" onClick={this.handleIncrement} disabled={this.props.isLoading}>{">"}</button>
                </div>
            </div>
        );
    };
};

class RotateXButton extends React.Component {
    handleIncrement = () => {
        if (this.props.rotationX < (2 * Math.PI - Math.PI / 12)) {
            this.props.setRotationX(Math.round((this.props.rotationX + Math.PI / 12) * 10) / 10);
        }
    };

    handleDecrement = () => {
        if (this.props.rotationX >= (0 + Math.PI / 12)) {
            this.props.setRotationX(Math.round((this.props.rotationX - Math.PI / 12) * 10) / 10);
        }
    };

    render() {
        return (
            <div className="increment-container">
                <div className="increment-container-part-left">
                    <BsGlobe2 />
                    <span className="increment-title">Rotation X: {this.props.rotationX} rad.</span>
                </div>
                <div className="increment-container-part-right">
                    <button className="increment-button left" onClick={this.handleDecrement} disabled={this.props.isLoading}>{"<"}</button>
                    <button className="increment-button right" onClick={this.handleIncrement} disabled={this.props.isLoading}>{">"}</button>
                </div>
            </div>
        );
    };
};

class RotateYButton extends React.Component {
    handleIncrement = () => {
        if (this.props.rotationY < (2 * Math.PI - Math.PI / 12)) {
            this.props.setRotationY(Math.round((this.props.rotationY + Math.PI / 12) * 10) / 10);
        }
    };

    handleDecrement = () => {
        if (this.props.rotationY >= (0 + Math.PI / 12)) {
            this.props.setRotationY(Math.round((this.props.rotationY - Math.PI / 12) * 10) / 10);
        }
    };

    render() {
        return (
            <div className="increment-container">
                <div className="increment-container-part-left">
                    <BsGlobe2 />
                    <span className="increment-title">Rotation Y: {this.props.rotationY} rad.</span>
                </div>
                <div className="increment-container-part-right">
                    <button className="increment-button left" onClick={this.handleDecrement} disabled={this.props.isLoading}>{"<"}</button>
                    <button className="increment-button right" onClick={this.handleIncrement} disabled={this.props.isLoading}>{">"}</button>
                </div>
            </div>
        );
    };
};

class RotateZButton extends React.Component {
    handleIncrement = () => {
        if (this.props.rotationZ < (2 * Math.PI - Math.PI / 12)) {
            this.props.setRotationZ(Math.round((this.props.rotationZ + Math.PI / 12) * 10) / 10);
        }
    };

    handleDecrement = () => {
        if (this.props.rotationZ >= (0 + Math.PI / 12)) {
            this.props.setRotationZ(Math.round((this.props.rotationZ - Math.PI / 12) * 10) / 10);
        }
    };

    render() {
        return (
            <div className="increment-container">
                <div className="increment-container-part-left">
                    <BsGlobe2 />
                    <span className="increment-title">Rotation X: {this.props.rotationZ} rad.</span>
                </div>
                <div className="increment-container-part-right">
                    <button className="increment-button left" onClick={this.handleDecrement} disabled={this.props.isLoading} >{"<"}</button>
                    <button className="increment-button right" onClick={this.handleIncrement} disabled={this.props.isLoading}>{">"}</button>
                </div>
            </div>
        );
    };
};


// =====================================================================================================================
// Component for Sidebar.
// =====================================================================================================================

const Sidebar = props => {    
    const [sdfString, setSdfString] = useState("");
    const [style, setStyle] = useState("Ball-and-stick");
    const [look, setLook] = useState("Cartoon");
    const [initialRender, setInitialRender] = useState(true);
    const [includeHydrogens, setIncludeHydrogens] = useState(false);
    const [resolution, setResolution] = useState(25);
    const [rotationX, setRotationX] = useState(0); // in radians
    const [rotationY, setRotationY] = useState(0); // in radians
    const [rotationZ, setRotationZ] = useState(0); // in radians

    const sidebarOpen = <BsChevronDoubleLeft className="cinemol-sidebar-toggle-icon" />;
    const sidebarClosed = <BsChevronDoubleRight className="cinemol-sidebar-toggle-icon" />;
    
    // Send SDF to API and retrieve SVG.
    const drawModel = () => {
        
        // set is loading to true
        props.setIsLoading(true);

        return new Promise((resolve, reject) => {
            fetch("/api/draw_model", {
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
                    "rotation_z": rotationZ
                })
            })
                .then(response => response.json())
                .then(json => {
                    resolve(json);
                    if (json.status === "success") {
                        props.setSvgString(json.payload.svg_string);
                    } else if (json.status === "warning") {
                        toast.warn(json.message, { autoClose: true });
                    } else if (json.status === "failure") {
                        toast.error(json.message, { autoClose: true });
                    } else {
                        toast.error("Unknown error!", { autoClose: true });
                    }
                    // set is loading to false
                    props.setIsLoading(false);
            })
            .catch(error => {
                const msg = "Could not parse SDF string!"
                toast.error(msg, { autoClose: true });
                console.log(msg)
                console.log(error);
                reject(error);
                // set is loading to false
                props.setIsLoading(false);
            });
        });
    };

    // useEffect to redraw the model when look or style changes, but skip the initial render.
    useEffect(() => {
        if (initialRender) {
            setInitialRender(false);
        } else {
            drawModel();
        }
    }, [look, style, includeHydrogens, resolution, rotationX, rotationY, rotationZ, sdfString]);

    return (
        <div className={"cinemol-sidebar" + (props.isOpen ? " open" : "")}>
            <button onClick={props.toggleSidebar} className="cinemol-sidebar-toggle">
                {props.isOpen ? sidebarOpen : sidebarClosed }
            </button>
            <div className="cinemol-sidebar-content">
                <div className="cinemol-sidebar-version">Version 0.1.0</div>
                <UploadButton isLoading={props.isLoading} setSdfString={setSdfString} />
                <DownloadSvgButton isLoading={props.isLoading} downloadSvgString={props.downloadSvgString} />
                <LoadExampleButton isLoading={props.isLoading} setSdfString={setSdfString} />
                <ToggleStyleButton isLoading={props.isLoading} style={style} setStyle={setStyle} />
                <ToggleLookButton isLoading={props.isLoading} look={look} setLook={setLook} />
                <IncludeHydrogensButton isLoading={props.isLoading} includeHydrogens={includeHydrogens} setIncludeHydrogens={setIncludeHydrogens} />
                <ToggleBackGroundButton isLoading={props.isLoading} mode={props.mode} setMode={props.setMode} />
                <OpgenGithubIssuesButton isLoading={props.isLoading} />
                <CounterButton isLoading={props.isLoading} resolution={resolution} setResolution={setResolution} />
                <RotateXButton isLoading={props.isLoading} rotationX={rotationX} setRotationX={setRotationX} />
                <RotateYButton isLoading={props.isLoading} rotationY={rotationY} setRotationY={setRotationY} />
                <RotateZButton isLoading={props.isLoading} rotationZ={rotationZ} setRotationZ={setRotationZ} />
            </div>
        </div>
    );
};

// =====================================================================================================================
// Render page.
// =====================================================================================================================

const CineMol = () => {
    const [svgString, setSvgString] = useState("");
    const [sidebarOpen, setSideBarOpen] = useState(true);
    const [mode, setMode] = useState("dark");
    const [isLoading, setIsLoading] = useState(false);

    const handleViewSidebar = () => { 
        setSideBarOpen(!sidebarOpen); 
    };

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

    return (
        <div className="cinemol-container">
            <Sidebar 
                isOpen={sidebarOpen} 
                toggleSidebar={handleViewSidebar} 
                svgString={svgString}
                setSvgString={setSvgString} 
                downloadSvgString={handleDownloadSvgString}
                mode={mode}
                setMode={setMode}
                setIsLoading={setIsLoading}
                isLoading={isLoading}
            />
            <div className={`cinemol-viewer-container ${mode}`}>
                <div className={`cinemol-viewer ${isLoading ? 'grayed-out' : ''}`} dangerouslySetInnerHTML={{ __html: svgString }} />
                {/* <div>
                    <img className="cinemol-viewer" src={`data:image/svg+xml;utf8,${encodeURIComponent(svgString)}`} />
                </div> */}
           </div>
        </div>
    );
};

export default CineMol;