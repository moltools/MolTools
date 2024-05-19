#!/usr/bin/env python
"""This module contains the main API server for the application."""
from flask import Flask, Response
import neo4j

from routes.common import Status, ResponseData
from routes.chem import blueprint_smiles_to_svg
from routes.cinemol import blueprint_draw_model, blueprint_fetch_cinemol_version
from routes.retromol import (
    NEO4J_URI,
    NEO4J_USER,
    NEO4J_PASSWORD,
    blueprint_submit_for_review,
    blueprint_bioactivity_labels,
    blueprint_producing_organisms,
    blueprint_parse_smiles,
    blueprint_parse_proto_cluster,
    blueprint_match_database
)

app = Flask(__name__)
app.register_blueprint(blueprint_smiles_to_svg)
app.register_blueprint(blueprint_draw_model)
app.register_blueprint(blueprint_fetch_cinemol_version)
app.register_blueprint(blueprint_submit_for_review)
app.register_blueprint(blueprint_bioactivity_labels)
app.register_blueprint(blueprint_producing_organisms)
app.register_blueprint(blueprint_parse_smiles)
app.register_blueprint(blueprint_parse_proto_cluster)
app.register_blueprint(blueprint_match_database)

@app.errorhandler(404)
def not_found(e) -> Response:
    """Return a custom 404 error.

    :param e: The error message.
    :type e: Exception
    :return: The response to return.
    :rtype: Response
    """
    print(e)
    return app.send_static_file("index.html")

@app.route("/")
def index() -> Response:
    """Return the index page.

    :return: The index page.
    :rtype: Response
    """
    return app.send_static_file("index.html")

@app.route("/api/ping_server", methods=["GET"])
def ping_server() -> Response:
    """API endpoint for checking the server status.

    :return: The response to return.
    :rtype: Response
    """
    return ResponseData(Status.Success).to_dict()

@app.route("/api/ping_database", methods=["GET"])
def ping_database() -> Response:
    """API endpoint for checking the database status.

    :return: The response to return.
    :rtype: Response
    """
    try:
        if NEO4J_USER and NEO4J_PASSWORD:
            driver = neo4j.GraphDatabase.driver(
                NEO4J_URI, 
                auth=(NEO4J_USER, NEO4J_PASSWORD)
            )
        else:
            driver = neo4j.GraphDatabase.driver(NEO4J_URI)

        with driver.session() as session:
            result = session.run("MATCH (n) RETURN count(n) AS count")
            count = result.single()["count"]

        driver.close()

        msg = f"Database is up and running with {count} nodes!"
        return ResponseData(Status.Success, message=msg).to_dict()

    except Exception as e:
        msg = f"No database connection: {e}"
        return ResponseData(Status.Failure, message=msg).to_dict()

def main() -> None:
    """Run the app locally for development."""
    # Main is only run when the file is run directly during local development.
    # When running in Docker, the app is run directly from the Dockerfile.
    app.run(host="localhost", port=4000, debug=True)

if __name__ == "__main__":
    main()
