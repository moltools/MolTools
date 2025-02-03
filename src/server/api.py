#!/usr/bin/env python
"""This module contains the main API server for the application."""
from flask import Flask, Response

from routes.common import success
from routes.cinemol import blueprint_draw_model, blueprint_fetch_cinemol_version
from routes.biosynfoni import blueprint_predict_biosynthetic_class, blueprint_biosynfoni_version

app = Flask(__name__)
app.register_blueprint(blueprint_draw_model)
app.register_blueprint(blueprint_fetch_cinemol_version)
app.register_blueprint(blueprint_predict_biosynthetic_class)
app.register_blueprint(blueprint_biosynfoni_version)


@app.route("/")
def index() -> Response:
    """Return the index page.

    :return: The index page.
    :rtype: Response
    """
    return app.send_static_file("index.html")


@app.route("/api/fetch_server_status", methods=["GET"])
def fetch_server_status() -> Response:
    """API endpoint for checking the server status.

    :return: The response to return.
    :rtype: Response
    """
    return success("Server is up and running!")


def main() -> None:
    """Run the app locally for development."""
    # Main is only run when the file is run directly during local development.
    # When running in Docker, the app is run directly from the Dockerfile.
    app.run(host="localhost", port=4000, debug=True)


if __name__ == "__main__":
    main()
