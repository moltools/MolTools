from flask import Flask, Response
import neo4j 

from routes.common import Status, ResponseData
from routes.cinemol import blueprint_draw_model, blueprint_fetch_cinemol_version
from routes.retromol import (
    blueprint_parse_retromol, 
    blueprint_embed_retromol, 
    blueprint_find_matches, 
    blueprint_query_database 
)

app = Flask(__name__)
app.register_blueprint(blueprint_draw_model)
app.register_blueprint(blueprint_fetch_cinemol_version)
app.register_blueprint(blueprint_parse_retromol)
app.register_blueprint(blueprint_embed_retromol)
app.register_blueprint(blueprint_find_matches)
app.register_blueprint(blueprint_query_database)

@app.errorhandler(404)
def not_found(e) -> Response:
    """
    Return a custom 404 error.
    """
    return app.send_static_file("index.html")

@app.route("/")
def index() -> Response:
    """
    Return the index page.
    """
    return app.send_static_file("index.html")

@app.route("/api/ping_server", methods=["GET"])   
def ping_server() -> Response:
    """
    API endpoint for checking the server status.
    """
    return ResponseData(Status.Success).to_dict()

@app.route("/api/ping_database", methods=["GET"])
def ping_database() -> Response:
    """
    API endpoint for checking the database status.
    """
    try:
        # driver = neo4j.GraphDatabase.driver("bolt://database:7687") # Change for Docker 
        driver = neo4j.GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))

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
    """
    Run the app locally for development.
    """
    # Main is only run when the file is run directly during local development.
    # When running in Docker, the app is run directly from the Dockerfile.
    app.run(host="localhost", port=4000, debug=True)

if __name__ == "__main__":
    main()