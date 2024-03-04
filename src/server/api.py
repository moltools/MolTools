from flask import Flask, Response

from routes.check_status import blueprint_ping_server, blueprint_ping_database
from routes.pikachu import blueprint_draw_smiles_with_pikachu
from routes.chemistry import blueprint_draw_smiles
from routes.biosynfoni import blueprint_predict_biosynthetic_class
from routes.cinemol import blueprint_draw_model, blueprint_fetch_cinemol_version
from routes.retromol import blueprint_parse_retromol, blueprint_embed_retromol, blueprint_find_matches, blueprint_query_database 

app = Flask(__name__)
app.register_blueprint(blueprint_ping_server)
app.register_blueprint(blueprint_ping_database)
app.register_blueprint(blueprint_draw_smiles_with_pikachu)
app.register_blueprint(blueprint_draw_smiles)
app.register_blueprint(blueprint_predict_biosynthetic_class)
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

def main() -> None:
    """
    Run the app locally for development.
    """
    # Main is only run when the file is run directly during local development.
    # When running in Docker, the app is run directly from the Dockerfile.
    app.run(host="localhost", port=4000, debug=True)

if __name__ == "__main__":
    main()