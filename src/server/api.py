from flask import Flask, Response

from routes.ping_server import blueprint_ping_server

app = Flask(__name__)
app.register_blueprint(blueprint_ping_server)

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