from flask import Blueprint, Response
import neo4j

from .common import Status, ResponseData

blueprint_ping_server = Blueprint("ping_server", __name__)
@blueprint_ping_server.route("/api/ping_server", methods=["GET"])   
def ping_server() -> Response:
    return ResponseData(Status.Success).to_dict()

blueprint_ping_database = Blueprint("ping_database", __name__)
@blueprint_ping_database.route("/api/ping_database", methods=["GET"])
def ping_database() -> Response:
    try:
        # driver = neo4j.GraphDatabase.driver("bolt://database:7687")
        driver = neo4j.GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "password"))
        with driver.session() as session:
            result = session.run("MATCH (n) RETURN count(n) AS count")
            count = result.single()["count"]
        driver.close()
        message = f"Database is up and running with {count} nodes!"
        return ResponseData(Status.Success, message=message).to_dict()
    except Exception as e:
        return ResponseData(Status.Failure, message=f"{e}").to_dict()