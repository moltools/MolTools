from flask import Blueprint, Response

from .common import Status, ResponseData

blueprint_ping_server = Blueprint("ping_server", __name__)
@blueprint_ping_server.route("/api/ping_server", methods=["GET"])   
def ping_server() -> Response:
    """
    Return a simple health check.

    :return: A simple health check.
    :rtype: Response
    """
    return ResponseData(Status.Success).to_dict()