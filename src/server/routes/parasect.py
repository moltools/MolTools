from flask import Blueprint, Response, request 

from .common import Status, ResponseData

blueprint_predict_parasect = Blueprint("predict_parasect", __name__)
@blueprint_predict_parasect.route("/api/predict_parasect", methods=["POST"])
def predict_parasect() -> Response:
    """
    API endpoint for making a prediction with PARAS.

    :return: Response data.
    :rtype: ResponseData
    """
    # Parse request data.
    data = request.get_json()

    print(data)
    
    # msg = "Backend for PARASECT not fully implemented yet!"
    # return ResponseData(Status.Warning, message=msg).to_dict()

    payload = {"predictions": {"x": 0.1, "y": 0.13, "z": 0.3}}
    print(payload)
    return ResponseData(Status.Success, message="Success!", payload=payload).to_dict()
