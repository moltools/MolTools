from flask import Blueprint, Response, request 

from .common import Status, ResponseData

blueprint_predict_paras = Blueprint("predict_paras", __name__)
@blueprint_predict_paras.route("/api/predict_paras", methods=["POST"])
def predict_paras() -> Response:
    """
    API endpoint for making a prediction with PARAS.

    :return: Response data.
    :rtype: ResponseData
    """
    # Parse request data.
    data = request.get_json()

    print(data)
    
    msg = "Backend for PARAS not fully implemented yet!"
    return ResponseData(Status.Warning, message=msg).to_dict()
