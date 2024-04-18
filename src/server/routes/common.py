"""This module contains common classes and functions used by the API endpoints."""
import typing as ty
from enum import Enum, auto

class Status(Enum):
    """The status of a response.

    :cvar Success: The response was successful.
    :cvar Warning: The response was successful, but with a warning.
    :cvar Failure: The response failed.
    """
    Success = auto()
    Warning = auto()
    Failure = auto()

    def __str__(self) -> str:
        return self.name.lower()

class ResponseData:
    """A class for creating a response object.
    
    :param status: The status of the response.
    :type status: Status
    :param payload: The payload of the response.
    :type payload: dict
    :param message: The message of the response.
    """
    def __init__(
        self,
        status: Status,
        payload: ty.Optional[dict] = None,
        message: ty.Optional[str] = None
    ) -> None:
        self.status = status
        self.payload = payload if payload is not None else dict()
        self.message = message if message is not None else "No message provided!"

    def to_dict(self) -> ty.Dict[str, ty.Any]:
        """Return the response as a dictionary.
        
        :return: The response as a dictionary.
        :rtype: ty.Dict[str, ty.Any]
        """
        return dict(
            status=str(self.status),
            payload=self.payload,
            message=self.message
        )
