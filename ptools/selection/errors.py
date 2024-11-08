"""Implement exceptions for the selection parser."""


class MisformattedExpressionError(Exception):
    message_fmt = "Misformatted expression at: '{}'"

    def __init__(self, token: str):
        self.token = token
        super().__init__(self.message_fmt.format(self.token))


class UnknownTokenError(MisformattedExpressionError):
    message_fmt = "Unknown token: '{}'"


class UnexpectedTokenError(Exception):
    """Error raised when a token was expected and found another one."""

    def __init__(self, expected: str, actual: str):
        self.expected = expected
        self.actual = actual
        message = f"Expected '{self.expected}', got '{self.actual}'"
        super().__init__(message)


class UnexpectedTrailingTokenError(Exception):
    def __init__(self, token: str):
        self.token = token
        message = f"Unexpected token at the end of the selection string: '{self.token}'"
        super().__init__(message)
