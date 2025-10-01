import sys

from loguru import logger


def configure_logger(debug: bool = False):
    logger.remove()  # Remove the default logger
    if debug:
        logger.add(
            sys.stderr,
            level="DEBUG",
            format="<green>{time:YYYY-MM-DD HH:mm:ss}</> | <cyan>{file}:{function}:{line}</>  | <lvl>{level:8s}</> | <lvl>{message}</>",
        )
    else:
        logger.add(
            sys.stderr,
            level="INFO",
            format="<green>{time:YYYY-MM-DD HH:mm:ss}</> | <lvl>{level:8s}</> | <lvl>{message}</>",
        )
