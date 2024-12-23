import sys

from loguru import logger


def configure_logger():
    logger.remove()  # Remove the default logger
    logger.add(sys.stderr, level="INFO")


configure_logger()
