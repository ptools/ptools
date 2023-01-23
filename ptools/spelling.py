"""ptools.spelling - Miscelleanous spelling functions."""


def pluralize(word: str) -> str:
    """Given a word, returns its plural form."""

    def format_output(word: str) -> str:
        output = word
        if tokens:
            output = "_".join(tokens) + "_" + word
        return output

    irregular = {
        "child": "children",
        "goose": "geese",
        "man": "men",
        "woman": "women",
        "tooth": "teeth",
        "foot": "feet",
        "mouse": "mice",
        "person": "people",
        "index": "indices",  # "indexes" actually is correct as well
    }

    tokens = []

    if "_" in word:
        tokens = word.split("_")
        word = tokens.pop(-1)  # makes plural form only on last word

    if len(word) == 1:
        return word + "s"

    if word in irregular:
        return format_output(irregular[word])

    if word in ("sheep", "series", "species", "data", "deer", "coordinates"):
        return format_output(word)

    if word.endswith("us"):
        return format_output(word[:-2] + "i")

    if word.endswith("is"):
        return format_output(word[:-2] + "es")

    if word.endswith("on"):
        return format_output(word[:-2] + "a")

    for suffix in ("s", "ss", "sh", "ch", "x", "z"):
        if word.endswith(suffix):
            return format_output(word + "es")

    for suffix in ("f", "fe"):
        exceptions: tuple[str, ...] = ("roof", "belief", "chef", "chief")
        if word.endswith(suffix) and word not in exceptions:
            return format_output(word[: -len(suffix)] + "ves")

    if word.endswith("y") and word[-2] not in "aeiouy":
        return format_output(word[:-1] + "ies")

    if word.endswith("o"):
        exceptions = ("photo", "piano", "halo")
        if word not in exceptions:
            return format_output(word + "es")

    if word.endswith("ex"):
        return format_output(word[:-2] + "ices")

    word = word + "s"

    return format_output(word)
