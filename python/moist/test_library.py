from moist.library import get_api_version


def test_api_version_format() -> None:
    version = get_api_version()
    parts = version.split(".")
    assert len(parts) == 3
    assert all(part.isdigit() for part in parts)
