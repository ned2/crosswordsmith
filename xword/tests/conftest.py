from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"
GOLDEN = Path(__file__).parent / "golden"


@pytest.fixture
def fixtures() -> Path:
    return FIXTURES


@pytest.fixture
def golden() -> Path:
    return GOLDEN
