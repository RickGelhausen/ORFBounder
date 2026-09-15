"""Distribution metadata contracts for installed commands and helper modules."""

from pathlib import Path
import tomllib


REPOSITORY = Path(__file__).resolve().parents[1]


def test_wheel_configuration_includes_commands_and_documented_helpers():
    configuration = tomllib.loads(
        (REPOSITORY / "pyproject.toml").read_text(encoding="utf-8")
    )
    readme = (REPOSITORY / "README.md").read_text(encoding="utf-8")

    assert configuration["project"]["scripts"] == {
        "orfbounder": "orfbounder:main",
        "orfbounder-batch": "orfbounder_batch:main",
    }
    assert set(configuration["tool"]["setuptools"]["packages"]) >= {"lib", "helpers"}
    for module in (
        "final_to_gff.py",
        "merge_tables.py",
        "recover_total_counts.py",
        "toy_example_generation.py",
    ):
        assert (REPOSITORY / "helpers" / module).is_file()
    for invocation in (
        "python -m helpers.toy_example_generation",
        "python -m helpers.recover_total_counts",
        "python -m helpers.final_to_gff",
        "python -m lib.merging",
    ):
        assert invocation in readme
    assert "helpers.merge_tables" in readme
    assert "compatibility entry point" in readme


def test_container_build_context_includes_packaged_helpers():
    dockerfile = (REPOSITORY / "Dockerfile").read_text(encoding="utf-8")
    dockerignore = (REPOSITORY / ".dockerignore").read_text(encoding="utf-8")

    assert "COPY helpers ./helpers" in dockerfile
    assert "!helpers/" in dockerignore
    assert "!helpers/*.py" in dockerignore
