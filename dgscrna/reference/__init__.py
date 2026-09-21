"""Public interface to the packaged R reference workflow."""


def run_reference(**kwargs):
    """Run the workflow; see reference.runner.run_reference for arguments."""
    from .runner import run_reference as run
    return run(**kwargs)


__all__ = ["run_reference"]
