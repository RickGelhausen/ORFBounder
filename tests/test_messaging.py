"""Terminal-message behavior for interactive users and captured batch logs."""

from io import StringIO

from lib import messaging


class _TerminalBuffer(StringIO):
    def isatty(self):
        return True


def test_noninteractive_output_has_no_ansi_sequences(capsys):
    messaging.message("reading")
    messaging.warning("check input")

    captured = capsys.readouterr()
    assert captured.out == "reading\n"
    assert captured.err == "check input\n"
    assert "\033[" not in captured.out + captured.err


def test_interactive_output_retains_color(monkeypatch):
    output = _TerminalBuffer()
    monkeypatch.setattr(messaging.sys, "stdout", output)
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.delenv("TERM", raising=False)

    messaging.success("done")

    assert output.getvalue() == f"{messaging.GREEN}done{messaging.ENDC}\n"


def test_no_color_environment_disables_interactive_color(monkeypatch):
    errors = _TerminalBuffer()
    monkeypatch.setattr(messaging.sys, "stderr", errors)
    monkeypatch.setenv("NO_COLOR", "")

    messaging.error("invalid")

    assert errors.getvalue() == "invalid\n"


def test_error_list_is_plain_in_noninteractive_logs(capsys):
    messaging.error_list(
        "Missing ", "Update the input.", "samples: ",
        ["TIS-A-1", "TIS-A-2"], ["TIS-A-1"],
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == (
        "Missing samples:  TIS-A-1\n"
        "          TIS-A-2\n"
        "Update the input.\n"
    )
    assert "\033[" not in captured.err
