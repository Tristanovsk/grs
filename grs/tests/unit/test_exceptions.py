import pytest

from grs.grs_exceptions import GRSException, GRS_IO_Exception


class TestGRSException:
    """Tests for the GRSException class."""

    def test_init_stores_message(self):
        exception = GRSException("Something went wrong")

        assert exception.msg == "Something went wrong"

    def test_init_stores_cause(self):
        cause = ValueError("Invalid value")

        exception = GRSException("Something went wrong", cause)

        assert exception.cause is cause

    def test_init_default_cause_is_none(self):
        exception = GRSException("Something went wrong")

        assert exception.cause is None

    def test_str_without_cause(self):
        exception = GRSException("Something went wrong")

        assert str(exception) == "'Something went wrong'\nNone"

    def test_str_with_cause(self):
        cause = ValueError("Invalid value")
        exception = GRSException("Something went wrong", cause)

        assert str(exception) == (
            "'Something went wrong'\n"
            "ValueError('Invalid value')"
        )

    def test_is_exception_subclass(self):
        exception = GRSException("Something went wrong")

        assert isinstance(exception, Exception)

    def test_can_be_raised_and_caught(self):
        with pytest.raises(GRSException) as exc_info:
            raise GRSException("Something went wrong")

        assert exc_info.value.msg == "Something went wrong"


class TestGRSIOException:
    """Tests for the GRS_IO_Exception class."""

    def test_init_formats_message_with_file_path(self):
        exception = GRS_IO_Exception("/tmp/data.txt", OSError("Permission denied"))

        assert exception.msg == (
            "I/O Exception, error while working with /tmp/data.txt"
        )

    def test_init_stores_original_cause(self):
        cause = OSError("Permission denied")

        exception = GRS_IO_Exception("/tmp/data.txt", cause)

        assert exception.cause is cause

    def test_init_inherits_from_grs_exception(self):
        exception = GRS_IO_Exception("/tmp/data.txt", OSError("File not found"))

        assert isinstance(exception, GRSException)
        assert isinstance(exception, Exception)

    def test_str_includes_message_and_cause(self):
        cause = OSError("Permission denied")
        exception = GRS_IO_Exception("/tmp/data.txt", cause)

        assert str(exception) == (
            "'I/O Exception, error while working with /tmp/data.txt'\n"
            "OSError('Permission denied')"
        )

    def test_can_be_raised_and_caught(self):
        cause = OSError("File not found")

        with pytest.raises(GRS_IO_Exception) as exc_info:
            raise GRS_IO_Exception("/tmp/data.txt", cause)

        assert exc_info.value.cause is cause

    def test_can_be_caught_as_grs_exception(self):
        exception = GRS_IO_Exception(
            "/tmp/data.txt",
            OSError("File not found"),
        )

        with pytest.raises(GRSException):
            raise exception


    @pytest.mark.parametrize(
        "file_path",
        [
            "data.txt",
            "/tmp/data.txt",
            "/path/with spaces/data.txt",
            "",
        ],
    )
    def test_init_supports_different_file_paths(self, file_path):
        exception = GRS_IO_Exception(file_path, OSError("I/O error"))

        assert exception.msg == (
            f"I/O Exception, error while working with {file_path}"
        )
