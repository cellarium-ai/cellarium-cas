import asyncio
import typing as t

from mockito import when


def async_context_manager_decorator(mocked_object: object) -> object:
    """
    Takes in a mock object and decorates it so that it can be used as an async context manager.
    E.g. as: `async with mock_object as obj: ...`

    The `__aenter__` method simpy returns the mock object, and the `__aexit__` method returns `None`.

    :param mocked_object: The mock object to decorate. Note: Mockito doesn't expose their mock object type.

    :return: The decorated mock object.
    """

    async def __aenter__(*_):
        return mocked_object

    async def __aexit__(*_):
        return None

    when(mocked_object).__aenter__().thenAnswer(__aenter__)
    when(mocked_object).__aexit__(...).thenAnswer(__aexit__)

    return mocked_object


def async_return(result: t.Any) -> asyncio.Future:
    """
    Return a result from an async function. Useful for mocking client calls.

    :param result: The result to return.
    :return: An asyncio Future that will return the result.
    """
    f = asyncio.Future()
    f.set_result(result)
    return f
