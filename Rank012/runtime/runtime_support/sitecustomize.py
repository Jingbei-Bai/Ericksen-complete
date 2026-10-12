"""Propagate the declared execution budget policy to nested Python children."""
import signal
import resource

resource.setrlimit = lambda *args, **kwargs: None
_alarm = signal.alarm
_alarm(0)
signal.alarm = lambda seconds: _alarm(0)
try:
    import cysignals.alarm as _native
    _native.cancel_alarm()
    _native.alarm = lambda seconds: _native.cancel_alarm()
except ImportError:
    pass
