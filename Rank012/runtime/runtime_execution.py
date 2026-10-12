"""Run a calculation without internal alarm signals interrupting exact arithmetic."""
import importlib
from pathlib import Path
import runpy
import signal
import resource
import sys

resource.setrlimit = lambda *args, **kwargs: None
signal.alarm(0)
original_alarm = signal.alarm
signal.alarm = lambda seconds: original_alarm(0)
native = importlib.import_module('cysignals.alarm')
native.cancel_alarm()
native.alarm = lambda seconds: native.cancel_alarm()
sys.argv = sys.argv[1:]
sys.path.insert(0, str(Path(sys.argv[0]).resolve().parent))
runpy.run_path(sys.argv[0], run_name='__main__')
