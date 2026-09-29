# Kept so `python SAGAconf.py ...` works as before. New interface: `sagaconf run` (sagaconf/cli.py).
from sagaconf.cli import legacy_run_main

if __name__ == "__main__":
    legacy_run_main()
