# Kept so `python SAGAconf_parser.py ...` works as before. New interface: `sagaconf parse` (sagaconf/cli.py).
from sagaconf.cli import legacy_parse_main

if __name__ == "__main__":
    legacy_parse_main()
