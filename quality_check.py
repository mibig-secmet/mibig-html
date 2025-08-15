#!/usr/bin/env python3
"""Check entry quality and set the quality level"""

from argparse import ArgumentParser, FileType

import orjson

from mibig.converters.shared.mibig import MibigEntry
from mibig.converters.shared.common import QualityLevel


def main():
    parser = ArgumentParser(description="Check entry quality and set the quality level")
    parser.add_argument("input_file", type=FileType("rb"), help="Input file in JSON format")
    parser.add_argument("output_file", type=FileType("wb"), help="Output file in JSON format")
    args = parser.parse_args()

    # Load the input data
    data = orjson.loads(args.input_file.read())

    # Process the data
    entry = MibigEntry.from_json(data)
    for level in (QualityLevel.HIGH,
                  QualityLevel.MEDIUM,
                  QualityLevel.LOW):
        errors = entry.validate(quality=level)
        if errors:
            continue
        entry.quality = level
        break


    # Write the output
    args.output_file.write(orjson.dumps(entry.to_json()))


if __name__ == "__main__":
    main()