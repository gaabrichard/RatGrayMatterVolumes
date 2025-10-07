#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re

def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('folder',
                   help='path containing files to be renamed')
    p.add_argument('pattern', help='pattern we want to replace')
    p.add_argument('substitution', help='string that will replace the pattern')
    p.add_argument('-p', '--prefix',
                   help='prefix for filename')
    p.add_argument('-s', '--suffix',
                   help='suffix for filename')
    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()
    
    for filename in sorted(os.scandir(args.folder), key=lambda e: e.name):
        if filename.is_file():
           old_full_name = filename.path
           name = filename.name 
           name = re.sub(args.pattern, args.substitution, name)
           new_full_name = os.path.join(args.folder, name)
           os.replace(old_full_name, new_full_name)

if __name__ == "__main__":
    main()

