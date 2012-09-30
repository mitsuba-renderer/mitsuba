#!/bin/bash
dependencies/bin/crash_report -S symbols64/ $@ | less
