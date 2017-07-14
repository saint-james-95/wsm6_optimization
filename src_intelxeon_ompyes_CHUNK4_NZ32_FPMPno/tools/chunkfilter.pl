#!/usr/bin/perl -p
# insert compile-time chunk sizes
s/\b(do\s+\d*\s*\w+\s*=\s*)ii*ts\s*,\s*ii*te\b/$1 1,4/gi; s/\bii*[mt]s\s*:\s*ii*[mt]e\b/1:4/gi;
