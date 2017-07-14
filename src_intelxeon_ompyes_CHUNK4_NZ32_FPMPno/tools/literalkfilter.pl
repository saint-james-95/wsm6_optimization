#!/usr/bin/perl -p
# insert compile-time vertical sizes
 s/\b(do\s+\d*\s*\w+\s*=\s*)kk*ts\s*,\s*kk*te\b/$1 1,32/gi;  s/\b(do\s+\d*\s*\w+\s*=\s*)kk*te\s*,\s*kk*ts\b/$1 32,1/gi;  s/\b(do\s+\d*\s*\w+\s*=\s*)kk*ts\s*,\s*kk*ts\b/$1 1,1/gi;  s/\b(do\s+\d*\s*\w+\s*=\s*)kk*ts\s*,\s*kk*te\b/$1 32,32/gi; s/\bkk*ts\s*:\s*kk*te\b/1:32/gi;s/\bkk*ms\s*:\s*kk*me\b/1:33/gi;
