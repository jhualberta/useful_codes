source /path/to/rat/env.sh
ratdb -s "postgres://snoplus:pwd@pgsql.snopl.us:5400/ratdb" list_all_objs_run 109137
ratdb -s "postgres://snoplus:pwd@pgsql.snopl.us:5400/ratdb" dump_table NOISE_RUN_LEVEL 109137
ratdb -s "postgres://snoplus:pwd@pgsql.snopl.us:5400/ratdb" dump_table CALIB_COMMON_RUN_LEVEL[MANIP] 109137

bunzip2 NOISE_RUN_LEVEL_100098-100098.ratdb.bz2
more NOISE_RUN_LEVEL_100098-100098.ratdb
