project(
  'arima_htslib', 'c',
  version : '0.1',
  default_options : [
    'buildtype=release', 'warning_level=3',
    'c_std=c99', 'default_library=static',
  ],
)

htslib = dependency('htslib', static: true)

filter_five_end = executable(
  'filter_five_end', 'src/filter_five_end/filter_five_end.c',
  install: true,
  dependencies: htslib,
)

get_stats = executable(
  'get_stats', 'src/get_stats/get_stats.c',
  install: true,
  dependencies: htslib,
)

two_read_bam_combiner = executable(
  'two_read_bam_combiner', 'src/two_read_bam_combiner/two_read_bam_combiner.c',
  install: true,
  dependencies: htslib,
)
