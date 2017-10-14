file(REMOVE_RECURSE
  "libDNAPregel-iniparser.pdb"
  "libDNAPregel-iniparser.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/DNAPregel-iniparser.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
