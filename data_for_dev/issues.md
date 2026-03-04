 The issue is that spatialcore 0.5.1 on PyPI is missing spatialcore/spatial/r_functions.R. The Python code at          
  spatialcore/spatial/domains.py:79 references it:                                            
                                                                                                                        
  R_FUNCTIONS_PATH = Path(__file__).parent / "r_functions.R"                                                            
                                                                                                                        
  And line 553 raises FileNotFoundError when it's not there. The R file isn't included in the package distribution —    
  likely missing from package_data or MANIFEST.in in the build config.
                                                                                                                        
  To fix: Include *.R files in the PyPI package. Depending on the build system:

  - pyproject.toml (hatchling/setuptools): add [tool.setuptools.package-data] with spatialcore = ["**/*.R"]
  - setup.py: add package_data={"spatialcore": ["**/*.R"]}
  - MANIFEST.in: add recursive-include spatialcore *.R

  Then cut a new release.