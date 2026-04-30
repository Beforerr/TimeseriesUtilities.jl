test name:
    julia-client --project=test -e 'using TestItemRunner; TestItemRunner.run_tests(pwd(); filter = ti -> ti.name == "{{name}}")'
