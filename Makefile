omit = --omit '*/external/*'
coverage = coverage run $(omit) --source mibig_html -m pytest
integration_flags = --override-ini=python_files=integration_*.py
integration_coverage = .coverage_integration
sanity_run = echo "sanity run missing"

unit:
	$(sanity_run)
	echo "simple reuse test missing"
	pytest --durations=3 mibig_html 

integration: clean
	python -m pytest --durations=3 mibig_html $(integration_flags)

clean:
	find . -name "*.h3?" -exec rm {} +
	find . -name '*.pyc' | xargs rm -f
	find . -name '__pycache__' | xargs rm -rf
	find . -name '*.dmnd' | grep -v test | xargs rm -f

squeakyclean: clean
	find . -name "*.tar.*" -exec rm {} +
	bash -c 'for d in $$(find . -maxdepth 2 -name "index.html"); do DIR=$$(dirname $$d); rm -r $$DIR; done'

cover: coverage

combined-coverage: coverage
	COVERAGE_FILE=$(integration_coverage) $(coverage) $(integration_flags)
	coverage combine -a $(integration_coverage)
	coverage html -d cover
	coverage report

coverage:
	$(sanity_run)
	rm -rf cover .coverage $(integration_coverage)
	coverage run $(omit),'*integration_*.py' --source mibig_html -m pytest mibig_html
	coverage html -d cover
	coverage report

.PHONY:	unit integration clean squeakyclean cover coverage combined-coverage
