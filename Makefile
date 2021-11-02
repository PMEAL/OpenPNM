.SILENT: clean nbtest test
.PHONY: clean nbtest test

clean:
	echo "Cleaning up the repo..."
	python bin/clean

nbtest:
	echo "Running example notebooks..."
	pytest examples --nbval-lax \
      --ignore="examples/paper_recreations/Blunt et al. (2013)" \
	  --ignore="examples/paper_recreations/Wu et al. (2010)"

test:
	echo "Running unit tests..."
	pytest --ignore=scripts

