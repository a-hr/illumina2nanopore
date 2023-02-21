clean:
	@rm -rf work/
	@rm -rf output/
	@rm -f slurm*.out
	@rm -f trace*.txt
	@rm -f .nextflow.log*

kill:
	squeue -u blazquL | grep $n | awk '{print $$1}' | xargs -n 1 scancel

pull:
	echo "Pulling containers ..."
	@mkdir -p containers
	@for container in `grep -oP "(?<=container = ').*(?=')" conf/cluster.config`; do \
		echo "Pulling $$container ..."; \
		containerName=`echo $$container | sed 's/\//-/g' | sed 's/:/-/g'`; \
		singularity -s pull --name $$containerName.img --dir containers/ docker://$$container; \
	done
	@echo "Done!"