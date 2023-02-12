git:
	git add .
	git commit -m "$m"
	git push -u origin main

clean:
	rm -rf work/
	rm -rf output/
	rm -f slurm*.out
	rm -f trace*.txt
	rm -f .nextflow.log*
	clear

kill:
	squeue -u blazquL | grep $n | awk '{print $1}' | xargs -n 1 scancel