default:
	@echo "Try: make html iter"

add:
	git add -f main_files/*

html:
	R --no-restore --no-save --verbose -e "rmarkdown::render('main.Rmd')"

iter: add
	git commit . -m "iteration" && git push

meminfo:
	ps -eo pmem,pcpu,pid,args | tail -n +2 | sort -rnk 1 | head

