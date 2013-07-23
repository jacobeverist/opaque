import pstats

if __name__ == "__main__":
	
	p = pstats.Stats('fooprof')
	#p.strip_dirs().sort_stats(-1).print_stats()
	#p.strip_dirs().sort_stats("cumulative").print_stats()
	p.strip_dirs().sort_stats("time").print_stats()

