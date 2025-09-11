from gene_order import ALPHA_ORDER, LOCUS_ORDER, PSEUDO_GENES

done = []
print('\n\n\nLOCUS_ORDER')
for x in LOCUS_ORDER:
	g = x.split('*')[0]
	if g not in done:
		print('"%s",' % g)
		done.append(g)
		
print('ALPHA_ORDER')
done = []
for x in ALPHA_ORDER:
	g = x.split('*')[0]
	if g not in done:
		print('"%s",' % g)
		done.append(g)
		
done = []
print('\n\n\nPSEUDO_GENES')
for x in PSEUDO_GENES:
	g = x.split('*')[0]
	if g not in done:
		print('"%s",' % g)
		done.append(g)
		
		