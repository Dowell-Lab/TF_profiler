dinucleotide
[outdir + '/generated_sequences/' + sample + '_position1_dinucleotide_probabilities.tsv',
outdir + '/generated_sequences/' + sample +'_conditional_probabilites_given'+ nuc + '.tsv',



mononucleotide
[outdir + '/generated_sequences/' + sample +'_mononucleotide_probabilites.tsv',


all
outdir + '/generated_sequences/' + str(sample) + '_' + seq_type + '.fa',
outdir + '/annotations/' + str(sample) + '_' + seq_type + '_centered.bed',
outdir + spef_dir + sample + '_' + seq_type + '_window.bed'

