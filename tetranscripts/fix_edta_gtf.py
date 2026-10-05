import sys

#convert gff3 to gtf
def fix_edta_to_gtf(gff_in, gtf_out):
  with open(gff_in, 'r') as infile, open(gtf_out, 'w') as outfile:
    count = 0
    for line in infile:
      if line.startswith('#') or not line.strip():
        continue

      parts = line.strip().split('\t')
      if len(parts) < 9:
        continue

      chrom = parts[0]
      source = 'EDTA'
      feature = 'exon'
      start = parts[3]
      end = parts[4]
      score = parts[5] if parts[5] != '.' else '0'
      strand = parts[6] if parts[6] in ['+', '-'] else '.'
      frame = '.'

      attr_list = parts[8].split(';')
      attrs = {}
      for item in attr_list:
        if '=' in item:
          k, v = item.split('=', 1)
          attrs[k] = v

      # Clear names (remove ":" to avoid conflict TEtranscripts)
      raw_name = (
          attrs.get('Name', f'TE_{count}').replace(':', '_').replace(';', '_')
      )
      classification = attrs.get('Classification', 'Unknown').replace(':', '_')
      
      # separate class and family (ex: LTR/Copia -> class_id=LTR, family_id=Copia)
      if '/' in classification:
        te_class, te_family = classification.split('/', 1)
      else:
        te_class = classification
        te_family = classification

      count += 1
      # Define unique IDs for each lócus to avoid mismm
      gene_id = f'{raw_name}_locus{count}'
      transcript_id = f'{raw_name}_tx{count}'

      gtf_attributes = (
          f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; family_id'
          f' "{te_family}"; class_id "{te_class}";'
      )

      outfile.write(
          f'{chrom}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{gtf_attributes}\n'
      )


if __name__ == '__main__':
  if len(sys.argv) != 3:
    print('Uso: python fix_edta_gtf.py <input.gff3> <output.gtf>')
    sys.exit(1)
  fix_edta_to_gtf(sys.argv[1], sys.argv[2])
#
#end